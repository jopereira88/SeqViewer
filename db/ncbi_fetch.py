#!/usr/bin/env python3
"""
Resumable NCBI sequence fetch utility.

Downloads FASTA sequences from NCBI for accessions that are missing from the
local FASTA file and writes them to an output FASTA file, conforming headers
to the most frequent header format of the original file.

The download can be stopped at any time (Ctrl-C) and resumed later: progress
is tracked in a JSON file and completed sequences are appended to the output
FASTA, so re-running the command skips what is already done.

Standalone usage (recommended for large downloads)
--------------------------------------------------
    python db/ncbi_fetch.py \\
        --db data/sequence_database.sqlite \\
        --output data/fetched_sequences.fasta \\
        --progress data/fetch_progress.json \\
        --email you@example.org \\
        [--api-key XXXXXXXX] \\
        [--batch-size 200]

    python db/ncbi_fetch.py \\
        --accessions-file missing_accessions.txt \\
        --output data/fetched_sequences.fasta \\
        --progress data/fetch_progress.json \\
        --email you@example.org

Accession sources (at least one required)
-----------------------------------------
    --db
        A SQLite database. The missing set is computed as
        (metadata UNION cluster_composition) MINUS sequences.
    --accessions-file
        A text file with one accession per line.
    --accessions
        Inline accessions separated by commas, semicolons or whitespace.

The fetched FASTA file can later be concatenated with the original FASTA and
passed to db/import_data.py as --fasta, or the ETL can be run with
--fetch-missing which reuses this module.
"""

import argparse
import json
import sqlite3
import time
from collections import Counter
from io import StringIO
from pathlib import Path

import os
import sys

# Allow running this module as a script (python db/ncbi_fetch.py) while
# importing the top-level 'utils' package.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import urllib.error

from Bio import Entrez, SeqIO

from utils.sequence_analysis import count_ambiguous_bases, degenerate_breakdown_json


# ============================================================
# HEADER FORMAT DETECTION / CONFORMING
# ============================================================

DEFAULT_DELIMITER = " |"


def _split_header(header):
    """
    Split one FASTA header into (accession, delimiter, description).

    Handles the common 'ACCESSION |description', 'ACCESSION|description',
    'ACCESSION | description' and plain 'ACCESSION description' formats.
    """
    header = str(header).lstrip(">").strip()

    if not header:
        return None, None, None

    pipe_index = header.find("|")

    if pipe_index != -1:
        accession = header[:pipe_index].rstrip()
        delimiter = header[len(accession):pipe_index + 1]
        description = header[pipe_index + 1:].strip()
        return accession, delimiter, description

    parts = header.split(None, 1)

    if len(parts) < 2:
        return header, "", ""

    return parts[0], " ", parts[1].strip()


def detect_fasta_header_delimiter(headers):
    """
    Detect the most frequent accession/description delimiter in raw FASTA
    headers.

    Parameters
    ----------
    headers : iterable[str]
        Raw header lines, with or without a leading '>'.

    Returns
    -------
    str
        Most frequent delimiter, e.g. ' |' or '|'.
    """
    delimiters = Counter()

    for header in headers:
        _accession, delimiter, _description = _split_header(header)

        if delimiter:
            delimiters[delimiter] += 1

    if not delimiters:
        return DEFAULT_DELIMITER

    return delimiters.most_common(1)[0][0]


def build_conformed_header(accession, description, delimiter=DEFAULT_DELIMITER):
    """
    Build a FASTA header conforming to the dominant header format.

    Returns
    -------
    str
        Header starting with '>', e.g. '>PZ405244.1 |description'.
    """
    description = str(description or "").strip().replace(";", "_")
    return f">{accession}{delimiter}{description}"


def normalise_sequence(sequence):
    """
    Apply the same normalisation used by the ETL seq_get parser.

    Uppercases, removes gaps, converts X to N, converts U to T and maps any
    other unexpected character to N.
    """
    iupac_nucleotides = set("ACGTURYSWKMBDHVN")

    sequence = str(sequence).upper().replace("-", "").replace("X", "N")
    sequence = sequence.replace("U", "T")

    return "".join(
        base if base in iupac_nucleotides else "N"
        for base in sequence
    )


# ============================================================
# ACCESSION SOURCES
# ============================================================

def load_accessions_from_db(db_path):
    """
    Compute accessions present in metadata or cluster_composition but missing
    from the sequences table.

    Parameters
    ----------
    db_path : str
        Path to a SQLite sequence database.

    Returns
    -------
    list[str]
        Sorted missing accession list.
    """
    conn = sqlite3.connect(db_path)

    try:
        query = """
            SELECT DISTINCT acc.accession
            FROM (
                SELECT accession FROM metadata
                UNION
                SELECT accession FROM cluster_composition
            ) acc
            WHERE NOT EXISTS (
                SELECT 1 FROM sequences s
                WHERE s.accession = acc.accession
            )
            ORDER BY acc.accession;
        """
        cursor = conn.execute(query)
        return [row[0] for row in cursor.fetchall()]
    finally:
        conn.close()


def load_accessions_from_file(path):
    """
    Load one accession per line from a text file.

    Comma/semicolon separated entries on a single line are tolerated too.
    """
    accessions = []

    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()

        if not line or line.startswith("#"):
            continue

        for separator in [",", ";", "\t"]:
            line = line.replace(separator, "\n")

        for entry in line.split("\n"):
            accession = entry.strip().lstrip(">").strip()
            if accession:
                accessions.append(accession)

    return accessions


def load_accessions_from_text(text):
    """Parse a comma/semicolon/whitespace separated accession list."""
    accessions = []

    for separator in [",", ";", "\n", "\t", " "]:
        text = text.replace(separator, "\n")

    seen = set()

    for entry in text.split("\n"):
        accession = entry.strip().lstrip(">").strip()

        if not accession or accession in seen:
            continue

        seen.add(accession)
        accessions.append(accession)

    return accessions


# ============================================================
# PROGRESS TRACKING
# ============================================================

def load_progress(progress_file):
    """
    Load the set of already-fetched accessions from a progress JSON file.

    Returns
    -------
    set[str]
    """
    path = Path(progress_file)

    if not path.exists():
        return set()

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return set(data.get("fetched", []))
    except (json.JSONDecodeError, OSError):
        return set()


def save_progress(progress_file, fetched):
    """
    Persist the set of already-fetched accessions to a progress JSON file.
    """
    path = Path(progress_file)
    path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "fetched": sorted(fetched),
        "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    tmp_path = path.with_suffix(path.suffix + ".tmp")

    tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def seed_fetched_from_fasta(output_fasta):
    """
    If no progress file exists, seed the fetched set from accessions already
    present in the output FASTA so an interrupted run never writes duplicates.
    """
    path = Path(output_fasta)

    if not path.exists():
        return set()

    fetched = set()

    try:
        with path.open("r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                fetched.add(str(record.id))
    except Exception:
        pass

    return fetched


# ============================================================
# NCBI FETCHING
# ============================================================

def _entrez_handle(accessions, email, api_key):
    Entrez.email = email

    if api_key:
        Entrez.api_key = api_key

    return Entrez.efetch(
        db="nucleotide",
        id=",".join(accessions),
        rettype="fasta",
        retmode="text",
    )


def fetch_batch(accessions, email, api_key, sleep_seconds, max_retries=5):
    """
    Fetch one batch of accessions from NCBI, retrying on transient errors.

    Returns
    -------
    dict
        accession -> (description, raw_sequence)
    """
    results = {}

    for attempt in range(max_retries):
        try:
            handle = _entrez_handle(accessions, email, api_key)
            text = handle.read()
            handle.close()
            break
        except (urllib.error.HTTPError, urllib.error.URLError) as exc:
            wait = 2 ** attempt

            if attempt == max_retries - 1:
                print(
                    f"  Failed batch of {len(accessions)} accessions after "
                    f"{max_retries} retries: {exc}"
                )
                return results

            print(f"  Retrying batch in {wait}s ({exc})")
            time.sleep(wait)
    else:
        return results

    for record in SeqIO.parse(StringIO(text), "fasta"):
        accession = str(record.id)
        results[accession] = (str(record.description), str(record.seq))

    if sleep_seconds:
        time.sleep(sleep_seconds)

    return results


def fetch_missing(
    accessions,
    output_fasta,
    progress_file,
    email,
    api_key=None,
    batch_size=200,
    sleep_seconds=None,
    delimiter=DEFAULT_DELIMITER,
):
    """
    Fetch missing accessions from NCBI, resumable across runs.

    Accessions already recorded in the progress file (or already present in
    the output FASTA when no progress file exists) are skipped. Fetched
    sequences are appended to the output FASTA with headers conformed to the
    given delimiter.

    Parameters
    ----------
    accessions : list[str]
        Accessions to fetch.
    output_fasta : str
        FASTA file to append fetched sequences to.
    progress_file : str
        JSON progress file.
    email : str
        NCBI Entrez email (required).
    api_key : str or None
        NCBI Entrez API key.
    batch_size : int
        Accessions per efetch request.
    sleep_seconds : float or None
        Delay between requests. Defaults to 0.34 without an API key and
        0.1 with one.
    delimiter : str
        Header delimiter used to conform fetched headers.

    Returns
    -------
    tuple
        (fetched_dict, still_missing)
        fetched_dict maps accession -> (description, normalised_sequence).
        still_missing is the list of accessions that could not be retrieved.
    """
    if not accessions:
        return {}, []

    fetched = load_progress(progress_file)

    if not fetched:
        fetched = seed_fetched_from_fasta(output_fasta)

    pending = [
        accession
        for accession in accessions
        if accession not in fetched
    ]

    if not pending:
        print("Nothing to fetch: all requested accessions are already done.")
        return {}, []

    if sleep_seconds is None:
        sleep_seconds = 0.1 if api_key else 0.34

    output_path = Path(output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fetched_dict = {}
    failed = []

    out_handle = output_path.open("a", encoding="utf-8")

    try:
        for start in range(0, len(pending), batch_size):
            batch = pending[start:start + batch_size]

            print(
                f"  Fetching {start + 1}–{start + len(batch)} of {len(pending)} "
                f"({batch[0]} ... {batch[-1]})"
            )

            batch_results = fetch_batch(
                accessions=batch,
                email=email,
                api_key=api_key,
                sleep_seconds=sleep_seconds,
            )

            for accession in batch:
                record = batch_results.get(accession)

                if record is None:
                    failed.append(accession)
                    continue

                description, raw_sequence = record
                sequence = normalise_sequence(raw_sequence)

                if not sequence:
                    failed.append(accession)
                    continue

                header = build_conformed_header(
                    accession,
                    description,
                    delimiter=delimiter,
                )

                out_handle.write(header + "\n")

                for i in range(0, len(sequence), 80):
                    out_handle.write(sequence[i:i + 80] + "\n")

                fetched.add(accession)
                fetched_dict[accession] = (description, sequence)

            save_progress(progress_file, fetched)

    except KeyboardInterrupt:
        print("\nInterrupted. Progress saved; re-run to resume.")
        save_progress(progress_file, fetched)
        raise

    finally:
        out_handle.close()

    if failed:
        print(
            f"  Could not fetch {len(failed)} accessions "
            f"(not found or network error). First examples: {failed[:10]}"
        )

    return fetched_dict, failed


def fetch_and_build_sequence_rows(
    accessions,
    output_fasta,
    progress_file,
    email,
    api_key=None,
    batch_size=200,
    sleep_seconds=None,
    delimiter=DEFAULT_DELIMITER,
):
    """
    Fetch missing accessions and return ETL-ready sequence row dictionaries.

    Convenience wrapper for db/import_data.py. Returns rows compatible with
    the sequence_rows format produced by transform_seqdict.
    """
    fetched_dict, _failed = fetch_missing(
        accessions=accessions,
        output_fasta=output_fasta,
        progress_file=progress_file,
        email=email,
        api_key=api_key,
        batch_size=batch_size,
        sleep_seconds=sleep_seconds,
        delimiter=delimiter,
    )

    rows = []

    for accession, (description, sequence) in fetched_dict.items():
        counts = count_ambiguous_bases(sequence)
        rows.append({
            "accession": accession,
            "description": description,
            "sequence": sequence,
            "n_N": counts["n_N"],
            "n_degenerate": counts["n_degenerate"],
            "degenerate_breakdown": degenerate_breakdown_json(sequence),
        })

    return rows


# ============================================================
# STANDALONE CLI
# ============================================================

def parse_args():
    """Parse command-line arguments for the standalone fetch CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Resumable NCBI FASTA download for accessions missing from a "
            "sequence database. Run again after an interruption to resume."
        )
    )

    source_group = parser.add_mutually_exclusive_group(required=True)

    source_group.add_argument(
        "--db",
        help=(
            "SQLite database. Missing accessions are computed as "
            "(metadata UNION cluster_composition) MINUS sequences."
        )
    )

    source_group.add_argument(
        "--accessions-file",
        help="Text file with one accession per line."
    )

    source_group.add_argument(
        "--accessions",
        help="Inline accessions separated by commas, semicolons or whitespace."
    )

    parser.add_argument(
        "--fasta",
        help=(
            "Original FASTA file, used only to detect the dominant header "
            "format so fetched headers conform to it."
        )
    )

    parser.add_argument(
        "--output",
        default="data/fetched_sequences.fasta",
        help="FASTA file fetched sequences are appended to."
    )

    parser.add_argument(
        "--progress",
        default="data/fetch_progress.json",
        help="JSON progress file used to resume an interrupted download."
    )

    parser.add_argument(
        "--email",
        required=True,
        help="Email address required by NCBI Entrez."
    )

    parser.add_argument(
        "--api-key",
        default=None,
        help="NCBI Entrez API key (optional, raises rate limits)."
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="Accessions per efetch request. Default: 200."
    )

    parser.add_argument(
        "--sleep",
        type=float,
        default=None,
        help=(
            "Seconds to sleep between requests. Defaults to 0.34 without an "
            "API key, 0.1 with one."
        )
    )

    parser.add_argument(
        "--delimiter",
        default=None,
        help="Force a header delimiter instead of auto-detecting it."
    )

    return parser.parse_args()


def _load_accessions(args):
    if args.db:
        return load_accessions_from_db(args.db)

    if args.accessions_file:
        return load_accessions_from_file(args.accessions_file)

    return load_accessions_from_text(args.accessions)


def main():
    args = parse_args()

    accessions = _load_accessions(args)
    print(f"Accessions to fetch: {len(accessions)}")

    if not accessions:
        print("No missing accessions to fetch.")
        return

    if args.delimiter is not None:
        delimiter = args.delimiter
    elif args.fasta:
        headers = []
        with open(args.fasta, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    headers.append(line.strip())
        delimiter = detect_fasta_header_delimiter(headers)
        print(f"Detected header delimiter: {delimiter!r}")
    else:
        delimiter = DEFAULT_DELIMITER
        print(f"Using default header delimiter: {delimiter!r}")

    fetched_dict, still_missing = fetch_missing(
        accessions=accessions,
        output_fasta=args.output,
        progress_file=args.progress,
        email=args.email,
        api_key=args.api_key,
        batch_size=args.batch_size,
        sleep_seconds=args.sleep,
        delimiter=delimiter,
    )

    print("Fetch summary")
    print("-------------")
    print(f"Fetched: {len(fetched_dict)}")
    print(f"Still missing: {len(still_missing)}")
    print(f"Output FASTA: {args.output}")
    print(f"Progress file: {args.progress}")

    if still_missing:
        print("Still-missing accessions:")
        for accession in still_missing[:50]:
            print(f"  {accession}")


if __name__ == "__main__":
    main()
