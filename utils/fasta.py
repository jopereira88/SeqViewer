"""
FASTA formatting utilities.
"""

import io
import re
import zipfile
from itertools import repeat


def format_fasta(accession, description, sequence, line_width=80, assembly=None):
    """
    Format one sequence as FASTA text.

    Parameters
    ----------
    accession : str
        Sequence accession identifier.

    description : str or None
        Sequence description.

    sequence : str or None
        Nucleotide sequence.

    line_width : int
        Maximum line width for sequence wrapping.

    assembly : str or None
        Assembly accession, appended to the header when present.

    Returns
    -------
    str
        Formatted FASTA string.
    """
    return _fasta_entry(
        accession=accession,
        description=description,
        sequence=sequence,
        line_width=line_width,
        assembly=assembly,
    )


def _fasta_entry(accession, description, sequence, line_width=80, assembly=None):
    """
    Build a single FASTA entry from already-extracted column values.

    Parameters
    ----------
    accession : str
        Sequence accession identifier.

    description : str or None
        Sequence description.

    sequence : str or None
        Nucleotide sequence.

    line_width : int
        Maximum line width for sequence wrapping.

    assembly : str or None
        Assembly accession, appended to the header when present.

    Returns
    -------
    str
        Formatted FASTA entry, empty when the sequence is missing.
    """
    if sequence is None:
        return ""

    header = f">{accession}"

    if description:
        header += f" | {description}"

    # NaN is truthy, so missing assemblies must be filtered explicitly.
    if assembly is not None:
        assembly = str(assembly).strip()

        if assembly and assembly.lower() not in ("nan", "none", "null"):
            header += f" | assembly={assembly}"

    sequence = str(sequence).replace("\n", "").replace(" ", "")

    return header + "\n" + "\n".join(
        sequence[i:i + line_width] for i in range(0, len(sequence), line_width)
    ) + "\n"


def format_multifasta(df, line_width=80):
    """
    Format a DataFrame of sequences as a multifasta string.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain at least 'accession' and 'sequence' columns.
        May optionally contain 'description' and 'assembly'.

    line_width : int
        Maximum line width for sequence wrapping.

    Returns
    -------
    str
        Concatenated FASTA entries.
    """
    if df is None or df.empty:
        return ""

    columns = df.columns

    accessions = df["accession"] if "accession" in columns else repeat("", len(df))
    descriptions = df["description"] if "description" in columns else repeat(None, len(df))
    sequences = df["sequence"] if "sequence" in columns else repeat(None, len(df))
    assemblies = df["assembly"] if "assembly" in columns else repeat(None, len(df))

    # zip over the columns directly: iterrows() builds a Series per row, which
    # costs about as much as the formatting itself on large exports.
    parts = []

    for accession, description, sequence, assembly in zip(
        accessions,
        descriptions,
        sequences,
        assemblies,
    ):
        entry = _fasta_entry(
            accession=accession,
            description=description,
            sequence=sequence,
            line_width=line_width,
            assembly=assembly,
        )
        if entry:
            parts.append(entry)

    return "".join(parts)


def _safe_cluster_filename(cluster_number):
    """Convert a cluster number string to a safe filename."""
    name = str(cluster_number).replace(">", "").replace(" ", "_")
    name = re.sub(r"[^\w\-.]", "_", name)
    return name


def _safe_organism_filename(organism):
    """Convert an organism name to a safe filename fragment."""
    name = str(organism).replace(">", "").replace(" ", "_")
    name = re.sub(r"[^\w\-.]", "_", name)
    return name


def _cluster_zip_filenames(keys):
    """
    Build a unique archive filename for every cluster in the export.

    Cluster numbers are only unique within an organism, so a multi-organism
    export needs the organism in the name or entries would overwrite each
    other in the archive. Single-organism exports keep the plain cluster
    name.

    Parameters
    ----------
    keys : iterable[tuple[str, str]]
        (organism, cluster_number) pairs in the order entries should be
        written. Every cluster in the export needs a key, whether or not it
        ends up with sequences, so a streaming caller can pass its cluster
        list directly.

    Returns
    -------
    tuple[list, list]
        Keys and the filename stem for each.
    """
    keys = list(keys)

    names = [_safe_cluster_filename(cluster_number) for _, cluster_number in keys]

    if len(set(names)) == len(names):
        return keys, names

    disambiguated = [
        f"{_safe_organism_filename(organism)}_{name}"
        for (organism, _), name in zip(keys, names)
    ]

    return keys, disambiguated


def build_cluster_zip(members_df, line_width=80):
    """
    Build an in-memory zip with one .fasta file per cluster.

    Parameters
    ----------
    members_df : pandas.DataFrame
        Must contain 'organism', 'cluster_number', 'accession',
        'description', 'sequence'.

    line_width : int
        Maximum line width for sequence wrapping.

    Returns
    -------
    bytes
        Zip file contents.
    """
    if members_df is None or members_df.empty:
        return b""

    keys, names = _cluster_zip_filenames(
        members_df.groupby(["organism", "cluster_number"]).groups
    )

    groups = {
        key: group
        for key, group in members_df.groupby(["organism", "cluster_number"])
    }

    buf = io.BytesIO()

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for key, fname in zip(keys, names):
            fasta_text = format_multifasta(groups[key], line_width=line_width)
            if not fasta_text:
                continue

            zf.writestr(f"{fname}.fasta", fasta_text)

    return buf.getvalue()


def build_cluster_zip_from_rows(cluster_keys, rows, line_width=80):
    """
    Build an in-memory zip with one .fasta file per cluster, streaming rows.

    Equivalent output to build_cluster_zip, but the sequences are consumed
    one at a time instead of through a dataframe. A whole-dataset export is
    over a million sequences, so the dataframe version needs several
    gigabytes while this one stays in the low hundreds of megabytes.

    Rows must be grouped by cluster and in cluster_keys order. Only the
    cluster currently being read is held in memory.

    Parameters
    ----------
    cluster_keys : iterable[tuple[str, str]]
        Every (organism, cluster_number) in the export, in the order the
        rows will arrive. Needed up front so filenames can be assigned
        before any sequence has been read.

    rows : iterable[tuple]
        (organism, cluster_number, accession, description, assembly,
        sequence) rows, ordered by cluster.

    line_width : int
        Maximum line width for sequence wrapping.

    Returns
    -------
    tuple[bytes, int]
        Zip file contents and the number of sequence rows written.
    """
    keys, names = _cluster_zip_filenames(cluster_keys)
    filenames = dict(zip(keys, names))

    def _filename_for(key):
        # A cluster missing from cluster_keys should still be exported under
        # its plain name rather than dropped.
        return filenames.get(key) or _safe_cluster_filename(key[1])

    buf = io.BytesIO()

    n_rows = 0

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        current_key = None
        parts = []

        for organism, cluster_number, accession, description, assembly, sequence in rows:
            key = (organism, cluster_number)

            if key != current_key:
                if parts:
                    zf.writestr(
                        f"{_filename_for(current_key)}.fasta",
                        "".join(parts),
                    )

                current_key = key
                parts = []

            n_rows += 1

            entry = _fasta_entry(
                accession=accession,
                description=description,
                sequence=sequence,
                line_width=line_width,
                assembly=assembly,
            )

            if entry:
                parts.append(entry)

        if parts:
            zf.writestr(
                f"{_filename_for(current_key)}.fasta",
                "".join(parts),
            )

    return buf.getvalue(), n_rows
