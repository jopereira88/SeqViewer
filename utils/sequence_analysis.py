"""
Sequence composition analysis utilities.

Pure functions (no Streamlit or database dependencies) shared by the ETL
pipeline and the viewer application.
"""

import json
from collections import Counter

# IUPAC ambiguity codes for nucleotides, excluding 'N' (unknown/any base).
DEGENERATE_BASES = "RYWSKMBDHV"

# Degenerate bases that are absent from a sequence still appear in the
# breakdown JSON so downstream code can rely on a stable set of keys.
_EMPTY_PER_BASE = {base: 0 for base in DEGENERATE_BASES}


def count_ambiguous_bases(sequence):
    """
    Count ambiguous / degenerate bases in one nucleotide sequence.

    Parameters
    ----------
    sequence : str or None
        Nucleotide sequence.

    Returns
    -------
    dict
        {
            "length": int,
            "n_N": int,
            "n_degenerate": int,
            "per_base": {base: count, ...}
        }

    'n_degenerate' counts IUPAC ambiguity codes R, Y, S, W, K, M, B, D,
    H, V. 'n_N' counts 'N' separately. 'per_base' always contains every
    degenerate base, so the JSON serialised form is stable.
    """
    if sequence is None:
        return {
            "length": 0,
            "n_N": 0,
            "n_degenerate": 0,
            "per_base": dict(_EMPTY_PER_BASE),
        }

    counts = Counter(str(sequence).upper())

    per_base = dict(_EMPTY_PER_BASE)
    for base in DEGENERATE_BASES:
        per_base[base] = int(counts.get(base, 0))

    return {
        "length": int(sum(counts.values())),
        "n_N": int(counts.get("N", 0)),
        "n_degenerate": int(sum(per_base.values())),
        "per_base": per_base,
    }


def degenerate_breakdown_json(sequence):
    """
    Serialise the per-base degenerate counts of a sequence to JSON text.

    Parameters
    ----------
    sequence : str or None
        Nucleotide sequence.

    Returns
    -------
    str
        JSON object string, e.g. '{"R": 0, "Y": 1, ...}'.
    """
    per_base = count_ambiguous_bases(sequence)["per_base"]
    return json.dumps(per_base, sort_keys=True)


def parse_accession_list(text):
    """
    Parse a user-supplied list of accessions.

    Accessions may be separated by semicolons (the recommended delimiter),
    commas, newlines or whitespace. Leading '>' characters (pasted FASTA
    headers) are stripped. Empty entries are dropped and order is preserved.

    Parameters
    ----------
    text : str or None
        Raw user input.

    Returns
    -------
    list[str]
        Cleaned, de-duplicated accession list.
    """
    if text is None:
        return []

    text = str(text)

    for separator in [";", ",", "\n", "\t", " "]:
        text = text.replace(separator, "\n")

    accessions = []
    seen = set()

    for raw in text.split("\n"):
        accession = raw.strip().lstrip(">").strip()

        if not accession:
            continue

        if accession in seen:
            continue

        seen.add(accession)
        accessions.append(accession)

    return accessions
