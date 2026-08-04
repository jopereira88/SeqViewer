"""
FASTA formatting utilities.
"""

import io
import re
import zipfile


def format_fasta(accession, description, sequence, line_width=80):
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

    Returns
    -------
    str
        Formatted FASTA string.
    """
    if sequence is None:
        return ""

    if description:
        header = f">{accession} | {description}"
    else:
        header = f">{accession}"

    sequence = str(sequence).replace("\n", "").replace(" ", "")

    lines = [header]

    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    return "\n".join(lines) + "\n"


def format_multifasta(df, line_width=80):
    """
    Format a DataFrame of sequences as a multifasta string.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain at least 'accession' and 'sequence' columns.
        May optionally contain 'description'.

    line_width : int
        Maximum line width for sequence wrapping.

    Returns
    -------
    str
        Concatenated FASTA entries.
    """
    parts = []

    for _, row in df.iterrows():
        entry = format_fasta(
            accession=row.get("accession", ""),
            description=row.get("description"),
            sequence=row.get("sequence"),
            line_width=line_width,
        )
        if entry:
            parts.append(entry)

    return "".join(parts)


def _safe_cluster_filename(cluster_number):
    """Convert a cluster number string to a safe filename."""
    name = str(cluster_number).replace(">", "").replace(" ", "_")
    name = re.sub(r"[^\w\-.]", "_", name)
    return name


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
    buf = io.BytesIO()

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for (organism, cluster_number), group in members_df.groupby(
            ["organism", "cluster_number"]
        ):
            fasta_text = format_multifasta(group, line_width=line_width)
            if not fasta_text:
                continue

            fname = f"{_safe_cluster_filename(cluster_number)}.fasta"
            zf.writestr(fname, fasta_text)

    return buf.getvalue()
