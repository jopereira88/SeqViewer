"""
Utility functions for the sequence viewer.
"""

from utils.fasta import (
    format_fasta,
    format_multifasta,
    build_cluster_zip,
    build_cluster_zip_from_rows,
)
from utils.tables import (
    df_to_delimited_bytes,
    df_to_tsv_bytes,
    df_to_csv_bytes,
    build_delimited_zip,
    build_delimited_zip_from_rows,
    cluster_metadata_zip_bytes,
    cluster_metadata_zip_bytes_from_rows,
    rows_to_tsv_stream,
)
from utils.sequence_analysis import (
    DEGENERATE_BASES,
    ambiguous_breakdown_json,
    count_ambiguous_bases,
    degenerate_breakdown_json,
    parse_accession_list,
)

__all__ = [
    "format_fasta",
    "format_multifasta",
    "build_cluster_zip",
    "build_cluster_zip_from_rows",
    "df_to_delimited_bytes",
    "df_to_tsv_bytes",
    "df_to_csv_bytes",
    "build_delimited_zip",
    "build_delimited_zip_from_rows",
    "cluster_metadata_zip_bytes",
    "cluster_metadata_zip_bytes_from_rows",
    "rows_to_tsv_stream",
    "DEGENERATE_BASES",
    "ambiguous_breakdown_json",
    "count_ambiguous_bases",
    "degenerate_breakdown_json",
    "parse_accession_list",
]
