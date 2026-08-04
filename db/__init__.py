"""
Database package.

Provides SQLite connection management and query functions for the
sequence viewer application.
"""

from db.connection import connect_sqlite
from db.queries import (
    get_table_counts,
    get_distinct_values,
    count_sequences,
    load_sequence_index,
    load_sequence_detail,
    load_sequences_by_accessions,
    load_sequence_records_by_accessions,
    load_cluster_members,
    load_cluster_member_sequences,
    search_clusters,
    count_clusters,
    load_cluster_member_sequences_for_cluster_subset,
    load_cluster_summary_for_cluster_subset,
    load_cluster_background_counts_for_cluster_subset,
    load_cluster_length_stats_for_cluster_subset,
    load_cluster_descriptor_table_for_cluster_subset,
)

__all__ = [
    "connect_sqlite",
    "get_table_counts",
    "get_distinct_values",
    "count_sequences",
    "load_sequence_index",
    "load_sequence_detail",
    "load_sequences_by_accessions",
    "load_sequence_records_by_accessions",
    "load_cluster_members",
    "load_cluster_member_sequences",
    "search_clusters",
    "count_clusters",
    "load_cluster_member_sequences_for_cluster_subset",
    "load_cluster_summary_for_cluster_subset",
    "load_cluster_background_counts_for_cluster_subset",
    "load_cluster_length_stats_for_cluster_subset",
    "load_cluster_descriptor_table_for_cluster_subset",
]
