"""
Services package.

Provides the business logic layer for the sequence viewer application.
"""

from services.sequence_service import (
    get_table_counts,
    count_sequences,
    load_sequence_index,
    load_sequence_detail,
    load_sequences_by_accessions,
    load_sequence_records_by_accessions,
)

from services.metadata_service import (
    get_distinct_values,
    get_filter_options,
)

from services.cluster_service import (
    load_cluster_members,
    load_cluster_member_sequences,
    search_clusters,
    count_clusters,
    load_cluster_summary_for_cluster_subset,
    load_cluster_descriptor_table_for_cluster_subset,
    load_cluster_background_counts_for_cluster_subset,
    load_cluster_length_stats_for_cluster_subset,
    load_cluster_member_sequences_for_cluster_subset,
    load_cluster_metadata_members_for_cluster_subset,
    load_cluster_metadata_clusters_for_cluster_subset,
    iter_cluster_member_sequences_for_cluster_subset,
    iter_cluster_metadata_members_for_cluster_subset,
)

from services.assembly_service import (
    get_assembly_links,
    load_assembly_link_index,
    load_cluster_assembly_values,
    AssemblyLinkIndex,
)

__all__ = [
    # Sequence
    "get_table_counts",
    "count_sequences",
    "load_sequence_index",
    "load_sequence_detail",
    "load_sequences_by_accessions",
    "load_sequence_records_by_accessions",
    # Metadata
    "get_distinct_values",
    "get_filter_options",
    # Cluster
    "load_cluster_members",
    "load_cluster_member_sequences",
    "search_clusters",
    "count_clusters",
    "load_cluster_summary_for_cluster_subset",
    "load_cluster_descriptor_table_for_cluster_subset",
    "load_cluster_background_counts_for_cluster_subset",
    "load_cluster_length_stats_for_cluster_subset",
    "load_cluster_member_sequences_for_cluster_subset",
    "load_cluster_metadata_members_for_cluster_subset",
    "load_cluster_metadata_clusters_for_cluster_subset",
    "iter_cluster_member_sequences_for_cluster_subset",
    "iter_cluster_metadata_members_for_cluster_subset",
    # Assembly
    "get_assembly_links",
    "load_assembly_link_index",
    "load_cluster_assembly_values",
    "AssemblyLinkIndex",
]
