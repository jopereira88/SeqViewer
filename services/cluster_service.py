"""
Cluster service layer.

Business logic for loading and processing cluster data.
"""

import db.queries as _q


def load_cluster_members(db_path, organism, cluster_number):
    """
    Load all members of a selected cluster.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    organism : str
        Cluster organism.

    cluster_number : str
        Cluster number.

    Returns
    -------
    pandas.DataFrame
        Cluster member dataframe.
    """
    return _q.load_cluster_members(
        db_path=db_path,
        organism=organism,
        cluster_number=cluster_number,
    )


def load_cluster_member_sequences(db_path, organism, cluster_number):
    """
    Load sequences for all members of a selected cluster.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    organism : str
        Cluster organism.

    cluster_number : str
        Cluster number.

    Returns
    -------
    pandas.DataFrame
        Cluster member sequences dataframe.
    """
    return _q.load_cluster_member_sequences(
        db_path=db_path,
        organism=organism,
        cluster_number=cluster_number,
    )


def search_clusters(
    db_path,
    organism_search="",
    cluster_search="",
    centroid_search="",
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All",
    min_n_sequences=None,
    max_n_sequences=None,
    limit=500,
    offset=0
):
    """
    Search clusters using precomputed derived tables.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    organism_search : str
        Search term for organism name.

    cluster_search : str
        Search term for cluster number.

    centroid_search : str
        Search term for centroid accession.

    segment_filter : str
        Segment filter or 'All'.

    genotype_filter : str
        Genotype filter or 'All'.

    ha_subtype_filter : str
        HA subtype filter or 'All'.

    na_subtype_filter : str
        NA subtype filter or 'All'.

    host_filter : str
        Host filter or 'All'.

    country_filter : str
        Country filter or 'All'.

    min_n_sequences : int or None
        Minimum number of sequences.

    max_n_sequences : int or None
        Maximum number of sequences.

    limit : int
        Maximum number of clusters to return.

    offset : int
        Number of rows to skip for pagination.

    Returns
    -------
    pandas.DataFrame
        Cluster search results.
    """
    return _q.search_clusters(
        db_path=db_path,
        organism_search=organism_search,
        cluster_search=cluster_search,
        centroid_search=centroid_search,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        min_n_sequences=min_n_sequences,
        max_n_sequences=max_n_sequences,
        limit=limit,
        offset=offset,
    )


def count_clusters(
    db_path,
    organism_search="",
    cluster_search="",
    centroid_search="",
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All",
    min_n_sequences=None,
    max_n_sequences=None,
):
    """
    Count clusters matching search filters without loading result rows.

    Returns
    -------
    int
        Total number of matching clusters.
    """
    return _q.count_clusters(
        db_path=db_path,
        organism_search=organism_search,
        cluster_search=cluster_search,
        centroid_search=centroid_search,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        min_n_sequences=min_n_sequences,
        max_n_sequences=max_n_sequences,
    )


def load_cluster_summary_for_cluster_subset(
    db_path,
    cluster_pairs,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Load cluster-level summary statistics for a subset of clusters.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    cluster_pairs : list[tuple[str, str]]
        List of (organism, cluster_number) pairs.

    segment_filter : str
        Segment filter or 'All'.

    genotype_filter : str
        Genotype filter or 'All'.

    ha_subtype_filter : str
        HA subtype filter or 'All'.

    na_subtype_filter : str
        NA subtype filter or 'All'.

    host_filter : str
        Host filter or 'All'.

    country_filter : str
        Country filter or 'All'.

    Returns
    -------
    pandas.DataFrame
        Cluster summary dataframe.
    """
    return _q.load_cluster_summary_for_cluster_subset(
        db_path=db_path,
        cluster_pairs=cluster_pairs,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
    )


def load_cluster_descriptor_table_for_cluster_subset(
    db_path,
    cluster_pairs,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Load a cluster-level descriptor table for a subset of clusters.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    cluster_pairs : list[tuple[str, str]]
        List of (organism, cluster_number) pairs.

    segment_filter : str
        Segment filter or 'All'.

    genotype_filter : str
        Genotype filter or 'All'.

    ha_subtype_filter : str
        HA subtype filter or 'All'.

    na_subtype_filter : str
        NA subtype filter or 'All'.

    host_filter : str
        Host filter or 'All'.

    country_filter : str
        Country filter or 'All'.

    Returns
    -------
    pandas.DataFrame
        Cluster descriptor dataframe.
    """
    return _q.load_cluster_descriptor_table_for_cluster_subset(
        db_path=db_path,
        cluster_pairs=cluster_pairs,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
    )


def load_cluster_background_counts_for_cluster_subset(
    db_path,
    cluster_pairs,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All",
    background_type="All"
):
    """
    Count background values for a subset of clusters.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    cluster_pairs : list[tuple[str, str]]
        List of (organism, cluster_number) pairs.

    segment_filter : str
        Segment filter or 'All'.

    genotype_filter : str
        Genotype filter or 'All'.

    ha_subtype_filter : str
        HA subtype filter or 'All'.

    na_subtype_filter : str
        NA subtype filter or 'All'.

    host_filter : str
        Host filter or 'All'.

    country_filter : str
        Country filter or 'All'.

    background_type : str
        Background type filter or 'All'.

    Returns
    -------
    pandas.DataFrame
        Background counts dataframe.
    """
    return _q.load_cluster_background_counts_for_cluster_subset(
        db_path=db_path,
        cluster_pairs=cluster_pairs,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        background_type=background_type,
    )


def load_cluster_length_stats_for_cluster_subset(
    db_path,
    cluster_pairs,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Load sequence length statistics for a subset of clusters.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    cluster_pairs : list[tuple[str, str]]
        List of (organism, cluster_number) pairs.

    segment_filter : str
        Segment filter or 'All'.

    genotype_filter : str
        Genotype filter or 'All'.

    ha_subtype_filter : str
        HA subtype filter or 'All'.

    na_subtype_filter : str
        NA subtype filter or 'All'.

    host_filter : str
        Host filter or 'All'.

    country_filter : str
        Country filter or 'All'.

    Returns
    -------
    pandas.DataFrame
        Length statistics dataframe.
    """
    return _q.load_cluster_length_stats_for_cluster_subset(
        db_path=db_path,
        cluster_pairs=cluster_pairs,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
    )


def load_cluster_member_sequences_for_cluster_subset(
    db_path,
    cluster_pairs,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Load member sequences for a subset of clusters.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    cluster_pairs : list[tuple[str, str]]
        List of (organism, cluster_number) pairs.

    segment_filter : str
        Segment filter or 'All'.

    genotype_filter : str
        Genotype filter or 'All'.

    ha_subtype_filter : str
        HA subtype filter or 'All'.

    na_subtype_filter : str
        NA subtype filter or 'All'.

    host_filter : str
        Host filter or 'All'.

    country_filter : str
        Country filter or 'All'.

    Returns
    -------
    pandas.DataFrame
        Member sequences dataframe.
    """
    return _q.load_cluster_member_sequences_for_cluster_subset(
        db_path=db_path,
        cluster_pairs=cluster_pairs,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
    )
