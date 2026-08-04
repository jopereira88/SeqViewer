"""
Sequence service layer.

Business logic for loading and processing sequence data.
"""

import db.queries as _q


def get_table_counts(db_path):
    """
    Count rows in the main and derived database tables.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    Returns
    -------
    dict
        Table name to row count mapping.
    """
    return _q.get_table_counts(db_path)


def count_sequences(
    db_path,
    accession_search="",
    organism_search="",
    species_search="",
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All",
    min_degenerate=0,
    max_degenerate=None,
    min_n=0,
    max_n=None,
    centroids_only=False,
):
    """
    Count sequences matching search filters.

    Returns an integer count.
    """
    return _q.count_sequences(
        db_path=db_path,
        accession_search=accession_search,
        organism_search=organism_search,
        species_search=species_search,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        min_degenerate=min_degenerate,
        max_degenerate=max_degenerate,
        min_n=min_n,
        max_n=max_n,
        centroids_only=centroids_only,
    )


def load_sequence_index(
    db_path,
    accession_search="",
    organism_search="",
    species_search="",
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All",
    min_degenerate=0,
    max_degenerate=None,
    min_n=0,
    max_n=None,
    limit=1000,
    offset=0,
    centroids_only=False,
):
    """
    Load the main sequence index table with filters applied.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    accession_search : str
        Search term for accession, description or organism name.

    organism_search : str
        Search term for organism name.

    species_search : str
        Search term for species.

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

    min_degenerate : int
        Minimum number of degenerate bases (RYSWKMBDHV).

    max_degenerate : int or None
        Maximum number of degenerate bases, or None for unbounded.

    min_n : int
        Minimum number of N bases.

    max_n : int or None
        Maximum number of N bases, or None for unbounded.

    limit : int
        Maximum number of rows.

    offset : int
        Number of rows to skip for pagination.

    centroids_only : bool
        If True, only return sequences that are cluster centroids.

    Returns
    -------
    pandas.DataFrame
        Sequence index dataframe.
    """
    return _q.load_sequence_index(
        db_path=db_path,
        accession_search=accession_search,
        organism_search=organism_search,
        species_search=species_search,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        min_degenerate=min_degenerate,
        max_degenerate=max_degenerate,
        min_n=min_n,
        max_n=max_n,
        limit=limit,
        offset=offset,
        centroids_only=centroids_only,
    )


def load_sequence_detail(db_path, accession):
    """
    Load one selected sequence and its metadata.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    accession : str
        Selected accession.

    Returns
    -------
    dict or None
        Selected sequence record, or None if not found.
    """
    return _q.load_sequence_detail(db_path=db_path, accession=accession)


def load_sequences_by_accessions(db_path, accessions):
    """
    Load sequences for a list of accessions.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    accessions : list[str]
        List of accession identifiers.

    Returns
    -------
    pandas.DataFrame
        Dataframe with accession, description and sequence.
    """
    return _q.load_sequences_by_accessions(db_path=db_path, accessions=accessions)


def load_sequence_records_by_accessions(
    db_path,
    accessions,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All",
    min_degenerate=0,
    max_degenerate=None,
    min_n=0,
    max_n=None,
):
    """
    Load full metadata and sequence records for an explicit accession list.

    Applies the sidebar metadata and degenerate/N filters.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    accessions : list[str]
        List of accession identifiers.

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

    min_degenerate : int
        Minimum number of degenerate bases.

    max_degenerate : int or None
        Maximum number of degenerate bases, or None for unbounded.

    min_n : int
        Minimum number of N bases.

    max_n : int or None
        Maximum number of N bases, or None for unbounded.

    Returns
    -------
    pandas.DataFrame
        Accession records with metadata, sequence and base counts.
    """
    return _q.load_sequence_records_by_accessions(
        db_path=db_path,
        accessions=accessions,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        min_degenerate=min_degenerate,
        max_degenerate=max_degenerate,
        min_n=min_n,
        max_n=max_n,
    )
