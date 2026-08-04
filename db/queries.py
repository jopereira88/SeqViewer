"""
Database query functions for the Streamlit sequence viewer.

This module contains only read/query functions.
No Streamlit UI code should be placed here.
"""

import sqlite3

import pandas as pd
import streamlit as st

from db.connection import connect_sqlite  # noqa: F401 — re-exported


@st.cache_data(show_spinner=False)
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
        Dictionary with table names as keys and row counts as values.
    """
    conn = sqlite3.connect(db_path)

    tables = [
        "metadata",
        "sequences",
        "clusters",
        "cluster_composition",
        "cluster_summary",
        "cluster_background_counts",
        "cluster_descriptors",
        "cluster_filter_counts",
    ]

    counts = {}

    for table in tables:
        query = """
            SELECT name
            FROM sqlite_master
            WHERE type = 'table'
              AND name = ?;
        """

        exists_df = pd.read_sql_query(query, conn, params=[table])

        if exists_df.empty:
            counts[table] = 0
            continue

        count_query = f"SELECT COUNT(*) AS n FROM {table};"
        counts[table] = int(pd.read_sql_query(count_query, conn)["n"].iloc[0])

    conn.close()

    return counts


@st.cache_data(show_spinner=False)
def get_distinct_values(db_path, table, column):
    """
    Get distinct non-null values from a table column.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    table : str
        Table name.

    column : str
        Column name.

    Returns
    -------
    list[str]
        List beginning with 'All', followed by distinct values.
    """
    allowed_tables = {
        "metadata",
        "sequences",
        "clusters",
        "cluster_composition",
        "cluster_summary",
        "cluster_background_counts",
        "cluster_descriptors",
        "cluster_filter_counts",
    }

    allowed_columns = {
        "metadata": {
            "accession",
            "organism_name",
            "genbank_refseq",
            "assembly",
            "release_date",
            "isolate",
            "species",
            "length",
            "nuc_completeness",
            "genotype",
            "ha_subtype",
            "na_subtype",
            "segment",
            "country",
            "host",
            "collection_date",
            "comment",
        },
        "sequences": {
            "sequence_id",
            "accession",
            "description",
            "sequence",
        },
        "clusters": {
            "organism",
            "cluster_number",
            "centroid",
        },
        "cluster_composition": {
            "organism",
            "cluster_number",
            "accession",
            "identity_to_centroid",
        },
        "cluster_summary": {
            "organism",
            "cluster_number",
            "centroid",
            "segment",
            "n_segments",
            "n_sequences",
            "n_hosts",
            "n_countries",
            "n_genotypes",
        },
        "cluster_background_counts": {
            "organism",
            "cluster_number",
            "centroid",
            "segment",
            "background_type",
            "background_value",
            "n_sequences",
        },
        "cluster_descriptors": {
            "organism",
            "cluster_number",
            "representative",
            "n_sequences",
            "dominant_genotype",
            "dominant_segment",
            "dominant_host",
            "dominant_country",
            "n_genotypes",
            "n_segments",
            "n_hosts",
            "n_countries",
            "is_singleton",
            "is_mixed_genotype",
            "is_mixed_segment",
            "is_multi_host",
            "is_multi_country",
        },
        "cluster_filter_counts": {
            "organism",
            "cluster_number",
            "centroid",
            "segment",
            "genotype",
            "ha_subtype",
            "na_subtype",
            "host",
            "country",
            "n_sequences",
        },
    }

    if table not in allowed_tables:
        raise ValueError(f"Invalid table name: {table}")

    if column not in allowed_columns[table]:
        raise ValueError(f"Invalid column name for {table}: {column}")

    conn = sqlite3.connect(db_path)

    query = f"""
        SELECT DISTINCT {column}
        FROM {table}
        WHERE {column} IS NOT NULL
        ORDER BY {column};
    """

    df = pd.read_sql_query(query, conn)
    conn.close()

    values = df[column].dropna().astype(str).tolist()

    return ["All"] + values


# ============================================================
# INTERNAL HELPERS
# ============================================================

def _normalise_filter_value(value):
    """
    Normalise sidebar filter values.

    The app uses 'All' for no filter.
    Derived filter-count tables use 'Unknown' for missing values.
    """
    if value is None:
        return "All"

    value = str(value).strip()

    if value == "":
        return "All"

    return value


def _metadata_filters_active(
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Return True if any sequence-level metadata filter is active.
    """
    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    ha_subtype_filter = _normalise_filter_value(ha_subtype_filter)
    na_subtype_filter = _normalise_filter_value(na_subtype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    return any([
        segment_filter != "All",
        genotype_filter != "All",
        ha_subtype_filter != "All",
        na_subtype_filter != "All",
        host_filter != "All",
        country_filter != "All",
    ])


def _clean_cluster_pairs(cluster_pairs):
    """
    Clean a list of (organism, cluster_number) pairs.

    Parameters
    ----------
    cluster_pairs : list[tuple[str, str]]
        Raw cluster pairs.

    Returns
    -------
    list[tuple[str, str]]
        Cleaned cluster pairs.
    """
    clean_pairs = []

    if not cluster_pairs:
        return clean_pairs

    for organism, cluster_number in cluster_pairs:
        if organism is None or pd.isna(organism):
            continue

        if cluster_number is None or pd.isna(cluster_number):
            continue

        clean_pairs.append((str(organism), str(cluster_number)))

    return clean_pairs


def _build_cluster_pair_where_clause(table_alias, clean_pairs):
    """
    Build a SQL WHERE fragment for cluster pair matching.

    Parameters
    ----------
    table_alias : str
        SQL table alias containing organism and cluster_number.

    clean_pairs : list[tuple[str, str]]
        Cleaned cluster pairs.

    Returns
    -------
    tuple[str, list]
        WHERE fragment and SQL parameters.
    """
    pair_conditions = []
    params = []

    for organism, cluster_number in clean_pairs:
        pair_conditions.append(
            f"({table_alias}.organism = ? AND CAST({table_alias}.cluster_number AS TEXT) = ?)"
        )
        params.extend([organism, cluster_number])

    where_clause = " OR ".join(pair_conditions)

    return where_clause, params


def _append_metadata_filters(
    query,
    params,
    table_alias,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Append metadata filters to a SQL query.

    This helper is used with base metadata or cluster_member_metadata.
    """
    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    ha_subtype_filter = _normalise_filter_value(ha_subtype_filter)
    na_subtype_filter = _normalise_filter_value(na_subtype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    if segment_filter != "All":
        query += f" AND {table_alias}.segment = ?"
        params.append(segment_filter)

    if genotype_filter != "All":
        query += f" AND {table_alias}.genotype = ?"
        params.append(genotype_filter)

    if ha_subtype_filter != "All":
        query += f" AND {table_alias}.ha_subtype = ?"
        params.append(ha_subtype_filter)

    if na_subtype_filter != "All":
        query += f" AND {table_alias}.na_subtype = ?"
        params.append(na_subtype_filter)

    if host_filter != "All":
        query += f" AND {table_alias}.host = ?"
        params.append(host_filter)

    if country_filter != "All":
        query += f" AND {table_alias}.country = ?"
        params.append(country_filter)

    return query, params


def _append_ambiguous_base_filters(
    query,
    params,
    min_degenerate=0,
    max_degenerate=None,
    min_n=0,
    max_n=None,
    table_alias="s",
):
    """
    Append degenerate-base and N-content filters to a query over the
    sequences table.

    The counts are precomputed by the ETL and stored on the sequences table.
    """
    if min_degenerate is None:
        min_degenerate = 0

    if min_n is None:
        min_n = 0

    if max_degenerate is not None:
        query += f" AND {table_alias}.n_degenerate BETWEEN ? AND ?"
        params.extend([int(min_degenerate), int(max_degenerate)])
    elif int(min_degenerate) > 0:
        query += f" AND {table_alias}.n_degenerate >= ?"
        params.append(int(min_degenerate))

    if max_n is not None:
        query += f" AND {table_alias}.n_N BETWEEN ? AND ?"
        params.extend([int(min_n), int(max_n)])
    elif int(min_n) > 0:
        query += f" AND {table_alias}.n_N >= ?"
        params.append(int(min_n))

    return query, params


def _append_filter_count_filters(
    query,
    params,
    table_alias,
    segment_filter="All",
    genotype_filter="All",
    ha_subtype_filter="All",
    na_subtype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Append filters to a query over cluster_filter_counts.

    cluster_filter_counts stores missing metadata values as 'Unknown'.
    """
    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    ha_subtype_filter = _normalise_filter_value(ha_subtype_filter)
    na_subtype_filter = _normalise_filter_value(na_subtype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    if segment_filter != "All":
        query += f" AND {table_alias}.segment = ?"
        params.append(segment_filter)

    if genotype_filter != "All":
        query += f" AND {table_alias}.genotype = ?"
        params.append(genotype_filter)

    if ha_subtype_filter != "All":
        query += f" AND {table_alias}.ha_subtype = ?"
        params.append(ha_subtype_filter)

    if na_subtype_filter != "All":
        query += f" AND {table_alias}.na_subtype = ?"
        params.append(na_subtype_filter)

    if host_filter != "All":
        query += f" AND {table_alias}.host = ?"
        params.append(host_filter)

    if country_filter != "All":
        query += f" AND {table_alias}.country = ?"
        params.append(country_filter)

    return query, params


def _empty_cluster_subset_summary_df():
    """
    Return an empty cluster summary dataframe with stable columns.
    """
    return pd.DataFrame(
        columns=[
            "organism",
            "cluster_number",
            "centroid",
            "segment",
            "n_segments",
            "n_sequences",
            "full_cluster_size",
            "n_hosts",
            "n_countries",
            "n_genotypes",
            "length_min",
            "length_q1",
            "length_mean",
            "length_median",
            "length_q3",
            "length_max",
            "length_std",
        ]
    )


def _empty_cluster_descriptor_df():
    """
    Return an empty cluster descriptor dataframe with stable columns.

    For app compatibility, JSON columns are exposed as:
    - genotypes
    - segments
    - hosts
    - countries
    """
    return pd.DataFrame(
        columns=[
            "organism",
            "cluster_number",
            "representative",
            "n_sequences",
            "genotypes",
            "segments",
            "hosts",
            "countries",
            "dominant_genotype",
            "dominant_segment",
            "dominant_host",
            "dominant_country",
            "n_genotypes",
            "n_segments",
            "n_hosts",
            "n_countries",
            "is_singleton",
            "is_mixed_genotype",
            "is_mixed_segment",
            "is_multi_host",
            "is_multi_country",
            "length_min",
            "length_mean",
            "length_median",
            "length_max",
            "length_std",
        ]
    )


# ============================================================
# SEQUENCE-LEVEL QUERIES
# ============================================================

def _build_sequence_index_where(
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
    Build the WHERE clause and parameters for sequence index queries.

    Returns the fragment appended after ``WHERE 1 = 1`` and the
    corresponding parameter list.
    """
    clauses = []
    params = []

    if accession_search:
        clauses.append(
            "AND (m.accession LIKE ? OR s.description LIKE ? OR m.organism_name LIKE ?)"
        )
        pattern = f"%{accession_search}%"
        params.extend([pattern, pattern, pattern])

    if organism_search:
        clauses.append("AND m.organism_name LIKE ?")
        params.append(f"%{organism_search}%")

    if species_search:
        clauses.append("AND m.species LIKE ?")
        params.append(f"%{species_search}%")

    meta_fragment, meta_params = _append_metadata_filters(
        query="",
        params=[],
        table_alias="m",
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
    )
    if meta_fragment:
        clauses.append(meta_fragment)
    params.extend(meta_params)

    ambig_fragment, ambig_params = _append_ambiguous_base_filters(
        query="",
        params=[],
        min_degenerate=min_degenerate,
        max_degenerate=max_degenerate,
        min_n=min_n,
        max_n=max_n,
    )
    if ambig_fragment:
        clauses.append(ambig_fragment)
    params.extend(ambig_params)

    if centroids_only:
        clauses.append("AND m.accession = c.centroid")

    return " ".join(clauses), params


_SEQUENCE_INDEX_BASE_FROM = """
    FROM metadata m
    LEFT JOIN sequences s
        ON m.accession = s.accession
    LEFT JOIN cluster_composition cc
        ON m.accession = cc.accession
    LEFT JOIN clusters c
        ON cc.organism = c.organism
        AND cc.cluster_number = c.cluster_number
"""


@st.cache_data(show_spinner=False)
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
    Count sequences matching the search filters.

    Returns an integer count.
    """
    where_fragment, params = _build_sequence_index_where(
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

    query = f"""
        SELECT COUNT(*) AS n
        {_SEQUENCE_INDEX_BASE_FROM}
        WHERE 1 = 1 {where_fragment};
    """

    conn = sqlite3.connect(db_path)
    result = pd.read_sql_query(query, conn, params=params)
    conn.close()

    return int(result.iloc[0, 0])


@st.cache_data(show_spinner=False)
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
    Load the main sequence index table.

    This joins metadata, sequences and cluster composition.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    accession_search : str
        Search term applied to accession, description and organism name.

    organism_search : str
        Search term applied to organism_name.

    species_search : str
        Search term applied to species.

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
        Minimum number of degenerate bases (RYSWKMBDHV) per sequence.

    max_degenerate : int or None
        Maximum number of degenerate bases per sequence, or None for unbounded.

    min_n : int
        Minimum number of N bases per sequence.

    max_n : int or None
        Maximum number of N bases per sequence, or None for unbounded.

    limit : int
        Maximum number of rows to return.

    offset : int
        Number of rows to skip for pagination.

    centroids_only : bool
        If True, only return sequences that are cluster centroids.

    Returns
    -------
    pandas.DataFrame
        Sequence index dataframe.
    """
    limit = int(limit)
    offset = int(offset)

    where_fragment, params = _build_sequence_index_where(
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

    query = f"""
        SELECT
            m.accession,
            m.organism_name,
            m.genbank_refseq,
            m.assembly,
            m.release_date,
            m.species,
            m.length,
            m.nuc_completeness,
            m.genotype,
            m.ha_subtype,
            m.na_subtype,
            m.segment,
            m.country,
            m.host,
            m.collection_date,
            s.description,
            s.n_N,
            s.n_degenerate,
            cc.organism AS cluster_organism,
            cc.cluster_number,
            cc.identity_to_centroid,
            c.centroid
        {_SEQUENCE_INDEX_BASE_FROM}
        WHERE 1 = 1 {where_fragment}
        ORDER BY
            m.organism_name,
            m.segment,
            m.accession
        LIMIT ? OFFSET ?;
    """

    params.extend([limit, offset])

    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(query, conn, params=params)
    conn.close()

    return df


@st.cache_data(show_spinner=False)
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
    conn = sqlite3.connect(db_path)

    query = """
        SELECT
            m.*,
            s.description,
            s.sequence
        FROM metadata m
        LEFT JOIN sequences s
            ON m.accession = s.accession
        WHERE m.accession = ?;
    """

    df = pd.read_sql_query(query, conn, params=[accession])
    conn.close()

    if df.empty:
        return None

    return df.iloc[0].to_dict()


@st.cache_data(show_spinner=False)
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
    if not accessions:
        return pd.DataFrame(
            columns=[
                "accession",
                "description",
                "sequence"
            ]
        )

    accessions = [
        str(accession)
        for accession in accessions
        if accession is not None and not pd.isna(accession)
    ]

    if not accessions:
        return pd.DataFrame(
            columns=[
                "accession",
                "description",
                "sequence"
            ]
        )

    chunk_size = 900

    dfs = []
    conn = sqlite3.connect(db_path)

    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i:i + chunk_size]
        placeholders = ",".join(["?"] * len(chunk))

        query = f"""
            SELECT
                accession,
                description,
                sequence
            FROM sequences
            WHERE accession IN ({placeholders})
            ORDER BY accession;
        """

        chunk_df = pd.read_sql_query(query, conn, params=chunk)
        dfs.append(chunk_df)

    conn.close()

    if not dfs:
        return pd.DataFrame(
            columns=[
                "accession",
                "description",
                "sequence"
            ]
        )

    return pd.concat(dfs, ignore_index=True)


@st.cache_data(show_spinner=False)
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
    Load metadata and sequences for an explicit list of accessions.

    Unlike load_sequences_by_accessions, this returns full metadata records
    joined with the sequence text and the precomputed ambiguous base counts,
    and applies the sidebar metadata and degenerate/N filters.

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
        Rows with accession, organism_name, segment, genotype, host, country,
        length, description, sequence, n_N, n_degenerate and
        degenerate_breakdown.
    """
    output_columns = [
        "accession",
        "organism_name",
        "segment",
        "genotype",
        "ha_subtype",
        "na_subtype",
        "host",
        "country",
        "length",
        "description",
        "sequence",
        "n_N",
        "n_degenerate",
        "degenerate_breakdown",
    ]

    accessions = [
        str(accession)
        for accession in accessions
        if accession is not None and not pd.isna(accession)
    ]

    if not accessions:
        return pd.DataFrame(columns=output_columns)

    chunk_size = 900

    dfs = []
    conn = sqlite3.connect(db_path)

    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i:i + chunk_size]
        placeholders = ",".join(["?"] * len(chunk))

        query = f"""
            SELECT
                m.accession,
                m.organism_name,
                m.segment,
                m.genotype,
                m.ha_subtype,
                m.na_subtype,
                m.host,
                m.country,
                m.length,
                s.description,
                s.sequence,
                s.n_N,
                s.n_degenerate,
                s.degenerate_breakdown
            FROM metadata m
            JOIN sequences s
                ON m.accession = s.accession
            WHERE m.accession IN ({placeholders})
        """

        params = list(chunk)

        query, params = _append_metadata_filters(
            query=query,
            params=params,
            table_alias="m",
            segment_filter=segment_filter,
            genotype_filter=genotype_filter,
            ha_subtype_filter=ha_subtype_filter,
            na_subtype_filter=na_subtype_filter,
            host_filter=host_filter,
            country_filter=country_filter,
        )

        query, params = _append_ambiguous_base_filters(
            query=query,
            params=params,
            min_degenerate=min_degenerate,
            max_degenerate=max_degenerate,
            min_n=min_n,
            max_n=max_n,
            table_alias="s",
        )

        query += """
            ORDER BY m.accession;
        """

        chunk_df = pd.read_sql_query(query, conn, params=params)
        dfs.append(chunk_df)

    conn.close()

    if not dfs:
        return pd.DataFrame(columns=output_columns)

    return pd.concat(dfs, ignore_index=True)


# ============================================================
# SELECTED CLUSTER QUERIES
# ============================================================

@st.cache_data(show_spinner=False)
def load_cluster_members(db_path, organism, cluster_number):
    """
    Load all members of a selected cluster.

    Uses the cluster_member_metadata view.
    """
    conn = sqlite3.connect(db_path)

    query = """
        SELECT
            organism,
            cluster_number,
            accession,
            identity_to_centroid,
            centroid,
            organism_name,
            segment,
            genotype,
            host,
            country,
            collection_date,
            length
        FROM cluster_member_metadata
        WHERE organism = ?
          AND CAST(cluster_number AS TEXT) = ?
        ORDER BY
            identity_to_centroid DESC,
            accession;
    """

    df = pd.read_sql_query(
        query,
        conn,
        params=[organism, str(cluster_number)]
    )

    conn.close()

    return df


@st.cache_data(show_spinner=False)
def load_cluster_member_sequences(db_path, organism, cluster_number):
    """
    Load sequences for all members of a selected cluster.

    Uses the cluster_member_metadata view.
    """
    conn = sqlite3.connect(db_path)

    query = """
        SELECT
            accession,
            description,
            sequence,
            identity_to_centroid
        FROM cluster_member_metadata
        WHERE organism = ?
          AND CAST(cluster_number AS TEXT) = ?
        ORDER BY
            identity_to_centroid DESC,
            accession;
    """

    df = pd.read_sql_query(
        query,
        conn,
        params=[organism, str(cluster_number)]
    )

    conn.close()

    return df


# ============================================================
# CLUSTER SEARCH
# ============================================================

@st.cache_data(show_spinner=False)
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

    If no metadata filter is active, this uses cluster_summary.

    If segment/genotype/host/country filters are active, this uses
    cluster_filter_counts and returns matching sequence counts.

    Important
    ---------
    n_sequences means:
    - full cluster size when no metadata filter is active;
    - matching sequence count when metadata filters are active.

    full_cluster_size always represents the complete cluster size.
    """
    if limit is not None:
        limit = int(limit)

    offset = int(offset)

    if min_n_sequences is not None:
        min_n_sequences = int(min_n_sequences)

    if max_n_sequences is not None:
        max_n_sequences = int(max_n_sequences)

    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    ha_subtype_filter = _normalise_filter_value(ha_subtype_filter)
    na_subtype_filter = _normalise_filter_value(na_subtype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    metadata_filters_active = _metadata_filters_active(
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter
    )

    conn = sqlite3.connect(db_path)

    params = []

    if metadata_filters_active:
        query = """
            SELECT
                fc.organism,
                fc.cluster_number,
                cs.centroid,

                CASE
                    WHEN COUNT(DISTINCT fc.segment) = 1 THEN MAX(fc.segment)
                    WHEN COUNT(DISTINCT fc.segment) = 0 THEN NULL
                    ELSE 'Mixed'
                END AS segment,

                COUNT(DISTINCT fc.segment) AS n_segments,
                SUM(fc.n_sequences) AS n_sequences,

                cs.n_hosts,
                cs.n_countries,
                cs.n_genotypes,

                cs.n_sequences AS full_cluster_size
            FROM cluster_filter_counts fc
            LEFT JOIN cluster_summary cs
                ON fc.organism = cs.organism
                AND fc.cluster_number = cs.cluster_number
            WHERE 1 = 1
        """

        if organism_search:
            query += " AND fc.organism LIKE ?"
            params.append(f"%{organism_search}%")

        if cluster_search:
            query += " AND CAST(fc.cluster_number AS TEXT) LIKE ?"
            params.append(f"%{cluster_search}%")

        if centroid_search:
            query += " AND fc.centroid LIKE ?"
            params.append(f"%{centroid_search}%")

        query, params = _append_filter_count_filters(
            query=query,
            params=params,
            table_alias="fc",
            segment_filter=segment_filter,
            genotype_filter=genotype_filter,
            ha_subtype_filter=ha_subtype_filter,
            na_subtype_filter=na_subtype_filter,
            host_filter=host_filter,
            country_filter=country_filter
        )

        query += """
            GROUP BY
                fc.organism,
                fc.cluster_number,
                cs.centroid,
                cs.n_hosts,
                cs.n_countries,
                cs.n_genotypes,
                cs.n_sequences
            HAVING 1 = 1
        """

        if min_n_sequences is not None:
            query += " AND SUM(fc.n_sequences) >= ?"
            params.append(min_n_sequences)

        if max_n_sequences is not None:
            query += " AND SUM(fc.n_sequences) <= ?"
            params.append(max_n_sequences)

        query += """
            ORDER BY
                n_sequences DESC,
                fc.organism,
                segment,
                fc.cluster_number
        """

    else:
        query = """
            SELECT
                cs.organism,
                cs.cluster_number,
                cs.centroid,
                cs.segment,
                cs.n_segments,
                cs.n_sequences,
                cs.n_hosts,
                cs.n_countries,
                cs.n_genotypes,
                cs.n_sequences AS full_cluster_size
            FROM cluster_summary cs
            WHERE 1 = 1
        """

        if organism_search:
            query += " AND cs.organism LIKE ?"
            params.append(f"%{organism_search}%")

        if cluster_search:
            query += " AND CAST(cs.cluster_number AS TEXT) LIKE ?"
            params.append(f"%{cluster_search}%")

        if centroid_search:
            query += " AND cs.centroid LIKE ?"
            params.append(f"%{centroid_search}%")

        query += """
            GROUP BY
                cs.organism,
                cs.cluster_number,
                cs.centroid,
                cs.segment,
                cs.n_segments,
                cs.n_sequences,
                cs.n_hosts,
                cs.n_countries,
                cs.n_genotypes
            HAVING 1 = 1
        """

        if min_n_sequences is not None:
            query += " AND cs.n_sequences >= ?"
            params.append(min_n_sequences)

        if max_n_sequences is not None:
            query += " AND cs.n_sequences <= ?"
            params.append(max_n_sequences)

        query += """
            ORDER BY
                cs.n_sequences DESC,
                cs.organism,
                cs.segment,
                cs.cluster_number
        """

    if limit is not None:
        query += " LIMIT ? OFFSET ?"
        params.append(limit)
        params.append(offset)
    elif offset > 0:
        query += " LIMIT -1 OFFSET ?"
        params.append(offset)

    query += ";"

    df = pd.read_sql_query(query, conn, params=params)
    conn.close()

    return df


@st.cache_data(show_spinner=False)
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

    Returns an integer count.
    """
    if min_n_sequences is not None:
        min_n_sequences = int(min_n_sequences)

    if max_n_sequences is not None:
        max_n_sequences = int(max_n_sequences)

    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    ha_subtype_filter = _normalise_filter_value(ha_subtype_filter)
    na_subtype_filter = _normalise_filter_value(na_subtype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    metadata_filters_active = _metadata_filters_active(
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter
    )

    conn = sqlite3.connect(db_path)

    params = []

    if metadata_filters_active:
        query = """
            SELECT COUNT(*) FROM (
                SELECT 1
                FROM cluster_filter_counts fc
                LEFT JOIN cluster_summary cs
                    ON fc.organism = cs.organism
                    AND fc.cluster_number = cs.cluster_number
                WHERE 1 = 1
        """

        if organism_search:
            query += " AND fc.organism LIKE ?"
            params.append(f"%{organism_search}%")

        if cluster_search:
            query += " AND CAST(fc.cluster_number AS TEXT) LIKE ?"
            params.append(f"%{cluster_search}%")

        if centroid_search:
            query += " AND fc.centroid LIKE ?"
            params.append(f"%{centroid_search}%")

        query, params = _append_filter_count_filters(
            query=query,
            params=params,
            table_alias="fc",
            segment_filter=segment_filter,
            genotype_filter=genotype_filter,
            ha_subtype_filter=ha_subtype_filter,
            na_subtype_filter=na_subtype_filter,
            host_filter=host_filter,
            country_filter=country_filter
        )

        query += """
            GROUP BY
                fc.organism,
                fc.cluster_number,
                cs.centroid,
                cs.n_hosts,
                cs.n_countries,
                cs.n_genotypes,
                cs.n_sequences
            HAVING 1 = 1
        """

        if min_n_sequences is not None:
            query += " AND SUM(fc.n_sequences) >= ?"
            params.append(min_n_sequences)

        if max_n_sequences is not None:
            query += " AND SUM(fc.n_sequences) <= ?"
            params.append(max_n_sequences)

        query += " )"

    else:
        query = """
            SELECT COUNT(*) FROM (
                SELECT 1
                FROM cluster_summary cs
                WHERE 1 = 1
        """

        if organism_search:
            query += " AND cs.organism LIKE ?"
            params.append(f"%{organism_search}%")

        if cluster_search:
            query += " AND CAST(cs.cluster_number AS TEXT) LIKE ?"
            params.append(f"%{cluster_search}%")

        if centroid_search:
            query += " AND cs.centroid LIKE ?"
            params.append(f"%{centroid_search}%")

        query += """
            GROUP BY
                cs.organism,
                cs.cluster_number,
                cs.centroid,
                cs.segment,
                cs.n_segments,
                cs.n_sequences,
                cs.n_hosts,
                cs.n_countries,
                cs.n_genotypes
            HAVING 1 = 1
        """

        if min_n_sequences is not None:
            query += " AND cs.n_sequences >= ?"
            params.append(min_n_sequences)

        if max_n_sequences is not None:
            query += " AND cs.n_sequences <= ?"
            params.append(max_n_sequences)

        query += " )"

    query += ";"

    result = pd.read_sql_query(query, conn, params=params)
    conn.close()

    return int(result.iloc[0, 0])


# ============================================================
# SUBSET-BASED CLUSTER QUERIES
# ============================================================

@st.cache_data(show_spinner=False)
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

    Metadata filters are applied at sequence level.
    Uses cluster_member_metadata.
    """
    output_columns = [
        "organism",
        "cluster_number",
        "accession",
        "description",
        "sequence"
    ]

    clean_pairs = _clean_cluster_pairs(cluster_pairs)

    if not clean_pairs:
        return pd.DataFrame(columns=output_columns)

    chunk_size = 400

    dfs = []
    conn = sqlite3.connect(db_path)

    for i in range(0, len(clean_pairs), chunk_size):
        chunk_pairs = clean_pairs[i:i + chunk_size]

        where_clause, params = _build_cluster_pair_where_clause(
            table_alias="cmm",
            clean_pairs=chunk_pairs
        )

        query = f"""
            SELECT
                cmm.organism,
                cmm.cluster_number,
                cmm.accession,
                cmm.description,
                cmm.sequence
            FROM cluster_member_metadata cmm
            WHERE ({where_clause})
        """

        query, params = _append_metadata_filters(
            query=query,
            params=params,
            table_alias="cmm",
            segment_filter=segment_filter,
            genotype_filter=genotype_filter,
            ha_subtype_filter=ha_subtype_filter,
            na_subtype_filter=na_subtype_filter,
            host_filter=host_filter,
            country_filter=country_filter
        )

        query += """
            ORDER BY
                cmm.organism,
                cmm.cluster_number,
                cmm.accession;
        """

        chunk_df = pd.read_sql_query(query, conn, params=params)
        dfs.append(chunk_df)

    conn.close()

    if not dfs:
        return pd.DataFrame(columns=output_columns)

    return pd.concat(dfs, ignore_index=True)


@st.cache_data(show_spinner=False)
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
    Load cluster-level summary statistics for an explicit subset of clusters.

    Uses precomputed cluster_summary and cluster_filter_counts.

    If metadata filters are active:
    - n_sequences = matching sequence count
    - full_cluster_size = complete cluster size
    - length statistics remain complete-cluster statistics
    """
    clean_pairs = _clean_cluster_pairs(cluster_pairs)

    if not clean_pairs:
        return _empty_cluster_subset_summary_df()

    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    metadata_filters_active = _metadata_filters_active(
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        host_filter=host_filter,
        country_filter=country_filter
    )

    chunk_size = 400

    dfs = []
    conn = sqlite3.connect(db_path)

    for i in range(0, len(clean_pairs), chunk_size):
        chunk_pairs = clean_pairs[i:i + chunk_size]

        if metadata_filters_active:
            where_clause, params = _build_cluster_pair_where_clause(
                table_alias="fc",
                clean_pairs=chunk_pairs
            )

            query = f"""
                SELECT
                    fc.organism,
                    fc.cluster_number,
                    cs.centroid,

                    CASE
                        WHEN COUNT(DISTINCT fc.segment) = 1 THEN MAX(fc.segment)
                        WHEN COUNT(DISTINCT fc.segment) = 0 THEN NULL
                        ELSE 'Mixed'
                    END AS segment,

                    COUNT(DISTINCT fc.segment) AS n_segments,
                    SUM(fc.n_sequences) AS n_sequences,
                    cs.n_sequences AS full_cluster_size,

                    cs.n_hosts,
                    cs.n_countries,
                    cs.n_genotypes,

                    cs.length_min,
                    cs.length_q1,
                    cs.length_mean,
                    cs.length_median,
                    cs.length_q3,
                    cs.length_max,
                    cs.length_std
                FROM cluster_filter_counts fc
                LEFT JOIN cluster_summary cs
                    ON fc.organism = cs.organism
                    AND fc.cluster_number = cs.cluster_number
                WHERE ({where_clause})
            """

            query, params = _append_filter_count_filters(
                query=query,
                params=params,
                table_alias="fc",
                segment_filter=segment_filter,
                genotype_filter=genotype_filter,
                ha_subtype_filter=ha_subtype_filter,
                na_subtype_filter=na_subtype_filter,
                host_filter=host_filter,
                country_filter=country_filter
            )

            query += """
                GROUP BY
                    fc.organism,
                    fc.cluster_number,
                    cs.centroid,
                    cs.n_sequences,
                    cs.n_hosts,
                    cs.n_countries,
                    cs.n_genotypes,
                    cs.length_min,
                    cs.length_q1,
                    cs.length_mean,
                    cs.length_median,
                    cs.length_q3,
                    cs.length_max,
                    cs.length_std
                ORDER BY
                    n_sequences ASC,
                    fc.organism,
                    segment,
                    fc.cluster_number;
            """

        else:
            where_clause, params = _build_cluster_pair_where_clause(
                table_alias="cs",
                clean_pairs=chunk_pairs
            )

            query = f"""
                SELECT
                    cs.organism,
                    cs.cluster_number,
                    cs.centroid,
                    cs.segment,
                    cs.n_segments,
                    cs.n_sequences,
                    cs.n_sequences AS full_cluster_size,
                    cs.n_hosts,
                    cs.n_countries,
                    cs.n_genotypes,
                    cs.length_min,
                    cs.length_q1,
                    cs.length_mean,
                    cs.length_median,
                    cs.length_q3,
                    cs.length_max,
                    cs.length_std
                FROM cluster_summary cs
                WHERE ({where_clause})
                ORDER BY
                    cs.n_sequences ASC,
                    cs.organism,
                    cs.segment,
                    cs.cluster_number;
            """

        chunk_df = pd.read_sql_query(query, conn, params=params)
        dfs.append(chunk_df)

    conn.close()

    if not dfs:
        return _empty_cluster_subset_summary_df()

    df = pd.concat(dfs, ignore_index=True)

    if df.empty:
        return _empty_cluster_subset_summary_df()

    return df


@st.cache_data(show_spinner=False)
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
    Count background values for an explicit subset of clusters.

    Uses precomputed cluster_background_counts when only segment filtering
    is active.

    If genotype/host/country filters are active, uses cluster_member_metadata
    to respect the full filter combination exactly.
    """
    allowed_backgrounds = ["host", "country", "genotype"]

    if background_type != "All" and background_type not in allowed_backgrounds:
        raise ValueError(f"Invalid background_type: {background_type}")

    output_columns = [
        "organism",
        "cluster_number",
        "centroid",
        "segment",
        "background_type",
        "background_value",
        "n_sequences",
    ]

    clean_pairs = _clean_cluster_pairs(cluster_pairs)

    if not clean_pairs:
        return pd.DataFrame(columns=output_columns)

    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    backgrounds_to_load = (
        allowed_backgrounds
        if background_type == "All"
        else [background_type]
    )

    non_segment_filters_active = any([
        genotype_filter != "All",
        ha_subtype_filter != "All",
        na_subtype_filter != "All",
        host_filter != "All",
        country_filter != "All",
    ])

    chunk_size = 400

    dfs = []
    conn = sqlite3.connect(db_path)

    for i in range(0, len(clean_pairs), chunk_size):
        chunk_pairs = clean_pairs[i:i + chunk_size]

        if not non_segment_filters_active:
            where_clause, params = _build_cluster_pair_where_clause(
                table_alias="cbc",
                clean_pairs=chunk_pairs
            )

            background_placeholders = ",".join(["?"] * len(backgrounds_to_load))

            query = f"""
                SELECT
                    cbc.organism,
                    cbc.cluster_number,
                    cbc.centroid,
                    cbc.segment,
                    cbc.background_type,
                    cbc.background_value,
                    cbc.n_sequences
                FROM cluster_background_counts cbc
                WHERE ({where_clause})
                  AND cbc.background_type IN ({background_placeholders})
            """

            params.extend(backgrounds_to_load)

            if segment_filter != "All":
                query += " AND cbc.segment = ?"
                params.append(segment_filter)

            query += """
                ORDER BY
                    cbc.organism,
                    cbc.segment,
                    cbc.cluster_number,
                    cbc.background_type,
                    cbc.n_sequences DESC,
                    cbc.background_value;
            """

            chunk_df = pd.read_sql_query(query, conn, params=params)
            dfs.append(chunk_df)

        else:
            where_clause, params = _build_cluster_pair_where_clause(
                table_alias="cmm",
                clean_pairs=chunk_pairs
            )

            query = f"""
                SELECT
                    cmm.organism,
                    cmm.cluster_number,
                    cmm.centroid,
                    cmm.accession,
                    cmm.segment,
                    cmm.host,
                    cmm.country,
                    cmm.genotype
                FROM cluster_member_metadata cmm
                WHERE ({where_clause})
            """

            query, params = _append_metadata_filters(
                query=query,
                params=params,
                table_alias="cmm",
                segment_filter=segment_filter,
                genotype_filter=genotype_filter,
                ha_subtype_filter=ha_subtype_filter,
                na_subtype_filter=na_subtype_filter,
                host_filter=host_filter,
                country_filter=country_filter
            )

            chunk_df = pd.read_sql_query(query, conn, params=params)

            if chunk_df.empty:
                continue

            output_dfs = []

            for field in backgrounds_to_load:
                field_df = chunk_df[
                    chunk_df[field].notna()
                    & (chunk_df[field].astype(str).str.strip() != "")
                ].copy()

                if field_df.empty:
                    continue

                grouped = (
                    field_df.groupby(
                        [
                            "organism",
                            "cluster_number",
                            "centroid",
                            "segment",
                            field,
                        ],
                        dropna=False
                    )
                    .agg(
                        n_sequences=("accession", "nunique")
                    )
                    .reset_index()
                    .rename(columns={field: "background_value"})
                )

                grouped["background_type"] = field

                grouped = grouped[
                    [
                        "organism",
                        "cluster_number",
                        "centroid",
                        "segment",
                        "background_type",
                        "background_value",
                        "n_sequences",
                    ]
                ]

                output_dfs.append(grouped)

            if output_dfs:
                dfs.append(pd.concat(output_dfs, ignore_index=True))

    conn.close()

    if not dfs:
        return pd.DataFrame(columns=output_columns)

    result = pd.concat(dfs, ignore_index=True)

    if result.empty:
        return pd.DataFrame(columns=output_columns)

    result = result.sort_values(
        by=[
            "organism",
            "segment",
            "cluster_number",
            "background_type",
            "n_sequences",
            "background_value",
        ],
        ascending=[True, True, True, True, False, True]
    )

    return result


@st.cache_data(show_spinner=False)
def load_cluster_length_stats_for_cluster_subset(
    db_path,
    cluster_pairs,
    segment_filter="All",
    genotype_filter="All",
    host_filter="All",
    country_filter="All"
):
    """
    Load sequence length statistics for an explicit subset of clusters.

    Uses precomputed cluster_summary.

    Important
    ---------
    Length statistics are complete-cluster statistics. If metadata filters
    are active, n_sequences is the matching sequence count, while
    full_cluster_size remains the complete cluster size.
    """
    summary_df = load_cluster_summary_for_cluster_subset(
        db_path=db_path,
        cluster_pairs=cluster_pairs,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        host_filter=host_filter,
        country_filter=country_filter
    )

    output_columns = [
        "organism",
        "cluster_number",
        "centroid",
        "segment",
        "n_sequences",
        "full_cluster_size",
        "length_min",
        "length_q1",
        "length_mean",
        "length_median",
        "length_q3",
        "length_max",
        "length_std",
    ]

    if summary_df.empty:
        return pd.DataFrame(columns=output_columns)

    return summary_df[output_columns].copy()


@st.cache_data(show_spinner=False)
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
    Load a cluster-level descriptor table for an explicit subset of clusters.

    Uses precomputed cluster_descriptors.

    For compatibility with the existing app, JSON columns are aliased as:
    - genotypes
    - segments
    - hosts
    - countries
    """
    output_columns = [
        "organism",
        "cluster_number",
        "representative",
        "n_sequences",
        "genotypes",
        "segments",
        "hosts",
        "countries",
        "dominant_genotype",
        "dominant_segment",
        "dominant_host",
        "dominant_country",
        "n_genotypes",
        "n_segments",
        "n_hosts",
        "n_countries",
        "is_singleton",
        "is_mixed_genotype",
        "is_mixed_segment",
        "is_multi_host",
        "is_multi_country",
        "length_min",
        "length_mean",
        "length_median",
        "length_max",
        "length_std",
    ]

    clean_pairs = _clean_cluster_pairs(cluster_pairs)

    if not clean_pairs:
        return _empty_cluster_descriptor_df()

    segment_filter = _normalise_filter_value(segment_filter)
    genotype_filter = _normalise_filter_value(genotype_filter)
    ha_subtype_filter = _normalise_filter_value(ha_subtype_filter)
    na_subtype_filter = _normalise_filter_value(na_subtype_filter)
    host_filter = _normalise_filter_value(host_filter)
    country_filter = _normalise_filter_value(country_filter)

    metadata_filters_active = _metadata_filters_active(
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter
    )

    chunk_size = 400

    dfs = []
    conn = sqlite3.connect(db_path)

    for i in range(0, len(clean_pairs), chunk_size):
        chunk_pairs = clean_pairs[i:i + chunk_size]

        if metadata_filters_active:
            where_clause, params = _build_cluster_pair_where_clause(
                table_alias="fc",
                clean_pairs=chunk_pairs
            )

            query = f"""
                SELECT DISTINCT
                    cd.organism,
                    cd.cluster_number,
                    cd.representative,
                    cd.n_sequences,
                    cd.genotypes_json AS genotypes,
                    cd.segments_json AS segments,
                    cd.hosts_json AS hosts,
                    cd.countries_json AS countries,
                    cd.dominant_genotype,
                    cd.dominant_segment,
                    cd.dominant_host,
                    cd.dominant_country,
                    cd.n_genotypes,
                    cd.n_segments,
                    cd.n_hosts,
                    cd.n_countries,
                    cd.is_singleton,
                    cd.is_mixed_genotype,
                    cd.is_mixed_segment,
                    cd.is_multi_host,
                    cd.is_multi_country,
                    cd.length_min,
                    cd.length_mean,
                    cd.length_median,
                    cd.length_max,
                    cd.length_std
                FROM cluster_filter_counts fc
                LEFT JOIN cluster_descriptors cd
                    ON fc.organism = cd.organism
                    AND fc.cluster_number = cd.cluster_number
                WHERE ({where_clause})
            """

            query, params = _append_filter_count_filters(
                query=query,
                params=params,
                table_alias="fc",
                segment_filter=segment_filter,
                genotype_filter=genotype_filter,
                ha_subtype_filter=ha_subtype_filter,
                na_subtype_filter=na_subtype_filter,
                host_filter=host_filter,
                country_filter=country_filter
            )

            query += """
                ORDER BY
                    cd.n_sequences ASC,
                    cd.organism,
                    cd.dominant_segment,
                    cd.cluster_number;
            """

        else:
            where_clause, params = _build_cluster_pair_where_clause(
                table_alias="cd",
                clean_pairs=chunk_pairs
            )

            query = f"""
                SELECT
                    cd.organism,
                    cd.cluster_number,
                    cd.representative,
                    cd.n_sequences,
                    cd.genotypes_json AS genotypes,
                    cd.segments_json AS segments,
                    cd.hosts_json AS hosts,
                    cd.countries_json AS countries,
                    cd.dominant_genotype,
                    cd.dominant_segment,
                    cd.dominant_host,
                    cd.dominant_country,
                    cd.n_genotypes,
                    cd.n_segments,
                    cd.n_hosts,
                    cd.n_countries,
                    cd.is_singleton,
                    cd.is_mixed_genotype,
                    cd.is_mixed_segment,
                    cd.is_multi_host,
                    cd.is_multi_country,
                    cd.length_min,
                    cd.length_mean,
                    cd.length_median,
                    cd.length_max,
                    cd.length_std
                FROM cluster_descriptors cd
                WHERE ({where_clause})
                ORDER BY
                    cd.n_sequences ASC,
                    cd.organism,
                    cd.dominant_segment,
                    cd.cluster_number;
            """

        chunk_df = pd.read_sql_query(query, conn, params=params)
        dfs.append(chunk_df)

    conn.close()

    if not dfs:
        return _empty_cluster_descriptor_df()

    df = pd.concat(dfs, ignore_index=True)

    if df.empty:
        return _empty_cluster_descriptor_df()

    return df[output_columns]