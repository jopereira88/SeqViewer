#!/usr/bin/env python3
"""
Streamlit interface for browsing sequence metadata, sequences and clusters.

Run with:

    streamlit run app.py
"""

from pathlib import Path
import json

import pandas as pd
import streamlit as st

from db import connect_sqlite
from services import (
    get_table_counts,
    get_distinct_values,
    count_sequences,
    load_sequence_index,
    load_sequence_detail,
    load_sequences_by_accessions,
    load_sequence_records_by_accessions,
    load_cluster_members,
    load_cluster_member_sequences,
    load_cluster_metadata_members_for_cluster_subset,
    load_cluster_metadata_clusters_for_cluster_subset,
    iter_cluster_member_sequences_for_cluster_subset,
    iter_cluster_metadata_members_for_cluster_subset,
    search_clusters,
    count_clusters,
    get_assembly_links,
    load_cluster_assembly_values,
)
from utils import (
    format_fasta,
    format_multifasta,
    build_cluster_zip,
    build_cluster_zip_from_rows,
    cluster_metadata_zip_bytes_from_rows,
    parse_accession_list,
)


# ============================================================
# CONFIG
# ============================================================

DEFAULT_DB_PATH = "data/sequence_database.sqlite"

st.set_page_config(
    page_title="Sequence Viewer",
    layout="wide"
)


# ============================================================
# HELPERS
# ============================================================

def _sequences_have_ambiguous_counts(db_path):
    """
    Return True if the sequences table has the ETL-precomputed ambiguous
    base count columns.
    """
    import sqlite3

    conn = sqlite3.connect(db_path)

    try:
        columns = {
            row[1]
            for row in conn.execute("PRAGMA table_info(sequences)").fetchall()
        }
    finally:
        conn.close()

    return {"n_N", "n_degenerate"} <= columns


def _render_sequence_detail(db_path, selected_accession, detail):
    """
    Render the full detail view for one selected sequence.

    Shows a summary column plus Metadata, Sequence, Ambiguous bases,
    Cluster and Cluster members tabs.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    selected_accession : str
        Accession of the selected sequence.

    detail : dict
        Full detail record from load_sequence_detail.
    """
    st.divider()
    st.subheader(f"Selected accession: {selected_accession}")

    def _int_value(key, default=0):
        value = detail.get(key)
        if value is None or pd.isna(value):
            return default
        return int(value)

    left_col, right_col = st.columns([1, 2])

    with left_col:
        st.markdown("### Summary")

        st.write(f"**Accession:** {detail.get('accession')}")
        st.write(f"**Organism:** {detail.get('organism_name')}")
        st.write(f"**Species:** {detail.get('species')}")
        st.write(f"**Segment:** {detail.get('segment')}")
        st.write(f"**Genotype:** {detail.get('genotype')}")
        st.write(f"**Length:** {detail.get('length')}")
        st.write(f"**N bases:** {_int_value('n_N')}")
        st.write(f"**Degenerate bases:** {_int_value('n_degenerate')}")
        st.write(f"**Host:** {detail.get('host')}")
        st.write(f"**Country:** {detail.get('country')}")
        st.write(
            f"**Collection date:** {detail.get('collection_date')}"
        )
        st.write(f"**Release date:** {detail.get('release_date')}")

        fasta_text = format_fasta(
            accession=detail.get("accession"),
            description=detail.get("description"),
            sequence=detail.get("sequence")
        )

        st.download_button(
            label="Download selected sequence as FASTA",
            data=fasta_text,
            file_name=f"{selected_accession}.fasta",
            mime="text/plain",
            disabled=(not bool(fasta_text))
        )

    with right_col:
        tabs = st.tabs([
            "Metadata",
            "Sequence",
            "Ambiguous bases",
            "Cluster",
            "Cluster members"
        ])

        with tabs[0]:
            st.markdown("### Metadata")

            metadata_fields = {
                key: value
                for key, value in detail.items()
                if key not in {
                    "sequence",
                    "n_N",
                    "n_degenerate",
                    "degenerate_breakdown",
                    "cluster_organism",
                    "cluster_number",
                    "identity_to_centroid",
                    "centroid",
                }
            }

            st.json(metadata_fields)

        with tabs[1]:
            st.markdown("### Sequence")

            sequence = detail.get("sequence")

            if sequence is None:
                st.warning("No sequence found for this accession.")
            else:
                st.write(
                    f"Sequence length in database: {len(sequence)}"
                )
                st.text_area(
                    label="Sequence",
                    value=sequence,
                    height=300
                )

        with tabs[2]:
            st.markdown("### Ambiguous bases")

            sequence = detail.get("sequence")
            n_n = _int_value("n_N")
            n_degenerate = _int_value("n_degenerate")

            raw_breakdown = detail.get("degenerate_breakdown")

            if isinstance(raw_breakdown, str):
                try:
                    per_base = json.loads(raw_breakdown)
                except (ValueError, TypeError):
                    per_base = {}
            else:
                per_base = dict(raw_breakdown or {})

            if sequence is not None:
                st.write(f"**Sequence length:** {len(sequence)}")
            st.write(f"**N bases:** {n_n}")
            st.write(f"**Degenerate bases:** {n_degenerate}")

            breakdown_rows = [
                {"Base": base, "Count": count}
                for base, count in per_base.items()
            ]
            breakdown_rows.append({"Base": "N", "Count": n_n})
            breakdown_rows.append(
                {"Base": "Total", "Count": n_degenerate + n_n}
            )

            st.dataframe(
                pd.DataFrame(breakdown_rows),
                use_container_width=True,
                hide_index=True,
            )

        with tabs[3]:
            st.markdown("### Cluster")

            cluster_organism = detail.get("cluster_organism")
            cluster_number = detail.get("cluster_number")
            centroid = detail.get("centroid")
            identity = detail.get("identity_to_centroid")

            if cluster_number is None or pd.isna(cluster_number):
                st.info(
                    "No cluster information available for this accession."
                )
            else:
                st.write(
                    f"**Cluster organism:** {cluster_organism}"
                )
                st.write(
                    f"**Cluster number:** {cluster_number}"
                )
                st.write(f"**Centroid:** {centroid}")
                st.write(
                    f"**Identity to centroid:** {identity}"
                )

        with tabs[4]:
            st.markdown("### Cluster members")

            cluster_organism = detail.get("cluster_organism")
            cluster_number = detail.get("cluster_number")

            if cluster_number is None or pd.isna(cluster_number):
                st.info("No cluster selected.")
            else:
                members_df = load_cluster_members(
                    db_path=db_path,
                    organism=cluster_organism,
                    cluster_number=cluster_number
                )

                st.write(f"**Cluster:** {cluster_number}")
                st.write(f"**Members:** {len(members_df)}")

                st.dataframe(
                    members_df,
                    use_container_width=True,
                    hide_index=True
                )


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("Sequence Viewer")

db_path = st.sidebar.text_input(
    "SQLite database path",
    value=DEFAULT_DB_PATH
)

if not Path(db_path).exists():
    st.error(f"Database not found: {db_path}")
    st.stop()

if not _sequences_have_ambiguous_counts(db_path):
    st.error(
        "The sequences table is missing the precomputed ambiguous-base "
        "columns (n_N, n_degenerate). Re-import the database with the ETL "
        "using the updated schema: "
        "`.viewer/bin/python db/import_data.py ... --recreate`"
    )
    st.stop()

counts = get_table_counts(db_path)

st.sidebar.subheader("Database summary")

for table_name, count in counts.items():
    st.sidebar.write(f"**{table_name}:** {count}")

st.sidebar.divider()

# -- Metadata filters (shared by both tabs) --

st.sidebar.subheader("Metadata filters")

segment_filter = st.sidebar.selectbox(
    "Segment",
    options=get_distinct_values(db_path, "metadata", "segment")
)

genotype_filter = st.sidebar.selectbox(
    "Genotype",
    options=get_distinct_values(db_path, "metadata", "genotype")
)

ha_subtype_filter = st.sidebar.selectbox(
    "HA Subtype",
    options=get_distinct_values(db_path, "metadata", "ha_subtype")
)

na_subtype_filter = st.sidebar.selectbox(
    "NA Subtype",
    options=get_distinct_values(db_path, "metadata", "na_subtype")
)

host_filter = st.sidebar.selectbox(
    "Host",
    options=get_distinct_values(db_path, "metadata", "host")
)

country_filter = st.sidebar.selectbox(
    "Country",
    options=get_distinct_values(db_path, "metadata", "country")
)

st.sidebar.divider()

# -- Sequence search filters --

st.sidebar.subheader("Sequence search")

multiple_accessions = st.sidebar.checkbox(
    "Search multiple accessions (semicolon-separated)",
    value=False,
    help=(
        "Paste a list of accessions separated by semicolons (;). "
        "Commas, newlines and whitespace are also accepted."
    )
)

if multiple_accessions:
    accession_list_text = st.sidebar.text_area(
        "Accession list",
        value="",
        placeholder="CY123456.1; CY123457.1; ...",
        height=140,
    )
    accession_search = ""
else:
    accession_list_text = ""
    accession_search = st.sidebar.text_input(
        "Search accession / description / organism",
        value=""
    )

organism_search = st.sidebar.text_input(
    "Search organism name",
    value=""
)

species_search = st.sidebar.text_input(
    "Search species",
    value=""
)

centroids_only = st.sidebar.checkbox(
    "Cluster centroids only",
    value=False,
    help="Only show sequences that are cluster centroids."
)

SEQ_PAGE_SIZES = (
    list(range(10, 100, 10))
    + list(range(100, 1001, 50))
)

seq_page_size = st.sidebar.selectbox(
    "Results per page",
    options=SEQ_PAGE_SIZES,
    index=SEQ_PAGE_SIZES.index(100),
)

selection_mode = st.sidebar.radio(
    "Selection mode",
    options=["Detail view", "Download (multi-select)"],
    help=(
        "Detail view: click a row to inspect metadata, sequence and "
        "cluster info. Download: select multiple rows to download a "
        "multifasta file."
    ),
)

st.sidebar.divider()

# -- Degenerate base / N content filters --

st.sidebar.subheader("Degenerate bases / N content")

min_degenerate = st.sidebar.number_input(
    "Min degenerate bases",
    min_value=0,
    max_value=30000,
    value=0,
    step=1,
    help="Minimum total of R, Y, S, W, K, M, B, D, H, V per sequence."
)

max_degenerate = st.sidebar.number_input(
    "Max degenerate bases",
    min_value=0,
    max_value=30000,
    value=30000,
    step=1,
    help="Maximum total of R, Y, S, W, K, M, B, D, H, V per sequence."
)

min_n = st.sidebar.number_input(
    "Min N bases",
    min_value=0,
    max_value=30000,
    value=0,
    step=1,
    help="Minimum number of N bases per sequence."
)

max_n = st.sidebar.number_input(
    "Max N bases",
    min_value=0,
    max_value=30000,
    value=30000,
    step=1,
    help="Maximum number of N bases per sequence."
)

st.sidebar.divider()

# -- Cluster search filters --

st.sidebar.subheader("Cluster search")

cluster_organism_options = get_distinct_values(db_path, "clusters", "organism")

cluster_organism_filter = st.sidebar.selectbox(
    "Cluster organism",
    options=cluster_organism_options
)

min_cluster_size = st.sidebar.number_input(
    "Minimum cluster size",
    min_value=1,
    max_value=10000,
    value=2,
    step=1,
    help="Hide singletons by default. Set to 1 to show all."
)


# ============================================================
# MAIN PAGE
# ============================================================

st.title("Sequence Viewer")

st.caption(
    "Interactive browser for sequence metadata, FASTA sequences and cluster membership."
)

tab_sequences, tab_clusters = st.tabs(["Sequences", "Clusters"])


# ============================================================
# SEQUENCES TAB
# ============================================================

with tab_sequences:

    parsed_accessions = (
        parse_accession_list(accession_list_text)
        if multiple_accessions and accession_list_text.strip()
        else []
    )

    if multiple_accessions and parsed_accessions:
        # ============================================================
        # MULTI-ACCESSION SEARCH
        # ============================================================

        st.subheader("Accession list search")

        acc_cache_key = (
            db_path,
            tuple(sorted(parsed_accessions)),
            segment_filter,
            genotype_filter,
            ha_subtype_filter,
            na_subtype_filter,
            host_filter,
            country_filter,
            min_degenerate,
            max_degenerate,
            min_n,
            max_n,
        )

        if st.session_state.get("acc_key") != acc_cache_key:
            st.session_state["acc_key"] = acc_cache_key
            st.session_state.pop("acc_data", None)

        if "acc_data" not in st.session_state:
            with st.status(
                "Loading accession records...", expanded=False
            ) as status:
                found_df = load_sequences_by_accessions(
                    db_path=db_path,
                    accessions=parsed_accessions,
                )
                found_accessions = set(found_df["accession"].astype(str))
                not_found = [
                    accession
                    for accession in parsed_accessions
                    if accession not in found_accessions
                ]

                records_df = load_sequence_records_by_accessions(
                    db_path=db_path,
                    accessions=parsed_accessions,
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

                filtered_out = len(found_df) - len(records_df)

                st.session_state["acc_data"] = {
                    "records_df": records_df,
                    "found_count": len(found_df),
                    "not_found": not_found,
                    "filtered_out": filtered_out,
                }
                status.update(
                    label="Accession records loaded.", state="complete"
                )

        acc_data = st.session_state["acc_data"]
        records_df = acc_data["records_df"].copy()
        not_found = acc_data["not_found"]
        filtered_out = acc_data["filtered_out"]

        if not records_df.empty:
            records_df["sequence_length"] = (
                records_df["sequence"].fillna("").str.len()
            )

        total_n = int(records_df["n_N"].sum()) if not records_df.empty else 0
        total_degenerate = (
            int(records_df["n_degenerate"].sum())
            if not records_df.empty
            else 0
        )

        metric_row = st.columns(4)
        metric_row[0].metric("Requested", len(parsed_accessions))
        metric_row[1].metric("Found in DB", acc_data["found_count"])
        metric_row[2].metric("Filtered out", filtered_out)
        metric_row[3].metric("Not found", len(not_found))

        count_row = st.columns(2)
        count_row[0].metric("Total Ns (displayed)", total_n)
        count_row[1].metric("Total degenerate bases (displayed)", total_degenerate)

        if not_found:
            st.warning(
                f"{len(not_found)} accession(s) not found in the database: "
                + ", ".join(not_found[:50])
                + (" ..." if len(not_found) > 50 else "")
            )

        if records_df.empty:
            st.info(
                "No records match the selected accessions and filters "
                "(all found accessions may have been filtered out)."
            )
        else:
            display_columns = [
                "accession",
                "organism_name",
                "segment",
                "genotype",
                "host",
                "country",
                "length",
                "sequence_length",
                "n_N",
                "n_degenerate",
            ]

            st.download_button(
                label=(
                    f"Download displayed sequences as multifasta "
                    f"({len(records_df)} sequences)"
                ),
                data=format_multifasta(records_df),
                file_name="accession_list.fasta",
                mime="text/plain",
            )

            acc_event = st.dataframe(
                records_df[display_columns],
                use_container_width=True,
                hide_index=True,
                on_select="rerun",
                selection_mode="single-row",
                key="acc_results_table",
            )

            acc_selected_rows = acc_event.selection.rows

            if not acc_selected_rows:
                st.info("Select a row to view sequence details.")
            else:
                acc_selected_row = records_df.iloc[acc_selected_rows[0]]
                selected_accession = str(acc_selected_row["accession"])

                detail = load_sequence_detail(
                    db_path=db_path,
                    accession=selected_accession
                )

                if detail is None:
                    st.error(
                        f"Could not load detail for accession: "
                        f"{selected_accession}"
                    )
                else:
                    _render_sequence_detail(
                        db_path=db_path,
                        selected_accession=selected_accession,
                        detail=detail,
                    )

    else:
        # ============================================================
        # SINGLE SEARCH INDEX
        # ============================================================

        # -- Pagination state --

        if "seq_page" not in st.session_state:
            st.session_state.seq_page = 0

        # -- Count total matching sequences --

        total_sequences = count_sequences(
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

        total_pages = (
            max(1, (total_sequences + seq_page_size - 1) // seq_page_size)
            if total_sequences > 0
            else 1
        )
        current_page = min(st.session_state.seq_page, total_pages - 1)

        # -- Reset page when filters change --

        seq_filter_key = (
            db_path, accession_search, organism_search, species_search,
            segment_filter, genotype_filter, ha_subtype_filter,
            na_subtype_filter, host_filter, country_filter,
            min_degenerate, max_degenerate, min_n, max_n,
            centroids_only, seq_page_size,
        )
        if st.session_state.get("seq_filter_key") != seq_filter_key:
            current_page = 0
            st.session_state.seq_page = 0
            st.session_state["seq_filter_key"] = seq_filter_key
            st.session_state.pop("seq_detail_accession", None)

        st.session_state.seq_page = current_page
        offset = current_page * seq_page_size

        st.subheader(
            f"Sequences found: {total_sequences}"
            f" (page {current_page + 1} of {total_pages})"
        )

        # -- Page navigation --

        nav_cols = st.columns([1, 1, 1])
        with nav_cols[0]:
            if st.button("Previous", disabled=(current_page == 0), key="seq_prev"):
                st.session_state.seq_page -= 1
                st.rerun()
        with nav_cols[2]:
            if st.button(
                "Next",
                disabled=(current_page >= total_pages - 1),
                key="seq_next",
            ):
                st.session_state.seq_page += 1
                st.rerun()

        # -- Load current page only --

        is_download_mode = selection_mode == "Download (multi-select)"

        sequence_index_df = load_sequence_index(
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
            limit=seq_page_size,
            offset=offset,
            centroids_only=centroids_only,
        )

        if sequence_index_df.empty:
            st.warning("No records found for the selected filters.")
        else:
            event = st.dataframe(
                sequence_index_df,
                use_container_width=True,
                hide_index=True,
                on_select="rerun",
                selection_mode="multi-row" if is_download_mode else "single-row",
            )

            selected_rows = list(event.selection.rows)

            # Preserve the selected accession across "Selection mode" toggles.
            # The dataframe widget id depends on selection_mode, so Streamlit
            # discards the row selection (and hides the detail panel) when the
            # mode changes. Restore it here when we are back in detail view.
            if not is_download_mode:
                if selected_rows:
                    st.session_state["seq_detail_accession"] = str(
                        sequence_index_df.iloc[selected_rows[0]]["accession"]
                    )
                elif "seq_detail_accession" in st.session_state:
                    restored_idx = sequence_index_df.index[
                        sequence_index_df["accession"].astype(str)
                        == st.session_state["seq_detail_accession"]
                    ]
                    if len(restored_idx):
                        selected_rows = [restored_idx[0]]
                    else:
                        st.session_state.pop("seq_detail_accession", None)

            if is_download_mode:
                # -- Multi-select download --

                download_accessions = (
                    sequence_index_df.iloc[selected_rows]["accession"]
                    .astype(str)
                    .tolist()
                    if selected_rows
                    else sequence_index_df["accession"].astype(str).tolist()
                )

                download_label = (
                    f"Download selected sequences as multifasta "
                    f"({len(download_accessions)} sequences)"
                    if selected_rows
                    else f"Download all sequences on this page as multifasta "
                         f"({len(download_accessions)} sequences)"
                )

                with st.status(
                    "Loading sequences for download...", expanded=False
                ) as status:
                    download_seqs = load_sequences_by_accessions(
                        db_path=db_path,
                        accessions=download_accessions,
                    )
                    status.update(
                        label="Sequences loaded.", state="complete"
                    )

                if download_seqs.empty:
                    st.warning("No sequence data found for selected rows.")
                else:
                    st.download_button(
                        label=download_label,
                        data=format_multifasta(download_seqs),
                        file_name="selected_sequences.fasta",
                        mime="text/plain",
                    )
            else:
                # -- Single-row detail view --

                if not selected_rows:
                    st.info("Select one row from the sequence index.")
                else:
                    selected_idx = selected_rows[0]
                    selected_row = sequence_index_df.iloc[selected_idx]

                    selected_accession = selected_row["accession"]

                    detail = load_sequence_detail(
                        db_path=db_path,
                        accession=selected_accession
                    )

                    if detail is None:
                        st.error(
                            f"Could not load detail for accession: "
                            f"{selected_accession}"
                        )
                    else:
                        _render_sequence_detail(
                            db_path=db_path,
                            selected_accession=selected_accession,
                            detail=detail,
                        )


# ============================================================
# CLUSTERS TAB
# ============================================================

with tab_clusters:

    st.subheader("Cluster Browser")

    st.caption(
        "Search clusters by organism, cluster number, centroid or assembly. "
        "Select a cluster to view members, its contributing assemblies, "
        "linked clusters and download exports."
    )

    col_cluster_search, col_centroid_search, col_assembly_search = st.columns(3)

    with col_cluster_search:
        cluster_search_term = st.text_input(
            "Search cluster number",
            value="",
            key="cluster_search_input"
        )

    with col_centroid_search:
        centroid_search_term = st.text_input(
            "Search centroid accession",
            value="",
            key="centroid_search_input"
        )

    with col_assembly_search:
        assembly_search_term = st.text_input(
            "Search assembly",
            value="",
            key="assembly_search_input",
            help="Show clusters containing sequences from assemblies "
                 "matching this accession."
        )

    # -- Pagination state --

    if "cluster_page" not in st.session_state:
        st.session_state.cluster_page = 0

    # -- Page size selector --

    page_size = st.selectbox("Results per page", options=[100, 250, 500], index=0)

    # -- Count total matching clusters --

    total_clusters = count_clusters(
        db_path=db_path,
        organism_search=cluster_organism_filter if cluster_organism_filter != "All" else "",
        cluster_search=cluster_search_term,
        centroid_search=centroid_search_term,
        segment_filter=segment_filter,
        genotype_filter=genotype_filter,
        ha_subtype_filter=ha_subtype_filter,
        na_subtype_filter=na_subtype_filter,
        host_filter=host_filter,
        country_filter=country_filter,
        min_n_sequences=min_cluster_size,
        assembly_search=assembly_search_term,
    )

    total_pages = max(1, (total_clusters + page_size - 1) // page_size) if total_clusters > 0 else 1
    current_page = min(st.session_state.cluster_page, total_pages - 1)

    # Reset page and stale downloads when filters change
    filter_key = (
        cluster_organism_filter, cluster_search_term, centroid_search_term,
        assembly_search_term,
        segment_filter, genotype_filter, ha_subtype_filter, na_subtype_filter,
        host_filter, country_filter,
        min_cluster_size, page_size
    )
    if st.session_state.get("cluster_filter_key") != filter_key:
        current_page = 0
        st.session_state.cluster_page = 0
        st.session_state["cluster_filter_key"] = filter_key
        st.session_state.pop("centroids_fasta", None)
        st.session_state.pop("centroids_count", None)
        st.session_state.pop("clusters_zip", None)
        st.session_state.pop("clusters_zip_info", None)
        st.session_state.pop("cluster_metadata_zip", None)
        st.session_state.pop("cluster_metadata_zip_info", None)

    if total_clusters == 0:
        st.warning("No clusters found for the selected filters.")
    else:
        st.session_state.cluster_page = current_page
        offset = current_page * page_size

        st.subheader(f"Clusters found: {total_clusters} (page {current_page + 1} of {total_pages})")

        # -- Page navigation --

        nav_cols = st.columns([1, 1, 1])
        with nav_cols[0]:
            if st.button("Previous", disabled=(current_page == 0), key="cluster_prev"):
                st.session_state.cluster_page -= 1
                st.rerun()
        with nav_cols[2]:
            if st.button("Next", disabled=(current_page >= total_pages - 1), key="cluster_next"):
                st.session_state.cluster_page += 1
                st.rerun()

        # -- Load current page only --

        cluster_results_df = search_clusters(
            db_path=db_path,
            organism_search=cluster_organism_filter if cluster_organism_filter != "All" else "",
            cluster_search=cluster_search_term,
            centroid_search=centroid_search_term,
            segment_filter=segment_filter,
            genotype_filter=genotype_filter,
            ha_subtype_filter=ha_subtype_filter,
            na_subtype_filter=na_subtype_filter,
            host_filter=host_filter,
            country_filter=country_filter,
            min_n_sequences=min_cluster_size,
            assembly_search=assembly_search_term,
            include_assembly_columns=True,
            limit=page_size,
            offset=offset,
        )

        # -- On-demand bulk downloads (all matching clusters, not just current page) --

        dl_cols = st.columns(3)

        with dl_cols[0]:
            if st.button("Generate centroids FASTA", key="gen_centroids_fasta"):
                all_clusters_df = search_clusters(
                    db_path=db_path,
                    organism_search=cluster_organism_filter if cluster_organism_filter != "All" else "",
                    cluster_search=cluster_search_term,
                    centroid_search=centroid_search_term,
                    segment_filter=segment_filter,
                    genotype_filter=genotype_filter,
                    ha_subtype_filter=ha_subtype_filter,
                    na_subtype_filter=na_subtype_filter,
                    host_filter=host_filter,
                    country_filter=country_filter,
                    min_n_sequences=min_cluster_size,
                    assembly_search=assembly_search_term,
                    limit=None,
                )
                centroid_accessions = all_clusters_df["centroid"].dropna().unique().tolist()
                centroid_seqs = load_sequences_by_accessions(
                    db_path=db_path,
                    accessions=centroid_accessions,
                )
                centroid_fasta = format_multifasta(centroid_seqs)
                st.session_state["centroids_fasta"] = centroid_fasta
                st.session_state["centroids_count"] = len(centroid_seqs)

            if "centroids_fasta" in st.session_state:
                st.download_button(
                    label=f"Download centroids ({st.session_state['centroids_count']} sequences)",
                    data=st.session_state["centroids_fasta"],
                    file_name="centroids.fasta",
                    mime="text/plain",
                )

        with dl_cols[1]:
            if st.button("Generate clusters ZIP", key="gen_clusters_zip"):
                all_clusters_df = search_clusters(
                    db_path=db_path,
                    organism_search=cluster_organism_filter if cluster_organism_filter != "All" else "",
                    cluster_search=cluster_search_term,
                    centroid_search=centroid_search_term,
                    segment_filter=segment_filter,
                    genotype_filter=genotype_filter,
                    ha_subtype_filter=ha_subtype_filter,
                    na_subtype_filter=na_subtype_filter,
                    host_filter=host_filter,
                    country_filter=country_filter,
                    min_n_sequences=min_cluster_size,
                    assembly_search=assembly_search_term,
                    limit=None,
                )
                all_cluster_pairs = list(
                    zip(all_clusters_df["organism"], all_clusters_df["cluster_number"])
                )
                # Streamed rather than loaded as a dataframe: a whole-dataset export is
                # over a million sequences and needs gigabytes as a dataframe.
                zip_bytes, n_sequences = build_cluster_zip_from_rows(
                    cluster_keys=all_cluster_pairs,
                    rows=iter_cluster_member_sequences_for_cluster_subset(
                        db_path=db_path,
                        cluster_pairs=all_cluster_pairs,
                        segment_filter=segment_filter,
                        genotype_filter=genotype_filter,
                        ha_subtype_filter=ha_subtype_filter,
                        na_subtype_filter=na_subtype_filter,
                        host_filter=host_filter,
                        country_filter=country_filter,
                    ),
                )
                st.session_state["clusters_zip"] = zip_bytes
                st.session_state["clusters_zip_info"] = f"{len(all_cluster_pairs)} clusters, {n_sequences} sequences"

            if "clusters_zip" in st.session_state:
                st.download_button(
                    label=f"Download all clusters as zip ({st.session_state['clusters_zip_info']})",
                    data=st.session_state["clusters_zip"],
                    file_name="clusters.zip",
                    mime="application/zip",
                )

        with dl_cols[2]:
            if st.button("Generate metadata TSV zip", key="gen_cluster_metadata_zip"):
                all_clusters_df = search_clusters(
                    db_path=db_path,
                    organism_search=cluster_organism_filter if cluster_organism_filter != "All" else "",
                    cluster_search=cluster_search_term,
                    centroid_search=centroid_search_term,
                    segment_filter=segment_filter,
                    genotype_filter=genotype_filter,
                    ha_subtype_filter=ha_subtype_filter,
                    na_subtype_filter=na_subtype_filter,
                    host_filter=host_filter,
                    country_filter=country_filter,
                    min_n_sequences=min_cluster_size,
                    assembly_search=assembly_search_term,
                    limit=None,
                )
                all_cluster_pairs = list(
                    zip(all_clusters_df["organism"], all_clusters_df["cluster_number"])
                )
                metadata_clusters_df = load_cluster_metadata_clusters_for_cluster_subset(
                    db_path=db_path,
                    cluster_pairs=all_cluster_pairs,
                    segment_filter=segment_filter,
                    genotype_filter=genotype_filter,
                    ha_subtype_filter=ha_subtype_filter,
                    na_subtype_filter=na_subtype_filter,
                    host_filter=host_filter,
                    country_filter=country_filter,
                )
                metadata_member_columns = [
                    "organism",
                    "cluster_number",
                    "centroid",
                    "accession",
                    "assembly",
                    "segment",
                    "genotype",
                    "ha_subtype",
                    "na_subtype",
                    "host",
                    "country",
                    "collection_date",
                    "length",
                    "identity_to_centroid",
                ]

                # members.tsv is streamed from a cursor rather than
                # built as a dataframe: a whole-dataset export is over a
                # million rows and the dataframe would cost about 1.5 GB.
                member_rows = iter_cluster_metadata_members_for_cluster_subset(
                    db_path=db_path,
                    cluster_pairs=all_cluster_pairs,
                    segment_filter=segment_filter,
                    genotype_filter=genotype_filter,
                    ha_subtype_filter=ha_subtype_filter,
                    na_subtype_filter=na_subtype_filter,
                    host_filter=host_filter,
                    country_filter=country_filter,
                )
                metadata_zip_bytes, n_members = cluster_metadata_zip_bytes_from_rows(
                    clusters_df=metadata_clusters_df,
                    member_columns=metadata_member_columns,
                    member_rows=member_rows,
                )
                st.session_state["cluster_metadata_zip"] = metadata_zip_bytes
                st.session_state["cluster_metadata_zip_info"] = (
                    f"{len(metadata_clusters_df)} clusters, "
                    f"{n_members} members"
                )

            if "cluster_metadata_zip" in st.session_state:
                st.download_button(
                    label=f"Download cluster + member metadata ({st.session_state['cluster_metadata_zip_info']})",
                    data=st.session_state["cluster_metadata_zip"],
                    file_name="cluster_metadata.zip",
                    mime="application/zip",
                )

        # -- Results table --

        cluster_event = st.dataframe(
            cluster_results_df,
            use_container_width=True,
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row",
            key="cluster_results_table"
        )

        cluster_selected_rows = cluster_event.selection.rows

        if not cluster_selected_rows:
            st.info("Select a cluster from the table above to view its members.")
        else:
            cluster_selected_idx = cluster_selected_rows[0]
            cluster_selected_row = cluster_results_df.iloc[cluster_selected_idx]

            sel_organism = cluster_selected_row["organism"]
            sel_cluster_number = cluster_selected_row["cluster_number"]
            sel_centroid = cluster_selected_row.get("centroid", "")
            sel_n_sequences = cluster_selected_row.get("n_sequences", 0)
            sel_full_size = cluster_selected_row.get("full_cluster_size", sel_n_sequences)

            st.divider()
            st.subheader(f"Cluster: {sel_cluster_number}")

            info_cols = st.columns(4)
            info_cols[0].metric("Organism", sel_organism)
            info_cols[1].metric("Cluster number", sel_cluster_number)
            info_cols[2].metric("Members (filtered)", sel_n_sequences)
            info_cols[3].metric("Members (total)", sel_full_size)

            sel_n_assemblies = cluster_selected_row.get("n_assemblies")
            sel_dominant_assembly = cluster_selected_row.get("dominant_assembly")
            sel_assemblies_preview = cluster_selected_row.get("assemblies_preview")

            if pd.notna(sel_n_assemblies):
                assembly_cols = st.columns(2)
                assembly_cols[0].metric("Assemblies", int(sel_n_assemblies))
                assembly_cols[1].metric(
                    "Dominant assembly",
                    sel_dominant_assembly if pd.notna(sel_dominant_assembly) else "-",
                )

            st.write(f"**Centroid:** {sel_centroid}")

            if isinstance(sel_assemblies_preview, str) and sel_assemblies_preview:
                st.caption(f"Assemblies: {sel_assemblies_preview}")

            members_seq_df = load_cluster_member_sequences(
                db_path=db_path,
                organism=sel_organism,
                cluster_number=sel_cluster_number,
                segment_filter=segment_filter,
                genotype_filter=genotype_filter,
                ha_subtype_filter=ha_subtype_filter,
                na_subtype_filter=na_subtype_filter,
                host_filter=host_filter,
                country_filter=country_filter,
            )

            st.write(f"**Sequences loaded:** {len(members_seq_df)}")

            if not members_seq_df.empty:
                multifasta_text = format_multifasta(members_seq_df)

                st.download_button(
                    label=f"Download multifasta ({len(members_seq_df)} sequences)",
                    data=multifasta_text,
                    file_name=f"cluster_{sel_cluster_number}.fasta",
                    mime="text/plain",
                )

                with st.expander("Preview sequences"):
                    preview_cols = [
                        column
                        for column in [
                            "accession",
                            "description",
                            "assembly",
                            "identity_to_centroid",
                        ]
                        if column in members_seq_df.columns
                    ]

                    st.dataframe(
                        members_seq_df[preview_cols].head(50),
                        use_container_width=True,
                        hide_index=True
                    )

            with st.expander("Contributing assemblies"):
                assemblies_df = load_cluster_assembly_values(
                    db_path=db_path,
                    organism=sel_organism,
                    cluster_number=sel_cluster_number,
                )

                if assemblies_df.empty:
                    st.info("No assembly data available for this cluster.")
                else:
                    if pd.notna(sel_n_assemblies) and len(assemblies_df) < int(
                        sel_n_assemblies
                    ):
                        st.caption(
                            f"Showing the {len(assemblies_df)} largest of "
                            f"{int(sel_n_assemblies)} assemblies."
                        )

                    st.dataframe(
                        assemblies_df,
                        use_container_width=True,
                        hide_index=True
                    )

            with st.expander("Linked clusters"):
                linked_df = get_assembly_links(
                    db_path=db_path,
                    organism=sel_organism,
                    cluster_number=sel_cluster_number,
                )

                if linked_df.empty:
                    st.info(
                        "No clusters are linked to this cluster through "
                        "shared assemblies."
                    )
                else:
                    st.caption(
                        "Clusters containing sequences from assemblies also "
                        "found in this cluster, strongest link first."
                    )
                    st.dataframe(
                        linked_df,
                        use_container_width=True,
                        hide_index=True
                    )
