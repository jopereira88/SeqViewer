"""
Assembly service layer.

Business logic for assembly-based cluster associations.

Clusters are associated through shared genome assemblies. Each unordered
cluster-cluster link is stored once in ``cluster_assembly_links``, so a
lookup for a given cluster must consider both the ``cluster_number`` side
and the ``linked_cluster_number`` side of every row. The in-memory index
below precomputes both directions once per database, which keeps repeated
detail-panel lookups cheap.
"""

import sys

import pandas as pd
import streamlit as st

import db.queries as _q


LINK_OUTPUT_COLUMNS = [
    "linked_cluster",
    "linked_cluster_segment",
    "n_connections",
]


class AssemblyLinkIndex:
    """
    Read-only, bidirectional cluster-to-cluster association index.

    Parameters
    ----------
    links_by_cluster : dict[tuple[str, str], list[tuple]]
        Mapping from ``(organism, cluster_number)`` to the clusters linked
        to it, each as ``(linked_organism, linked_cluster,
        linked_cluster_segment, n_connections)`` ordered by
        ``n_connections`` descending.

    n_links : int
        Number of distinct unordered cluster-cluster links indexed.
    """

    def __init__(self, links_by_cluster, n_links):
        self._links_by_cluster = links_by_cluster
        self.n_links = n_links

    def __len__(self):
        """
        Return the number of indexed clusters.
        """
        return len(self._links_by_cluster)

    def has_links(self, organism, cluster_number):
        """
        Check whether a cluster has any assembly-based associations.

        Parameters
        ----------
        organism : str
            Cluster organism.

        cluster_number : str
            Cluster number.

        Returns
        -------
        bool
            True if the cluster has at least one linked cluster.
        """
        return bool(
            self._links_by_cluster.get((organism, str(cluster_number)))
        )

    def links_for(self, organism, cluster_number):
        """
        Return the clusters linked to one cluster, strongest first.

        Parameters
        ----------
        organism : str
            Cluster organism.

        cluster_number : str
            Cluster number.

        Returns
        -------
        pandas.DataFrame
            Columns: linked_cluster, linked_cluster_segment, n_connections.
        """
        rows = self._links_by_cluster.get((organism, str(cluster_number)))

        if not rows:
            return pd.DataFrame(columns=LINK_OUTPUT_COLUMNS)

        return pd.DataFrame.from_records(
            rows,
            columns=[
                "linked_organism",
                "linked_cluster",
                "linked_cluster_segment",
                "n_connections",
            ],
        )[LINK_OUTPUT_COLUMNS]


@st.cache_resource(show_spinner=False)
def load_assembly_link_index(db_path):
    """
    Build and cache the in-memory cluster-to-cluster association index.

    The index is built once per database and shared across sessions. It
    degrades to an empty index when the association tables are absent.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    Returns
    -------
    AssemblyLinkIndex
        Bidirectional association index.
    """
    df = _q.load_cluster_assembly_table(db_path)

    if df.empty:
        return AssemblyLinkIndex({}, 0)

    # The association table stores only the linked cluster's segment, so
    # the present cluster's segment is resolved from cluster_summary.
    # That is what makes the reverse direction correct.
    segment_df = _q.load_cluster_segment_lookup(db_path)
    segment_by_cluster = {
        (sys.intern(str(organism)), sys.intern(str(cluster_number))):
            sys.intern(str(segment))
        for organism, cluster_number, segment
        in segment_df.itertuples(index=False, name=None)
    }

    links_by_cluster = {}
    n_links = 0

    for org, cluster, link_org, link_cluster, link_segment, n_connections in (
        df.itertuples(index=False, name=None)
    ):
        org_key = sys.intern(str(org))
        cluster_key = sys.intern(str(cluster))
        link_org_key = sys.intern(str(link_org))
        link_cluster_key = sys.intern(str(link_cluster))

        # Reverse direction resolves its partner segment from
        # cluster_summary rather than reusing the stored column, which
        # belongs to the linked cluster.
        present_segment = segment_by_cluster.get((org_key, cluster_key))
        linked_segment = segment_by_cluster.get(
            (link_org_key, link_cluster_key),
        )

        if linked_segment is None:
            linked_segment = (
                None if link_segment is None else sys.intern(str(link_segment))
            )

        n_links += 1
        n_connections = int(n_connections)

        links_by_cluster.setdefault((org_key, cluster_key), []).append(
            (
                link_org_key,
                link_cluster_key,
                linked_segment,
                n_connections,
            )
        )

        links_by_cluster.setdefault((link_org_key, link_cluster_key), []).append(
            (
                org_key,
                cluster_key,
                present_segment,
                n_connections,
            )
        )

    for rows in links_by_cluster.values():
        rows.sort(key=lambda row: row[3], reverse=True)

    return AssemblyLinkIndex(links_by_cluster, n_links)


def get_assembly_links(db_path, organism, cluster_number):
    """
    Load the clusters linked to one cluster through shared assemblies.

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
        Columns: linked_cluster, linked_cluster_segment, n_connections,
        ordered by n_connections descending.
    """
    return load_assembly_link_index(db_path).links_for(organism, cluster_number)


def load_cluster_assembly_values(db_path, organism, cluster_number, limit=200):
    """
    Load the assemblies contributed to one cluster.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    organism : str
        Cluster organism.

    cluster_number : str
        Cluster number.

    limit : int or None
        Maximum number of assemblies to return, strongest first.

    Returns
    -------
    pandas.DataFrame
        Columns: assembly, segment, n_sequences, is_centroid_assembly.
    """
    return _q.load_cluster_assembly_values(
        db_path=db_path,
        organism=organism,
        cluster_number=cluster_number,
        limit=limit,
    )