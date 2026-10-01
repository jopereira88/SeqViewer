"""
Tabular export utilities.

Helpers for turning dataframes into downloadable TSV/CSV payloads and for
bundling the cluster-level and member-level metadata exports into a single
zip archive.
"""

import csv
import io
import zipfile


def df_to_delimited_bytes(df, sep="\t", na_rep=""):
    """
    Serialise a dataframe to delimited text bytes.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to serialise.

    sep : str
        Field separator.

    na_rep : str
        Replacement for missing values.

    Returns
    -------
    bytes
        UTF-8 encoded delimited text, without the index column.
    """
    if df is None or df.empty:
        return b""

    return df.to_csv(
        sep=sep,
        index=False,
        na_rep=na_rep,
        lineterminator="\n",
    ).encode("utf-8")


def df_to_tsv_bytes(df):
    """
    Serialise a dataframe to tab-separated bytes.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to serialise.

    Returns
    -------
    bytes
        UTF-8 encoded TSV, without the index column.
    """
    return df_to_delimited_bytes(df, sep="\t")


def df_to_csv_bytes(df):
    """
    Serialise a dataframe to comma-separated bytes.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to serialise.

    Returns
    -------
    bytes
        UTF-8 encoded CSV, without the index column.
    """
    return df_to_delimited_bytes(df, sep=",")


def build_delimited_zip(entries):
    """
    Build an in-memory zip from named delimited tables.

    Parameters
    ----------
    entries : dict[str, pandas.DataFrame]
        Mapping of archive filename to dataframe. Empty dataframes are
        skipped.

    Returns
    -------
    bytes
        Zip file contents.
    """
    buf = io.BytesIO()

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for filename, df in entries.items():
            if df is None or df.empty:
                continue

            payload = df_to_tsv_bytes(df)

            if not payload:
                continue

            zf.writestr(filename, payload)

    return buf.getvalue()


def cluster_metadata_zip_bytes(clusters_df, members_df):
    """
    Build an in-memory zip with cluster-level and member-level metadata.

    The cluster-level table carries one row per cluster, including its
    assembly summary. The member-level table carries one row per sequence,
    including its assembly.

    Parameters
    ----------
    clusters_df : pandas.DataFrame
        Cluster-level metadata, including assembly summary columns.

    members_df : pandas.DataFrame
        Per-sequence metadata, including the assembly column.

    Returns
    -------
    bytes
        Zip file contents.
    """
    return build_delimited_zip({
        "clusters.tsv": clusters_df,
        "members.tsv": members_df,
    })


def cluster_metadata_zip_bytes_from_rows(
    clusters_df,
    member_columns,
    member_rows
):
    """
    Build the cluster + member metadata zip with a streamed members table.

    The member table has one row per sequence, so a whole-dataset export is
    over a million rows. Writing it straight from a cursor keeps the archive
    byte-identical to cluster_metadata_zip_bytes without holding those rows
    as a dataframe.

    Parameters
    ----------
    clusters_df : pandas.DataFrame
        Cluster-level metadata, including assembly summary columns.

    member_columns : list[str]
        Member-level column names, in the order member_rows provides them.

    member_rows : iterable[tuple]
        Per-sequence metadata rows, including the assembly column.

    Returns
    -------
    tuple[bytes, int]
        Zip file contents and the number of member rows written.
    """
    buf = io.BytesIO()

    member_count = 0

    def _counted_rows():
        nonlocal member_count

        for row in member_rows:
            member_count += 1
            yield row

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        clusters_payload = df_to_tsv_bytes(clusters_df)

        if clusters_payload:
            zf.writestr("clusters.tsv", clusters_payload)

        with zf.open("members.tsv", "w") as entry:
            for fragment in rows_to_tsv_stream(
                columns=member_columns,
                rows=_counted_rows(),
            ):
                entry.write(fragment.encode("utf-8"))

    return buf.getvalue(), member_count


def rows_to_tsv_stream(
    columns,
    rows,
    sep="\t",
    na_rep=""
):
    """
    Serialise row tuples to tab-separated text without a dataframe.

    Parameters
    ----------
    columns : list[str]
        Column names for the header row.

    rows : iterable[tuple]
        Row values in column order.

    sep : str
        Field separator.

    na_rep : str
        Replacement for missing values.

    Yields
    ------
    str
        Text fragments to concatenate. Rows are yielded individually so the
        caller can stream them into a file or archive entry.
    """
    buf = io.StringIO()

    writer = csv.writer(buf, delimiter=sep, lineterminator="\n")

    writer.writerow(list(columns))

    yield buf.getvalue()

    for row in rows:
        buf.seek(0)
        buf.truncate(0)

        writer.writerow([na_rep if value is None else value for value in row])

        yield buf.getvalue()


def build_delimited_zip_from_rows(filename, columns, rows, sep="\t", na_rep=""):
    """
    Build an in-memory zip with a single streamed delimited table.

    Parameters
    ----------
    filename : str
        Archive entry name.

    columns : list[str]
        Column names for the header row.

    rows : iterable[tuple]
        Row values in column order.

    sep : str
        Field separator.

    na_rep : str
        Replacement for missing values.

    Returns
    -------
    bytes
        Zip file contents.
    """
    buf = io.BytesIO()

    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        with zf.open(filename, "w") as entry:
            for fragment in rows_to_tsv_stream(
                columns=columns,
                rows=rows,
                sep=sep,
                na_rep=na_rep,
            ):
                entry.write(fragment.encode("utf-8"))

    return buf.getvalue()