"""
Metadata service layer.

Business logic for filter management and metadata queries.
"""

import db.queries as _q


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
    return _q.get_distinct_values(db_path=db_path, table=table, column=column)


def get_filter_options(db_path, table="metadata"):
    """
    Return a dictionary of filter option lists for the sidebar.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    table : str
        Table to query distinct values from.

    Returns
    -------
    dict
        Mapping of column name to list of distinct values.
    """
    return {
        "segment": get_distinct_values(db_path, table, "segment"),
        "genotype": get_distinct_values(db_path, table, "genotype"),
        "ha_subtype": get_distinct_values(db_path, table, "ha_subtype"),
        "na_subtype": get_distinct_values(db_path, table, "na_subtype"),
        "host": get_distinct_values(db_path, table, "host"),
        "country": get_distinct_values(db_path, table, "country"),
    }
