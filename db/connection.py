"""
SQLite connection management.

Provides a single, cached connection factory used by all query modules.
"""

import sqlite3

import streamlit as st


@st.cache_resource
def connect_sqlite(db_path):
    """
    Open and cache a SQLite connection.

    Parameters
    ----------
    db_path : str
        Path to SQLite database.

    Returns
    -------
    sqlite3.Connection
        SQLite connection with Row factory.
    """
    conn = sqlite3.connect(db_path, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn
