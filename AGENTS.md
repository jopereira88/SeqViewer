# AGENTS.md — SeqViewer

## What this is

Streamlit monolith for browsing influenza virus sequence metadata, FASTA
sequences, and CD-HIT clusters (IAV ~1M + IBV ~82K sequences). Python 3.13+,
deps: `streamlit`, `pandas`, `biopython` only. No frontend build, no API
layer, no test suite, no linter config.

## Commands

The virtualenv in this checkout is `.view/` (the README says `.venv`; use
whatever exists):

```bash
.view/bin/streamlit run app.py                       # run the app
.view/bin/python -m py_compile app.py db/queries.py  # syntax check
.view/bin/python db/import_data.py --help            # ETL (see README)
```

- Data lives in `data/sequence_database.sqlite` (~4.3 GB, gitignored).
- ETL: **one species per run, one DB.** First species with `--recreate`,
  every later species without it (it appends; `--recreate` unlinks the DB).
- `requirements.txt` has exactly: streamlit, pandas, biopython.

## Architecture

```
app.py              Entire Streamlit UI: sidebar filters, Sequences tab, Clusters tab
db/schema.sql       Tables, cluster_member_metadata view, indexes
db/import_data.py   ETL: semicolon CSV + FASTA + .clstr -> SQLite, then derived tables
db/queries.py       ALL read queries + ALL filter SQL; @st.cache_data on entry points
services/*.py       Thin pass-through wrappers over db.queries (no logic)
utils/              Pure helpers: fasta.py, tables.py, sequence_analysis.py
input_data/         Raw metadata CSVs, FASTA, CD-HIT .clstr (gitignored)
```

Data flow: sidebar widget → Python local var → `*_filter` kwarg threaded
through `services/*` → `db/queries.py` → SQLite.

## Adding a metadata filter (the core pattern)

Follow the `genbank_refseq_filter` / `assembly_filter` precedent exactly:

1. **UI** (`app.py`, sidebar "Metadata filters" block): add the widget.
   Free-text lists use `utils.sequence_analysis.parse_accession_list`
   (`; , \n \t` + spaces, de-duped, order-preserving).
2. **SQL helper** (`db/queries.py`): add the kwarg to
   `_append_metadata_filters` (equality `= ?`, or `IN (...)` for lists) and/or
   `_build_sequence_index_where`. Give it an inactive default (`"All"` or
   `()`) so **cluster callers that never pass it are unaffected** — this is
   the load-bearing convention; several cluster queries share these helpers.
3. **Thread it** through the sequence-level query functions in
   `db/queries.py`: `count_sequences`, `load_sequence_index`,
   `iter_sequence_index_sequences`, `load_sequence_records_by_accessions`
   (all take `@st.cache_data` kwargs — use hashable defaults).
4. **Service wrappers** (`services/sequence_service.py`): same kwarg in
   signature + forward verbatim.
5. **App call sites** (`app.py`): pass it at every call of the above
   (count, paged index, multi-FASTA export, multi-accession branch).
6. **Cache keys** (`app.py`): add the value to `seq_filter_key` and
   `acc_cache_key` (and `cluster_filter_key` if it affects clusters) or
   pagination/stale downloads will not reset.
7. **CHANGELOG.md**: entry under `[Unreleased] > Added` (Keep a Changelog).

## Cluster filtering gotchas

- `cluster_filter_counts` is a precomputed rollup of filter combinations for
  segment/genotype/ha/na/host/country only. `_metadata_filters_active` and
  `_append_filter_count_filters` branch on it.
- **`assembly` and `genbank_refseq` are deliberately absent** from that table
  (sequence-level cardinality would explode it ~114K → ~27M rows; see
  schema.sql). Assembly is instead applied on cluster queries as an
  `EXISTS` semi-join against `cluster_assembly` via `_append_assembly_filter`
  (substring `LIKE`), which composes with both filtered/unfiltered branches.
  The Clusters tab's in-tab assembly search (`app.py`, key
  `assembly_search_input`) is single-value substring; the sidebar assembly
  filter is exact and Sequences-tab only.
- Adding a *new* dimension to `cluster_filter_counts` requires an ETL
  rebuild (derived tables are built in `db/import_data.py`).

## Caching

`@st.cache_data` on every query entry point; changing any kwarg creates a new
cache entry. Session-state payloads (`seq_fasta_zip`, `acc_data`, …) are
invalidated manually by comparing filter tuples in `app.py`.

## Verification (no test suite exists)

```bash
.view/bin/python -m py_compile app.py db/queries.py services/*.py utils/*.py
```

- Query logic: small read-only script asserting `count_sequences(...)` equals
  a direct `SELECT COUNT(*)` for the same predicate (connect to the real DB,
  it is read-only for checks).
- UI: `from streamlit.testing.v1 import AppTest` —
  `AppTest.from_file("app.py", default_timeout=300)` (absolute path; chdir to
  repo root first), then `at.run()` and assert `not at.exception`. First run
  takes ~15-30 s (cold query cache), later runs are fast.
- Never run the ETL against `data/sequence_database.sqlite` casually
  (`--recreate` unlinks it). Test rebuilds go to `--db data/tmp.sqlite`.

## Conventions

- Numpy-style docstrings on public functions; inline comments explain *why*.
- Do not commit `data/`, `input_data/`, `.view/`, or `SESSION_CONTEXT.md`
  (gitignored).
- Update README.md when user-facing behaviour changes; CHANGELOG.md for every
  notable change.
