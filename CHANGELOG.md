# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [2026-10-02]

### Fixed

- **README import recipes referenced files that do not exist.** The IAV recipe
  pointed at `input_data/iav_metadata_0426.csv` and `input_data/iav_clusters.clstr`;
  the actual inputs are `input_data/iav_metadata_0526.csv` and
  `input_data/iav_sequences.clstr`. Corrected in both the recipe and the
  "Input data format" listing.
- **Both README recipes used `--recreate` against the same database file, so the
  second one destroyed the first.** `initialise_database` unlinks the target
  before loading, and every loader is a plain `INSERT INTO`, so the correct
  sequence is: first species with `--recreate`, subsequent species appended
  without it. The IAV/IBV recipes are fixed and the constraint is now stated
  explicitly, including that one ETL run handles a single organism.
- **A re-import silently dropped NCBI-fetched sequences.** Sequences downloaded
  by `db/ncbi_fetch.py` are written to `data/fetched_sequences.fasta`, which is
  not the `--fasta` input, so rebuilding from `input_data/iav_seq.fasta` alone
  loses them and re-reports them as missing. Added a "Reusing previously fetched
  sequences" section covering both the concatenation and `--fetch-missing`
  routes, cross-linked from the fetch documentation.
- **The "database must be built with the current schema" warning named the wrong
  thing.** It told users to re-import when `sequences` lacks `n_N` /
  `n_degenerate`, but those columns were already present; what was actually
  missing were the derived `cluster_assembly` and `cluster_assembly_links`
  tables. Replaced with a one-line `sqlite3` check for those tables.

### Changed

- `seq_get` now streams the FASTA line by line and normalises each record as it
  completes, instead of calling `file.readlines()` and then walking the resulting
  list. The previous form held every line of the file in memory as a separate
  string object simultaneously with the growing result dictionary.
- Sequence normalisation is performed by one compiled regex
  (`_NON_NUCLEOTIDE`) replacing a Python-level loop that iterated every
  character of every sequence and, for each non-IUPAC residue, rebuilt the whole
  string via `str.replace`.
- Added `build_sequence_table`, which fuses `decompose_fasta_headers` +
  `accession_to_seq_dict` + `join_fasta_dicts` into a single pass so the two
  intermediate dictionaries are no longer materialised. Duplicate accessions
  still resolve to the last occurrence, matching the previous chain exactly.
- `transform_seqdict` reuses the `per_base` counts it has already computed
  instead of calling `degenerate_breakdown_json`, which recomputed them. Added
  `utils.sequence_analysis.ambiguous_breakdown_json` so the JSON serialisation
  convention still lives in one place; `degenerate_breakdown_json` now delegates
  to it. `db.ncbi_fetch` is unchanged.
- `main` releases the FASTA-derived table before the cluster tables are built so
  the two are not resident at the same time.
- README now documents the ETL's memory profile: peak RSS is driven by the FASTA
  and by the number of rows held at once, not by the output database size, and a
  full IAV import needs roughly 6 GB of available RAM.

### Performance

Measured on a full IAV import (994,840 sequences, 1.78 GB FASTA, 12 cores), same
inputs, original versus optimised code:

| Metric | Before | After | Change |
|---|---|---|---|
| Peak RSS | 6.20 GB | 5.02 GB | -19% |
| Wall clock | 6:10 | 4:43 | -23% |
| User CPU | 285 s | 204 s | -29% |

Isolated `seq_get` benchmark on 270 MB / 150,000 records: 0.68 GB -> 0.40 GB
peak RSS and 7.3 s -> 1.6 s.

The remaining peak is no longer the FASTA parser. It is the load phase, where
every transformed row list is resident simultaneously — `validate_cross_references`
requires metadata, sequences, clusters and cluster composition at once. Reducing
it further needs the ETL restructured to stream rows into SQLite, which was
deliberately left out of scope here.

### Verification

- Differential harness comparing the new `seq_get` against the previous
  implementation imported from git: byte-identical on 150,000 real records /
  249M bases, and on a synthetic file covering semicolon headers, duplicate
  headers, embedded `>`, empty and blank-only records, CRLF line endings, every
  IUPAC code, non-ASCII residues (`ß` expanding to `SS`, `Ω` mapped to `N`) and
  internal whitespace.
- Full IAV re-import with the optimised code produced a database whose schema
  and all ten tables are content-identical to one built by the previous code,
  compared by streaming SHA-256 over every row of 994,840 rows per large table.
- `PRAGMA integrity_check` returns `ok`; the ETL's own record- and
  assembly-completeness checks report zero gaps for both species.

### Notes

- The local `data/sequence_database.sqlite` was rebuilt during this work because
  it predated `cluster_assembly` and `cluster_assembly_links`. Database files are
  gitignored, so this does not appear as a change in the repository.