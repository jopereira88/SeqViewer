# Sequence Viewer

Streamlit application for browsing influenza virus sequence metadata, FASTA sequences, and CD-HIT cluster membership.

Supports **Influenza A virus (IAV)** (~1M sequences) and **Influenza B virus (IBV)** (~82K sequences).

Features:
- Search the sequence index by accession / organism / species, filtered by metadata and by precomputed degenerate-base and N content.
- **Multi-accession search**: enable the checkbox in the sidebar and paste a semicolon-separated list of accessions. Results show per-sequence `n_N` and `n_degenerate` counts and a per-base breakdown, and can be downloaded as a multifasta.
- Cluster browser with pagination and bulk FASTA/ZIP downloads.

## Requirements

- Python 3.13+

## Setup

```bash
python3.13 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Run the app

```bash
.venv/bin/streamlit run app.py
```

The app reads from `data/sequence_database.sqlite` by default. To use a different database, change the path in the sidebar.

> The database must be built with the current schema: the `sequences` table stores
> the ETL-precomputed ambiguous base counts (`n_N`, `n_degenerate`). Re-import the
> database if it predates these columns.

## Update the database

The ETL pipeline reads raw metadata (CSV), sequences (FASTA), and CD-HIT clusters (.clstr), transforms them into a normalised schema, and loads everything into SQLite.

### Influenza A virus (IAV)

```bash
.venv/bin/python db/import_data.py \
    --db data/sequence_database.sqlite \
    --schema db/schema.sql \
    --metadata input_data/iav_metadata_0426.csv \
    --fasta input_data/iav_seq.fasta \
    --clstr input_data/iav_clusters.clstr \
    --load-mode transactional \
    --recreate
```

### Influenza B virus (IBV)

```bash
.venv/bin/python db/import_data.py \
    --db data/sequence_database.sqlite \
    --schema db/schema.sql \
    --metadata input_data/ibv_sequences.csv \
    --fasta input_data/ibv_sequences.fasta \
    --clstr input_data/ibv_sequences.clstr \
    --load-mode transactional \
    --recreate
```

### CLI options

| Flag | Description |
|---|---|
| `--db` | Path to SQLite database file |
| `--schema` | Path to `schema.sql` |
| `--metadata` | Path to metadata CSV (semicolon-delimited) |
| `--fasta` | Path to FASTA sequences file |
| `--clstr` | Path to CD-HIT `.clstr` cluster file |
| `--load-mode` | `transactional` (default, all-or-nothing) or `basic` (table-by-table commits) |
| `--recreate` | Delete and recreate the database before loading |
| `--metadata-delimiter` | Delimiter for metadata CSV (default: `;`) |
| `--record-completeness` | `fail` aborts the load if any sequence lacks metadata, a cluster composition entry or a FASTA sequence (per species); `warn` prints the report and continues (default: `warn`) |
| `--fetch-missing` | Download sequences missing from the FASTA file from NCBI before loading (requires `--entrez-email`; resumable via `--fetch-output`/`--fetch-progress`) |
| `--entrez-email` | Email address required by NCBI Entrez |
| `--entrez-api-key` | NCBI Entrez API key (optional, raises rate limits) |
| `--fetch-output` | FASTA file fetched sequences are appended to (default: `data/fetched_sequences.fasta`) |
| `--fetch-progress` | JSON progress file used to resume an interrupted download (default: `data/fetch_progress.json`) |
| `--fetch-batch-size` | Accessions per NCBI `efetch` request (default: 200) |
| `--fetch-sleep` | Seconds between NCBI requests (default: 0.34 without an API key, 0.1 with one) |

### Fetching missing sequences from NCBI (resumable)

Sequences that are present in metadata and cluster composition but missing from
the FASTA file can be downloaded from NCBI. Because large downloads take a long
time, this is a standalone, resumable CLI: stop it any time (Ctrl-C) and re-run
the same command to continue where it left off.

```bash
.venv/bin/python db/ncbi_fetch.py \
    --db data/sequence_database.sqlite \
    --output data/fetched_sequences.fasta \
    --progress data/fetch_progress.json \
    --email you@example.org \
    [--api-key XXXXXXXX]
```

Accessions can also be supplied with `--accessions-file` (one per line) or
`--accessions` (comma/semicolon separated). Pass `--fasta <original fasta>` to
detect the dominant header format so fetched headers conform to it. Fetched
headers use the same `>ACCESSION |description` format as the source files.

The fetched FASTA can be concatenated with the original file and passed to the
ETL as `--fasta`, or the ETL can be run with `--fetch-missing` (which reuses the
same resumable download).

### Record completeness check

The ETL checks that every sequence is present in all three sources *of the same
species*: metadata, cluster composition and the FASTA file. Gaps are reported
per species. Example output for the shipped database:

```
Record completeness gaps
------------------------
Accessions missing a FASTA sequence: 4176
  Species 'Influenza A virus':
    missing from FASTA (4176): PZ405244.1, PZ405245.1, ...
```

### Input data format

```
input_data/
├── iav_metadata_0426.csv    # IAV metadata (semicolon-delimited)
├── iav_seq.fasta            # IAV sequences
├── iav_clusters.clstr       # IAV CD-HIT clusters
├── ibv_sequences.csv        # IBV metadata
├── ibv_sequences.fasta      # IBV sequences
└── ibv_sequences.clstr      # IBV CD-HIT clusters
```

### Database tables

| Table | Contents |
|---|---|
| `metadata` | Sequence metadata (accession, organism, segment, genotype, host, country, etc.) |
| `sequences` | FASTA sequences linked by accession, with precomputed ambiguous base counts (`n_N`, `n_degenerate`) and a JSON `degenerate_breakdown` |
| `clusters` | Cluster centroids |
| `cluster_composition` | Member-to-cluster mapping with identity percentages |
| `cluster_summary` | Precomputed cluster statistics |
| `cluster_background_counts` | Host/country/genotype counts per cluster |
| `cluster_descriptors` | JSON value distributions, dominant values, and flags |
| `cluster_filter_counts` | Precomputed filter combination counts |
| `cluster_member_metadata` | Convenience view joining all member data |

## Project structure

```
seqviewer/
├── app.py                      # Streamlit UI
├── requirements.txt
├── db/
│   ├── __init__.py
│   ├── connection.py           # Shared SQLite connection
│   ├── queries.py              # Database query functions
│   ├── schema.sql              # Table and index definitions
│   ├── import_data.py          # ETL pipeline
│   └── ncbi_fetch.py           # Resumable NCBI sequence download CLI
├── services/
│   ├── __init__.py
│   ├── cluster_service.py      # Cluster business logic
│   ├── metadata_service.py     # Metadata/filter queries
│   └── sequence_service.py     # Sequence queries
├── utils/
│   ├── __init__.py
│   ├── fasta.py                # FASTA formatting
│   └── sequence_analysis.py    # Ambiguous base counting, accession parsing
├── data/
│   └── sequence_database.sqlite
└── input_data/
    └── *.csv, *.fasta, *.clstr
```
