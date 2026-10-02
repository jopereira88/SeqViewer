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

> The database must be built with the current schema. A quick way to check is to
> confirm `cluster_assembly` and `cluster_assembly_links` exist — they are
> derived at load time, so a database predating them cannot be used with the
> cluster association features and must be re-imported:
>
> ```bash
> sqlite3 data/sequence_database.sqlite \
>     "select name from sqlite_master where type='table' and name like 'cluster_assembly%';"
> ```

## Update the database

The ETL pipeline reads raw metadata (CSV), sequences (FASTA), and CD-HIT clusters (.clstr), transforms them into a normalised schema, and loads everything into SQLite.

> **One species per run, one database.** Each run handles a single organism —
> `infer_single_organism_from_metadata` raises if the metadata contains more than
> one. Run the first species with `--recreate` and every subsequent species
> **without** it: loaders use plain `INSERT INTO`, so the run appends to the
> existing database. Adding `--recreate` to a second species deletes the database
> (`initialise_database` unlinks the file) and discards the species loaded first.

> **Memory.** The pipeline transforms each source fully in memory before loading,
> so peak RSS is driven by the FASTA and by the number of rows held at once, not
> by the size of the output database. An IAV import (~995K sequences, 1.8 GB
> FASTA) peaks at roughly 5 GB and takes about 5 minutes on 12 cores; budget
> ~6 GB of available RAM or the run will fall back to swap and slow to a crawl.
> The IBV import is far lighter (~83K sequences) and needs well under 1 GB.
> Building into a throwaway path first (`--db data/tmp.sqlite`) is the safe way
> to test a rebuild, since `--recreate` unlinks the target before loading.

### Influenza A virus (IAV)

`--recreate` creates the database. Run this one first.

```bash
.venv/bin/python db/import_data.py \
    --db data/sequence_database.sqlite \
    --schema db/schema.sql \
    --metadata input_data/iav_metadata_0526.csv \
    --fasta input_data/iav_seq.fasta \
    --clstr input_data/iav_sequences.clstr \
    --load-mode transactional \
    --recreate
```

### Influenza B virus (IBV)

No `--recreate`: this appends IBV to the IAV database created above.

```bash
.venv/bin/python db/import_data.py \
    --db data/sequence_database.sqlite \
    --schema db/schema.sql \
    --metadata input_data/ibv_sequences.csv \
    --fasta input_data/ibv_sequences.fasta \
    --clstr input_data/ibv_sequences.clstr \
    --load-mode transactional
```

### Reusing previously fetched sequences

Sequences downloaded from NCBI by `db/ncbi_fetch.py` are written to
`data/fetched_sequences.fasta`, **not** to the original `--fasta` input. A
re-import that points at `input_data/iav_seq.fasta` alone silently loses them and
re-reports them as missing. Either concatenate them into the FASTA passed to
`--fasta`:

```bash
cat input_data/iav_seq.fasta data/fetched_sequences.fasta > data/iav_full.fasta
```

or re-run the IAV import with `--fetch-missing --entrez-email you@example.org`,
which resumes from `data/fetch_progress.json` and re-uses the already-downloaded
sequences without re-fetching them.

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

When `--db` is used, each successfully fetched sequence is also inserted into
the `sequences` table of that database (upsert) as it is downloaded, so the
database reflects the fetch immediately. Sequences fetched by an earlier run
(present in the fetched FASTA but missing from the database) are backfilled
before fetching starts, so re-running the same command closes the gap without
re-downloading anything. Accessions that cannot be inserted (e.g. missing from
the `metadata` table, or longer than the `sequences` table CHECK allows) are
skipped and reported. Pass `--no-update-db` to keep the original FASTA-only
behaviour.

The fetched FASTA can be concatenated with the original file and passed to the
ETL as `--fasta`, or the ETL can be run with `--fetch-missing` (which reuses the
same resumable download). See
[Reusing previously fetched sequences](#reusing-previously-fetched-sequences).

### Record completeness check

The ETL checks that every sequence is present in all three sources *of the same
species*: metadata, cluster composition and the FASTA file. Gaps are reported
per species. Example output when accessions are missing from the FASTA file:

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
├── iav_metadata_0526.csv    # IAV metadata (semicolon-delimited)
├── iav_seq.fasta            # IAV sequences
├── iav_sequences.clstr      # IAV CD-HIT clusters
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
| `cluster_assembly` | Assemblies contributing to each cluster, with per-assembly member counts |
| `cluster_assembly_links` | Precomputed cluster-to-cluster associations through shared assemblies |
| `cluster_member_metadata` | Convenience view joining all member data |

### Assembly associations

`assembly` records the genome assembly each sequence came from. Two clusters
are considered associated when they contain sequences from the same assembly,
which is the signal that they may represent the same biological entity across
separate submissions.

Two derived tables back this:

- `cluster_assembly` — one row per (cluster, segment, assembly) with
  `n_sequences` and `is_centroid_assembly`.
- `cluster_assembly_links` — one row per *unordered* cluster-cluster pair,
  stored once with `n_connections` (the number of shared assemblies). The
  smaller cluster key is stored first, so a lookup must consider both the
  `cluster_number` and the `linked_cluster_number` side of every row.

The ETL checks the invariant that every assembly contributes one cluster per
segment; violations are reported by the record completeness check and are fatal
under `--record-completeness fail`.

**Why `assembly` is not in `cluster_filter_counts`.** That table's composite
primary key would grow from ~114k to ~27M rows, because assembly cardinality is
sequence-level. The assembly filter is therefore applied as an `EXISTS`
semi-join against `cluster_assembly`, which composes with both the filtered and
unfiltered cluster-search branches and leaves the filter chain untouched.

In the Clusters tab, the **Search assembly** box restricts results to clusters
containing sequences from assemblies matching the accession (partial matches
allowed). Each result row carries `n_assemblies`, `dominant_assembly` and
`assemblies_preview`. Selecting a cluster shows its contributing assemblies and
the clusters linked to it, strongest link first. Link lookups are served from an
in-memory index built once per database via `@st.cache_resource`.

The **Generate metadata TSV zip** download covers all matching clusters with the
active filters applied and produces two files: `clusters.tsv` (one row per
cluster, including the assembly summary) and `members.tsv` (one row per
sequence, including its assembly). FASTA headers from the cluster exports also
carry `assembly=<accession>`.

### Bulk download performance

Exports that cover the whole dataset are large: over a million sequences and
1.8 G bases. Two things keep them tractable.

Cluster pair restrictions (the explicit list of clusters an export covers) are
staged in a temp table and joined, rather than expanded into an `OR` chain over
`(organism, cluster_number)`. The chain form reads as an index-friendly
equality test but gives the planner a disjunction it can only satisfy with a
full scan, so it cost minutes per export.

Sequence and member rows are streamed straight from the cursor into the archive
entry, instead of being collected into a dataframe first. `clusters.tsv` is
still small enough (one row per cluster) to serialise as a dataframe. Streaming
is what makes an unfiltered export possible: the dataframe form of the sequence
export needs several gigabytes of RAM, against roughly 0.3 GB streamed.
Compression is DEFLATE at the default level, which is slower than a low level
but a good deal smaller (95 MB versus 137 MB for a full export).

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
│   ├── assembly_service.py      # Assembly associations + cached link index
│   ├── cluster_service.py      # Cluster business logic
│   ├── metadata_service.py     # Metadata/filter queries
│   └── sequence_service.py     # Sequence queries
├── utils/
│   ├── __init__.py
│   ├── fasta.py                # FASTA formatting
│   ├── sequence_analysis.py    # Ambiguous base counting, accession parsing
│   └── tables.py               # TSV/CSV export and zipped table downloads
├── data/
│   └── sequence_database.sqlite
└── input_data/
    └── *.csv, *.fasta, *.clstr
```
