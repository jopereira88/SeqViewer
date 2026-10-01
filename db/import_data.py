#!/usr/bin/env python3
"""
Import sequence, metadata and cluster data into a SQLite database.

Expected project structure
--------------------------
sequence_explorer/
│
├── data/
│   └── sequence_database.sqlite
│
├── db/
│   ├── schema.sql
│   └── import_data.py
│
└── input_data/
    ├── metadata.csv
    ├── sequences.fasta
    └── clusters.clstr

Usage
-----
python db/import_data.py \
    --db data/sequence_database.sqlite \
    --schema db/schema.sql \
    --metadata input_data/metadata.csv \
    --fasta input_data/sequences.fasta \
    --clstr input_data/clusters.clstr \
    --load-mode transactional \
    --recreate
"""

import argparse
import csv
import json
import sqlite3
from pathlib import Path
import re
import statistics
from collections import Counter, defaultdict

import os
import sys

# Allow running this module as a script (python db/import_data.py) while
# importing the top-level 'utils' and 'db' packages.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.sequence_analysis import count_ambiguous_bases, degenerate_breakdown_json
from db import ncbi_fetch

########## PARSER AND TRANSFORMER FUNCS ##################
def seq_get(filename):
    '''Parses fasta file and returns accession numbers and sequences
    Accepts: filename -  a string for the path of the fasta file
    Returns: a dictonary with accession number and sequence
    Raises:FileNotFoundError if file path is wrong'''
    seqs={}
    iupac_nucleotides = {'A', 'C', 'G', 'T', 'U', 'R', \
                         'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'}   
    with open(filename,'r') as file:
        fasta=file.readlines()
    for i in range(len(fasta)):
        if '>' in fasta[i]:
            name=fasta[i].strip()
            name=name.replace(';','_')
            seqs[name]=''
        else:
            seqs[name]+=fasta[i].strip().upper().replace('-','').replace('X','N')
    for key in seqs:
        seqs[key]=seqs[key].replace('U','T')
        for nuc in seqs[key]:
            if nuc not in iupac_nucleotides:
                seqs[key]=seqs[key].replace(nuc,'N')
    return seqs

def parse_clstr(filename, access_only=True):
    '''parses a .clstr file and retrieves the accessions within each cluster
    Accepts: filename(str) path to the clstr file
    access_only (bool) True for accession only, false for the entire row
    Returns: clusters(dict)'''
    clusters = {}
    current_cluster = None

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = line.strip()
                clusters[current_cluster] = []
            else:
                if not access_only:
                    clusters[current_cluster].append(line.strip())
                else:
                    line=line.split(' ')
                    line=line[1]
                    line=line.replace('...','')
                    line=line.replace('>','')
                    line=line.split('_|')[0]
                    clusters[current_cluster].append(line)
    return clusters

def split_fasta_header(header, header_desc_sep=' |'):
    """
    Split a FASTA header into accession and description.

    Expected preferred format
    -------------------------
    >ACCESSION | description

    Fallback format
    ---------------
    >ACCESSION description

    Returns
    -------
    tuple
        accession, description
    """
    header = str(header).strip()

    if header.startswith(">"):
        header = header[1:]

    if header_desc_sep in header:
        accession, description = header.split(header_desc_sep, 1)
        return accession.strip(), description.strip()

    parts = header.split(maxsplit=1)

    accession = parts[0].strip()

    if len(parts) == 1:
        description = None
    else:
        description = parts[1].strip()

    return accession, description


def decompose_fasta_headers(fasta_dict, header_desc_sep=' |'):
    """
    Return a dictionary with accession: description from FASTA headers.
    """
    acc_desc = {}

    for header in fasta_dict.keys():
        accession, description = split_fasta_header(
            header,
            header_desc_sep=header_desc_sep
        )

        acc_desc[accession] = description

    return acc_desc


def accession_to_seq_dict(fasta_dict, header_desc_sep=' |'):
    """
    Return a dictionary with accession: sequence from FASTA records.
    """
    acc_seq = {}

    for header, seq in fasta_dict.items():
        accession, description = split_fasta_header(
            header,
            header_desc_sep=header_desc_sep
        )

        acc_seq[accession] = seq

    return acc_seq

def join_fasta_dicts(acc_desc,acc_seq):
    table={}
    for acc in acc_desc:
        table[acc]=(acc_desc[acc],acc_seq[acc])
    return table

def parse_clstr_member_line(line, header_desc_sep=' |'):
    """
    Parse one member line from a CD-HIT .clstr file.

    Expected examples
    -----------------
    0	982nt, >NC_026431.1... *
    1	980nt, >ABC123.1... at +/99.50%

    Also tolerates identifiers like:
    0	982nt, >NC_026431.1 | description... *

    Returns
    -------
    tuple
        accession, identity_to_centroid
    """
    identifier_match = re.search(r">(.+?)\.\.\.", line)

    if identifier_match is None:
        raise ValueError(f"Could not extract identifier from .clstr line: {line!r}")

    raw_identifier = identifier_match.group(1).strip()

    accession, _description = split_fasta_header(
        raw_identifier,
        header_desc_sep=header_desc_sep
    )

    if "*" in line:
        identity = 100.0
    else:
        identity_match = re.search(r"([0-9]+(?:\.[0-9]+)?)%", line)

        if identity_match is None:
            raise ValueError(f"Could not extract identity from .clstr line: {line!r}")

        identity = float(identity_match.group(1))

    return accession, identity

def mine_clstr_table(clstr_file, organism, header_desc_sep=' |'):
    """
    Parses a .clstr file retrieving cluster number and cluster representative.

    Returns
    -------
    dict
        {
            organism: {
                cluster_number: centroid_accession
            }
        }
    """
    clusters = parse_clstr(clstr_file, False)
    clus_centroid = {}

    for cl, members in clusters.items():
        centroids = []

        for member_line in members:
            accession, identity = parse_clstr_member_line(member_line,
                                                          header_desc_sep)

            if "*" in member_line:
                centroids.append(accession)

        if len(centroids) == 0:
            raise ValueError(f"No centroid found for cluster {cl}")

        if len(centroids) > 1:
            raise ValueError(f"Multiple centroids found for cluster {cl}: {centroids}")

        clus_centroid[cl] = centroids[0]

    table = {
        organism: {
            cluster_number: centroid
            for cluster_number, centroid in clus_centroid.items()
        }
    }

    return table


def mine_clstr_elements(clstr_file, organism, header_desc_sep= ' |'):
    """
    Parses a .clstr file retrieving cluster elements.

    Returns
    -------
    dict
        {
            organism: {
                accession: (cluster_number, identity_to_centroid)
            }
        }
    """
    clusters = parse_clstr(clstr_file, False)

    table = {}

    for cluster_number, members in clusters.items():
        for member_line in members:
            accession, identity = parse_clstr_member_line(member_line,
                                                          header_desc_sep)

            if accession in table:
                raise ValueError(
                    f"Accession appears in more than one cluster: {accession}"
                )

            table[accession] = (cluster_number, identity)

    return {organism: table}


METADATA_COLUMN_MAP = {
    "ACCESSION": "accession",
    "ORGANISM_NAME": "organism_name",
    "GENBANK_REFSEQ": "genbank_refseq",
    "ASSEMBLY": "assembly",
    "RELEASE_DATE": "release_date",
    "ISOLATE": "isolate",
    "SPECIES": "species",
    "LENGTH": "length",
    "NUC_COMPLETENESS": "nuc_completeness",
    "GENOTYPE": "genotype",
    "SEGMENT": "segment",
    "COUNTRY": "country",
    "HOST": "host",
    "COLLECTION_DATE": "collection_date",
    "COMMENT": "comment",
}


REQUIRED_METADATA_COLUMNS = [
    "ACCESSION",
    "ORGANISM_NAME",
    "GENBANK_REFSEQ",
    "ASSEMBLY",
    "RELEASE_DATE",
    "ISOLATE",
    "SPECIES",
    "LENGTH",
    "NUC_COMPLETENESS",
    "GENOTYPE",
    "SEGMENT",
    "COUNTRY",
    "HOST",
    "COLLECTION_DATE",
    "COMMENT",
]


def clean_empty(value):
    """
    Convert empty strings to None.

    SQLite stores Python None as SQL NULL.
    """
    if value is None:
        return None

    value = str(value).strip()

    if value == "":
        return None

    return value


def parse_integer(value, field_name="integer field"):
    """
    Convert a non-empty value to integer.
    """
    value = clean_empty(value)

    if value is None:
        return None

    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(
            f"Could not convert {field_name} to integer: {value!r}"
        ) from exc


def int_to_iupac_IAV(value):
    '''
    Converts a segment integer value into the IUPAC acronym.

    Args:
        value (str or int): segment value between 1 and 8.

    Returns:
        str: IUPAC acronym for the segment (Influenza A Virus only).
    '''
    iupac = {
        1: 'PB2',
        2: 'PB1',
        3: 'PA',
        4: 'HA',
        5: 'NP',
        6: 'NA',
        7: 'MP',
        8: 'NS'
    }

    try:
        return iupac[int(value)]
    except KeyError:
        print(f"Invalid segment value: {value}, expected 1 to 8")
        return None
    except ValueError:
        print(f"Invalid segment value: {value}, expected integer-like value")
        return None


def int_to_iupac_IBV(value):
    '''
    Converts a segment integer value into the IUPAC acronym.

    Args:
        value (str or int): segment value between 1 and 8.

    Returns:
        str: IUPAC acronym for the segment (Influenza B Virus only).
    '''
    iupac = {
        2: 'PB2',
        1: 'PB1',
        3: 'PA',
        4: 'HA',
        5: 'NP',
        6: 'NA',
        7: 'MP',
        8: 'NS'
    }

    try:
        return iupac[int(value)]
    except KeyError:
        print(f"Invalid segment value: {value}, expected 1 to 8")
        return None
    except ValueError:
        print(f"Invalid segment value: {value}, expected integer-like value")
        return None
####NEW AGENTS MAPPING FUNCTIONS HERE############

################################################3

#### EDIT THIS FUNCTION ON UPDATE #####
def get_segment_converter(species=None, organism_name=None):
    '''
    Selects the appropriate segment conversion function based on species or
    organism_name.

    Args:
        species (str or None): species field from metadata.
        organism_name (str or None): organism_name field from metadata.

    Returns:
        function or None:
            Segment conversion function if organism is recognised.
            None if no converter exists for the organism.
    '''

    species = "" if species is None else str(species).strip().lower()
    organism_name = "" if organism_name is None else str(organism_name).strip().lower()

    organism_text = f"{species} {organism_name}"

    if (
        "influenza a virus" in organism_text
        or "alphainfluenzavirus influenzae" in organism_text
    ):
        return int_to_iupac_IAV

    if (
        "influenza b virus" in organism_text
        or "betainfluenzavirus influenzae" in organism_text
    ):
        return int_to_iupac_IBV

    ##### NEW RULES HERE ####

    return None


def parse_segment(value, species=None, organism_name=None):
    '''
    Converts a raw segment value into the appropriate IUPAC segment acronym,
    using species or organism_name to choose the correct conversion logic.

    Args:
        value (str or int):
            Raw segment value.

        species (str or None):
            Species field from metadata.

        organism_name (str or None):
            Organism name field from metadata. Can include strain information,
            for example:
            "Influenza A virus (A/California/07/2009(H1N1))"

    Returns:
        str or None:
            Converted segment acronym, or None if conversion fails.
    '''

    if value is None:
        return None

    value = str(value).strip()

    if value == "":
        return None

    converter = get_segment_converter(
        species=species,
        organism_name=organism_name
    )

    if converter is None:
        print(
            "No segment converter found for "
            f"species={species!r}, organism_name={organism_name!r}"
        )
        return value

    return converter(value)


def parse_subtypes(genotype, species=None):
    """
    Extract HA and NA subtypes from genotype for Influenza A sequences.

    For Influenza A (Alphainfluenzavirus), genotype contains HxNy patterns.
    For all other organisms, subtypes are not applicable.

    Parameters
    ----------
    genotype : str or None
        Raw genotype value (e.g. 'H5N1', 'H3', 'N2', 'Victoria').
    species : str or None
        Species field used to determine if this is Influenza A.

    Returns
    -------
    tuple
        (ha_subtype, na_subtype) — each is str or None.
    """
    genotype = clean_empty(genotype)

    if genotype is None:
        return None, None

    # Only parse subtypes for Influenza A virus
    if species is None or "Alphainfluenzavirus" not in str(species):
        return None, None

    # Full pair: H5N1
    match = re.match(r'^(H\d+)(N\d+)$', genotype, re.IGNORECASE)
    if match:
        return match.group(1).upper(), match.group(2).upper()

    # HA-only: H5
    match = re.match(r'^(H\d+)$', genotype, re.IGNORECASE)
    if match:
        return match.group(1).upper(), None

    # NA-only: N1
    match = re.match(r'^(N\d+)$', genotype, re.IGNORECASE)
    if match:
        return None, match.group(1).upper()

    return None, None


def validate_metadata_header(fieldnames):
    """
    Validate that the metadata CSV contains all expected columns.
    """
    if fieldnames is None:
        raise ValueError("Metadata file has no header.")

    observed = set(fieldnames)
    expected = set(REQUIRED_METADATA_COLUMNS)

    missing = sorted(expected - observed)

    if missing:
        raise ValueError(f"Metadata file is missing columns: {missing}")

def extract_organism_key(organism_name):
    """
    Extract an operational organism key from ORGANISM_NAME.

    Generic rule
    ------------
    If the word 'virus' appears, return everything from the start of the string
    up to and including the first occurrence of 'virus'.

    The returned key normalises the final word to 'virus' so that:
    - 'Influenza A Virus'
    - 'Influenza A virus'

    collapse to the same key.

    This is used only downstream for cluster grouping.
    """
    if organism_name is None:
        return None

    organism_name = str(organism_name).strip()

    if organism_name == "":
        return None

    organism_name = re.sub(r"\s+", " ", organism_name)

    match = re.search(r"\bvirus\b", organism_name, flags=re.IGNORECASE)

    if match is not None:
        prefix = organism_name[:match.start()].strip()
        organism_key = f"{prefix} virus".strip()
    else:
        organism_key = organism_name.split("(", 1)[0].strip()

    organism_key = re.sub(r"\s+", " ", organism_key)

    if organism_key == "":
        return None

    return organism_key

def infer_single_organism_from_metadata(metadata_rows):
    """
    Infer the operational organism key from parsed metadata rows.

    This assumes that one ETL run corresponds to one organism after downstream
    organism-key normalisation.

    Args:
        metadata_rows (list[dict]):
            Parsed metadata rows.

    Returns:
        str:
            Operational organism key.

    Raises:
        ValueError:
            If no organism can be inferred, or if multiple normalised organisms
            are found.
    """
    organism_key_to_raw_values = {}

    for row in metadata_rows:
        raw_organism_name = row.get("organism_name")
        organism_key = extract_organism_key(raw_organism_name)

        if organism_key is None:
            continue

        if organism_key not in organism_key_to_raw_values:
            organism_key_to_raw_values[organism_key] = set()

        organism_key_to_raw_values[organism_key].add(raw_organism_name)

    if len(organism_key_to_raw_values) == 0:
        raise ValueError(
            "Could not infer organism from metadata: no valid organism_name found."
        )

    if len(organism_key_to_raw_values) > 1:
        details = {
            key: sorted(values)
            for key, values in organism_key_to_raw_values.items()
        }

        raise ValueError(
            "Multiple organisms found in metadata after organism-key extraction. "
            f"Found: {details}. "
            "This ETL currently expects one organism per run."
        )

    return next(iter(organism_key_to_raw_values))

def parse_metadata_row(raw_row, line_number=None):
    """
    Convert one raw CSV row into a dictionary compatible with the metadata table.
    """
    row = {}

    for source_col, target_col in METADATA_COLUMN_MAP.items():
        row[target_col] = clean_empty(raw_row.get(source_col))

    if row["accession"] is None:
        raise ValueError(f"Missing ACCESSION at line {line_number}")

    row["length"] = parse_integer(
        row["length"],
        field_name=f"LENGTH at line {line_number}"
    )

    row["segment"] = parse_segment(
        value=row["segment"],
        species=row["species"],
        organism_name=row["organism_name"]
    )

    ha_subtype, na_subtype = parse_subtypes(
        genotype=row["genotype"],
        species=row["species"]
    )
    row["ha_subtype"] = ha_subtype
    row["na_subtype"] = na_subtype

    return row


def parse_metadata_csv(metadata_csv_path, delimiter=";"):
    """
    Parse metadata CSV/semicolon-separated file into ETL-ready dictionaries.

    Parameters
    ----------
    metadata_csv_path : str
        Path to the metadata file.

    delimiter : str
        Field delimiter. Your example uses ';'.

    Returns
    -------
    list[dict]
        List of dictionaries ready for insertion into the metadata table.
    """
    metadata_rows = []

    seen_accessions = {}
    accession_to_index = {}

    with open(metadata_csv_path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)

        validate_metadata_header(reader.fieldnames)

        for line_number, raw_row in enumerate(reader, start=2):
            row = parse_metadata_row(raw_row, line_number=line_number)

            accession = row["accession"]

            if accession in seen_accessions:
                prev = seen_accessions[accession]
                if prev == row:
                    continue
                print(
                    f"  Warning: duplicate ACCESSION with different data "
                    f"at line {line_number}: {accession} — keeping last occurrence"
                )
                idx = accession_to_index[accession]
                metadata_rows[idx] = row
                seen_accessions[accession] = row
                continue

            accession_to_index[accession] = len(metadata_rows)
            seen_accessions[accession] = row
            metadata_rows.append(row)

    return metadata_rows

def transform_seqdict(seqdict):
    '''
    Transforms the sequence dict into the form to use for sqlite load

    Each row also carries precomputed ambiguous base counts (n_N,
    n_degenerate) and a JSON per-base breakdown.
    '''
    sequence_rows = []

    for accession, (description, sequence) in seqdict.items():
        counts = count_ambiguous_bases(sequence)

        sequence_rows.append({
            "accession": accession,
            "description": description,
            "sequence": sequence,
            "n_N": counts["n_N"],
            "n_degenerate": counts["n_degenerate"],
            "degenerate_breakdown": degenerate_breakdown_json(sequence)
        })
    return sequence_rows

def transform_clustdict(clustdict):
    '''
    Transforms the cluster dict into the form to use for sqlite load
    '''
    cluster_rows = []

    for organism, clusters_for_organism in clustdict.items():
        for cluster_number, centroid in clusters_for_organism.items():
            cluster_rows.append({
                "organism": organism,
                "cluster_number": cluster_number,
                "centroid": centroid
                })
    return cluster_rows

def transform_compdict(compdict):
    '''
    Transforms the cluster composition dict into the form to use for sqlite load
    '''
    cluster_composition_rows = []

    for organism, members_for_organism in compdict.items():
        for accession, (cluster_number, identity) in members_for_organism.items():
            cluster_composition_rows.append({
                "organism": organism,
                "cluster_number": cluster_number,
                "accession": accession,
                "identity_to_centroid": identity
            })
    return cluster_composition_rows


def deduplicate_rows(rows, label="rows"):
    """
    Remove exact duplicate dictionaries from a list, keeping the first occurrence.

    Parameters
    ----------
    rows : list[dict]
        List of row dictionaries.

    label : str
        Label for log messages.

    Returns
    -------
    list[dict]
        Deduplicated list.
    """
    seen = set()
    deduped = []
    dropped = 0

    for row in rows:
        key = tuple(sorted(row.items()))
        if key in seen:
            dropped += 1
            continue
        seen.add(key)
        deduped.append(row)

    if dropped:
        print(f"  Dropped {dropped} exact duplicate {label}")

    return deduped


# DERIVED CLUSTER TABLE BUILDERS

def clean_derived_value(value, unknown_value="Unknown"):
    """
    Normalise values for derived aggregate tables.

    Base tables keep missing values as SQL NULL.
    Derived tables use a controlled string for grouping/key stability,
    especially in composite PRIMARY KEY columns.
    """
    if value is None:
        return unknown_value

    value = str(value).strip()

    if value == "":
        return unknown_value

    return value


def clean_optional_value(value):
    """
    Return None for empty values, otherwise stripped string.

    Used for descriptor dominant values where NULL is acceptable.
    """
    if value is None:
        return None

    value = str(value).strip()

    if value == "":
        return None

    return value


def percentile(sorted_values, q):
    """
    Compute a percentile using linear interpolation.

    Parameters
    ----------
    sorted_values : list[float]
        Sorted numeric values.

    q : float
        Percentile between 0.0 and 1.0.

    Returns
    -------
    float or None
        Percentile value.
    """
    if not sorted_values:
        return None

    if len(sorted_values) == 1:
        return sorted_values[0]

    position = (len(sorted_values) - 1) * q
    lower_index = int(position)
    upper_index = min(lower_index + 1, len(sorted_values) - 1)

    lower_value = sorted_values[lower_index]
    upper_value = sorted_values[upper_index]

    fraction = position - lower_index

    return lower_value + ((upper_value - lower_value) * fraction)


def build_lookup_tables(metadata_rows, cluster_rows, cluster_composition_rows):
    """
    Build lookup dictionaries used by derived table generation.

    Returns
    -------
    tuple
        metadata_by_accession,
        centroid_by_cluster,
        members_by_cluster
    """
    metadata_by_accession = {
        row["accession"]: row
        for row in metadata_rows
    }

    centroid_by_cluster = {
        (row["organism"], row["cluster_number"]): row["centroid"]
        for row in cluster_rows
    }

    members_by_cluster = defaultdict(list)

    for row in cluster_composition_rows:
        cluster_key = (row["organism"], row["cluster_number"])
        members_by_cluster[cluster_key].append(row["accession"])

    return metadata_by_accession, centroid_by_cluster, members_by_cluster


def get_counter_from_member_metadata(member_metadata_rows, field):
    """
    Build a Counter for one metadata field from member metadata rows.

    Missing/empty values are ignored for descriptor dictionaries.
    """
    values = []

    for row in member_metadata_rows:
        value = clean_optional_value(row.get(field))

        if value is not None:
            values.append(value)

    return Counter(values)


def dominant_from_counter(counter):
    """
    Return the most common value from a Counter.

    Ties are resolved by Counter.most_common() order.
    """
    if not counter:
        return None

    return counter.most_common(1)[0][0]


def infer_cluster_segment(segment_counter):
    """
    Infer the cluster-level segment label.

    Returns
    -------
    str or None
        - the unique segment if exactly one segment exists
        - 'Mixed' if multiple segments exist
        - None if no segment exists
    """
    if len(segment_counter) == 0:
        return None

    if len(segment_counter) == 1:
        return next(iter(segment_counter.keys()))

    return "Mixed"


def build_cluster_summary_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build rows for cluster_summary.

    One row per complete cluster.
    """
    metadata_by_accession, centroid_by_cluster, members_by_cluster = build_lookup_tables(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    summary_rows = []

    for cluster_key, accessions in members_by_cluster.items():
        organism, cluster_number = cluster_key
        centroid = centroid_by_cluster[cluster_key]

        member_metadata_rows = [
            metadata_by_accession[accession]
            for accession in accessions
            if accession in metadata_by_accession
        ]

        segment_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "segment"
        )

        host_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "host"
        )

        country_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "country"
        )

        genotype_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "genotype"
        )

        lengths = [
            row.get("length")
            for row in member_metadata_rows
            if row.get("length") is not None
        ]

        lengths = sorted(float(length) for length in lengths)

        if lengths:
            length_min = min(lengths)
            length_q1 = percentile(lengths, 0.25)
            length_mean = statistics.mean(lengths)
            length_median = statistics.median(lengths)
            length_q3 = percentile(lengths, 0.75)
            length_max = max(lengths)

            if len(lengths) > 1:
                length_std = statistics.stdev(lengths)
            else:
                length_std = None
        else:
            length_min = None
            length_q1 = None
            length_mean = None
            length_median = None
            length_q3 = None
            length_max = None
            length_std = None

        summary_rows.append({
            "organism": organism,
            "cluster_number": cluster_number,
            "centroid": centroid,
            "segment": infer_cluster_segment(segment_counter),
            "n_segments": len(segment_counter),
            "n_sequences": len(set(accessions)),
            "n_hosts": len(host_counter),
            "n_countries": len(country_counter),
            "n_genotypes": len(genotype_counter),
            "length_min": length_min,
            "length_q1": length_q1,
            "length_mean": length_mean,
            "length_median": length_median,
            "length_q3": length_q3,
            "length_max": length_max,
            "length_std": length_std,
        })

    return summary_rows


def build_cluster_background_count_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build rows for cluster_background_counts.

    Long/tidy format:
    organism, cluster_number, centroid, segment,
    background_type, background_value, n_sequences.
    """
    metadata_by_accession, centroid_by_cluster, members_by_cluster = build_lookup_tables(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    aggregate = Counter()

    background_fields = [
        ("host", "host"),
        ("country", "country"),
        ("genotype", "genotype"),
    ]

    for cluster_key, accessions in members_by_cluster.items():
        organism, cluster_number = cluster_key
        centroid = centroid_by_cluster[cluster_key]

        for accession in accessions:
            metadata = metadata_by_accession.get(accession)

            if metadata is None:
                continue

            segment = clean_derived_value(metadata.get("segment"))

            for field_name, background_type in background_fields:
                background_value = clean_optional_value(metadata.get(field_name))

                if background_value is None:
                    continue

                key = (
                    organism,
                    cluster_number,
                    centroid,
                    segment,
                    background_type,
                    background_value
                )

                aggregate[key] += 1

    rows = []

    for (
        organism,
        cluster_number,
        centroid,
        segment,
        background_type,
        background_value
    ), n_sequences in aggregate.items():

        rows.append({
            "organism": organism,
            "cluster_number": cluster_number,
            "centroid": centroid,
            "segment": segment,
            "background_type": background_type,
            "background_value": background_value,
            "n_sequences": n_sequences,
        })

    return rows


def build_cluster_descriptor_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build rows for cluster_descriptors.

    One row per complete cluster.
    JSON fields store value-count dictionaries.
    """
    metadata_by_accession, centroid_by_cluster, members_by_cluster = build_lookup_tables(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    descriptor_rows = []

    for cluster_key, accessions in members_by_cluster.items():
        organism, cluster_number = cluster_key
        representative = centroid_by_cluster[cluster_key]

        member_metadata_rows = [
            metadata_by_accession[accession]
            for accession in accessions
            if accession in metadata_by_accession
        ]

        genotype_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "genotype"
        )

        segment_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "segment"
        )

        host_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "host"
        )

        country_counter = get_counter_from_member_metadata(
            member_metadata_rows,
            "country"
        )

        lengths = [
            row.get("length")
            for row in member_metadata_rows
            if row.get("length") is not None
        ]

        lengths = sorted(float(length) for length in lengths)

        if lengths:
            length_min = min(lengths)
            length_mean = statistics.mean(lengths)
            length_median = statistics.median(lengths)
            length_max = max(lengths)

            if len(lengths) > 1:
                length_std = statistics.stdev(lengths)
            else:
                length_std = None
        else:
            length_min = None
            length_mean = None
            length_median = None
            length_max = None
            length_std = None

        n_sequences = len(set(accessions))
        n_genotypes = len(genotype_counter)
        n_segments = len(segment_counter)
        n_hosts = len(host_counter)
        n_countries = len(country_counter)

        descriptor_rows.append({
            "organism": organism,
            "cluster_number": cluster_number,
            "representative": representative,
            "n_sequences": n_sequences,
            "genotypes_json": json.dumps(dict(genotype_counter), ensure_ascii=False),
            "segments_json": json.dumps(dict(segment_counter), ensure_ascii=False),
            "hosts_json": json.dumps(dict(host_counter), ensure_ascii=False),
            "countries_json": json.dumps(dict(country_counter), ensure_ascii=False),
            "dominant_genotype": dominant_from_counter(genotype_counter),
            "dominant_segment": dominant_from_counter(segment_counter),
            "dominant_host": dominant_from_counter(host_counter),
            "dominant_country": dominant_from_counter(country_counter),
            "n_genotypes": n_genotypes,
            "n_segments": n_segments,
            "n_hosts": n_hosts,
            "n_countries": n_countries,
            "is_singleton": 1 if n_sequences == 1 else 0,
            "is_mixed_genotype": 1 if n_genotypes > 1 else 0,
            "is_mixed_segment": 1 if n_segments > 1 else 0,
            "is_multi_host": 1 if n_hosts > 1 else 0,
            "is_multi_country": 1 if n_countries > 1 else 0,
            "length_min": length_min,
            "length_mean": length_mean,
            "length_median": length_median,
            "length_max": length_max,
            "length_std": length_std,
        })

    return descriptor_rows


def build_cluster_filter_count_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build rows for cluster_filter_counts.

    One row per observed combination:
    organism, cluster_number, centroid, segment, genotype, host, country.
    """
    metadata_by_accession, centroid_by_cluster, members_by_cluster = build_lookup_tables(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    aggregate = Counter()

    for cluster_key, accessions in members_by_cluster.items():
        organism, cluster_number = cluster_key
        centroid = centroid_by_cluster[cluster_key]

        for accession in accessions:
            metadata = metadata_by_accession.get(accession)

            if metadata is None:
                continue

            segment = clean_derived_value(metadata.get("segment"))
            genotype = clean_derived_value(metadata.get("genotype"))
            ha_subtype = clean_derived_value(metadata.get("ha_subtype"))
            na_subtype = clean_derived_value(metadata.get("na_subtype"))
            host = clean_derived_value(metadata.get("host"))
            country = clean_derived_value(metadata.get("country"))

            key = (
                organism,
                cluster_number,
                centroid,
                segment,
                genotype,
                ha_subtype,
                na_subtype,
                host,
                country
            )

            aggregate[key] += 1

    rows = []

    for (
        organism,
        cluster_number,
        centroid,
        segment,
        genotype,
        ha_subtype,
        na_subtype,
        host,
        country
    ), n_sequences in aggregate.items():

        rows.append({
            "organism": organism,
            "cluster_number": cluster_number,
            "centroid": centroid,
            "segment": segment,
            "genotype": genotype,
            "ha_subtype": ha_subtype,
            "na_subtype": na_subtype,
            "host": host,
            "country": country,
            "n_sequences": n_sequences,
        })

    return rows


def build_cluster_assembly_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build rows for cluster_assembly.

    One row per observed combination:
    organism, cluster_number, centroid, assembly, segment.

    Every sequence must belong to an assembly; sequences whose metadata
    assembly is missing are skipped here and reported by
    compute_assembly_gaps().
    """
    metadata_by_accession, centroid_by_cluster, members_by_cluster = build_lookup_tables(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    aggregate = Counter()
    centroid_assembly_keys = set()

    for cluster_key, accessions in members_by_cluster.items():
        organism, cluster_number = cluster_key
        centroid = centroid_by_cluster[cluster_key]

        for accession in accessions:
            metadata = metadata_by_accession.get(accession)

            if metadata is None:
                continue

            assembly = clean_derived_value(metadata.get("assembly"))
            segment = clean_derived_value(metadata.get("segment"))

            key = (
                organism,
                cluster_number,
                centroid,
                assembly,
                segment
            )

            aggregate[key] += 1

            if accession == centroid:
                centroid_assembly_keys.add(key)

    rows = []

    for (
        organism,
        cluster_number,
        centroid,
        assembly,
        segment
    ), n_sequences in aggregate.items():

        key = (organism, cluster_number, centroid, assembly, segment)

        rows.append({
            "organism": organism,
            "cluster_number": cluster_number,
            "centroid": centroid,
            "assembly": assembly,
            "segment": segment,
            "n_sequences": n_sequences,
            "is_centroid_assembly": 1 if key in centroid_assembly_keys else 0,
        })

    return rows


def build_cluster_assembly_link_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build rows for cluster_assembly_links.

    One row per unique cluster-cluster link, where n_connections is the
    number of assemblies shared by both clusters.

    Clusters are linked when they hold sequences from a common assembly.
    Because an assembly contributes at most one cluster per segment, and
    because each assembly spans one cluster per segment, an assembly
    links every pair of its clusters together.

    Each link is stored once with the lexicographically smaller cluster
    key first, so the primary key is symmetric-safe. Self-links are
    skipped.
    """
    cluster_assembly_rows = build_cluster_assembly_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    segment_counter_by_cluster = defaultdict(Counter)

    for row in cluster_assembly_rows:
        cluster_key = (row["organism"], row["cluster_number"])

        segment_counter_by_cluster[cluster_key][row["segment"]] += (
            row["n_sequences"]
        )

    # Mirror cluster_summary.segment so the two tables never disagree.
    cluster_segment_by_key = {
        cluster_key: infer_cluster_segment(segment_counter)
        for cluster_key, segment_counter in segment_counter_by_cluster.items()
    }

    clusters_by_assembly = defaultdict(set)

    for row in cluster_assembly_rows:
        cluster_key = (row["organism"], row["cluster_number"])

        clusters_by_assembly[row["assembly"]].add(cluster_key)

    link_counts = Counter()

    for cluster_keys in clusters_by_assembly.values():
        if len(cluster_keys) < 2:
            continue

        sorted_keys = sorted(cluster_keys)

        for i in range(len(sorted_keys)):
            first_key = sorted_keys[i]

            for j in range(i + 1, len(sorted_keys)):
                link_counts[(first_key, sorted_keys[j])] += 1

    rows = []

    for (first_key, second_key), n_connections in link_counts.items():
        rows.append({
            "organism": first_key[0],
            "cluster_number": first_key[1],
            "linked_organism": second_key[0],
            "linked_cluster_number": second_key[1],
            "linked_cluster_segment": cluster_segment_by_key.get(
                second_key,
                "Unknown"
            ),
            "n_connections": n_connections,
        })

    return rows


def build_derived_cluster_rows(
    metadata_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Build all derived cluster tables in memory.

    Returns
    -------
    dict
        {
            "cluster_summary": [...],
            "cluster_background_counts": [...],
            "cluster_descriptors": [...],
            "cluster_filter_counts": [...],
            "cluster_assembly": [...],
            "cluster_assembly_links": [...]
        }
    """
    cluster_summary_rows = build_cluster_summary_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    cluster_background_count_rows = build_cluster_background_count_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    cluster_descriptor_rows = build_cluster_descriptor_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    cluster_filter_count_rows = build_cluster_filter_count_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    cluster_assembly_rows = build_cluster_assembly_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    cluster_assembly_link_rows = build_cluster_assembly_link_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=cluster_composition_rows
    )

    return {
        "cluster_summary": cluster_summary_rows,
        "cluster_background_counts": cluster_background_count_rows,
        "cluster_descriptors": cluster_descriptor_rows,
        "cluster_filter_counts": cluster_filter_count_rows,
        "cluster_assembly": cluster_assembly_rows,
        "cluster_assembly_links": cluster_assembly_link_rows,
    }


def validate_cross_references(
    metadata_rows,
    sequence_rows,
    cluster_rows,
    cluster_composition_rows
):
    """
    Validate cross-table references before loading into SQLite.

    This gives clearer errors than waiting for SQLite foreign key failures.
    """
    metadata_accessions = {
        row["accession"]
        for row in metadata_rows
    }

    sequence_accessions = {
        row["accession"]
        for row in sequence_rows
    }

    cluster_keys = {
        (row["organism"], row["cluster_number"])
        for row in cluster_rows
    }

    missing_sequence_accessions = sorted(sequence_accessions - metadata_accessions)

    if missing_sequence_accessions:
        raise ValueError(
            "Some FASTA sequence accessions are missing from metadata. "
            f"First examples: {missing_sequence_accessions[:10]}"
        )

    centroid_accessions = {
        row["centroid"]
        for row in cluster_rows
    }

    missing_centroids = sorted(centroid_accessions - metadata_accessions)

    if missing_centroids:
        raise ValueError(
            "Some cluster centroids are missing from metadata. "
            f"First examples: {missing_centroids[:10]}"
        )

    composition_accessions = {
        row["accession"]
        for row in cluster_composition_rows
    }

    missing_composition_accessions = sorted(
        composition_accessions - metadata_accessions
    )

    if missing_composition_accessions:
        raise ValueError(
            "Some cluster composition accessions are missing from metadata. "
            f"First examples: {missing_composition_accessions[:10]}"
        )

    composition_cluster_keys = {
        (row["organism"], row["cluster_number"])
        for row in cluster_composition_rows
    }

    missing_cluster_keys = sorted(composition_cluster_keys - cluster_keys)

    if missing_cluster_keys:
        raise ValueError(
            "Some cluster composition rows reference missing clusters. "
            f"First examples: {missing_cluster_keys[:10]}"
        )


def species_key_for_row(row):
    """
    Derive the operational species/organism key for one metadata row.

    The three data sources do not share a single species column:
    - metadata has organism_name and species;
    - cluster data has an organism key;
    - FASTA sequences have neither.

    To compare accessions "of the same species" we normalise every table to
    the operational organism key extracted from organism_name, falling back
    to species and finally to 'Unknown'.
    """
    organism_key = extract_organism_key(row.get("organism_name"))

    if organism_key is not None:
        return organism_key

    species = row.get("species")

    if species is None:
        return "Unknown"

    species = str(species).strip()

    return species if species else "Unknown"


def compute_record_gaps(
    metadata_rows,
    sequence_rows,
    cluster_composition_rows,
    organism=None,
):
    """
    Find accessions that are not present in all source tables, per species.

    A complete record must appear in three places:
    - metadata (the metadata table);
    - cluster composition (the .clstr file);
    - FASTA (the sequence file).

    Returns
    -------
    dict
        {
            species_key: {
                "missing_from_metadata": [accession, ...],
                "missing_from_fasta": [accession, ...],
                "missing_from_composition": [accession, ...]
            }
        }
    """
    by_species = defaultdict(lambda: defaultdict(set))

    for row in metadata_rows:
        key = species_key_for_row(row)
        by_species[key]["metadata"].add(row["accession"])

    for row in sequence_rows:
        key = organism if organism is not None else "Unknown"
        by_species[key]["fasta"].add(row["accession"])

    for row in cluster_composition_rows:
        key = row["organism"] if row.get("organism") else "Unknown"
        by_species[key]["composition"].add(row["accession"])

    gaps = {}

    for key, sources in sorted(by_species.items()):
        metadata_set = sources.get("metadata", set())
        fasta_set = sources.get("fasta", set())
        composition_set = sources.get("composition", set())

        present_somewhere = metadata_set | fasta_set | composition_set

        missing_from_metadata = sorted(present_somewhere - metadata_set)
        missing_from_fasta = sorted(present_somewhere - fasta_set)
        missing_from_composition = sorted(present_somewhere - composition_set)

        if any([missing_from_metadata, missing_from_fasta, missing_from_composition]):
            gaps[key] = {
                "missing_from_metadata": missing_from_metadata,
                "missing_from_fasta": missing_from_fasta,
                "missing_from_composition": missing_from_composition,
            }

    return gaps


def report_record_gaps(gaps):
    """
    Print a human-readable summary of record completeness gaps.
    """
    total_missing_fasta = sum(
        len(details["missing_from_fasta"])
        for details in gaps.values()
    )

    if not gaps:
        print("Record completeness: every sequence has metadata, a cluster "
              "composition entry and a FASTA sequence.")
        return

    print("Record completeness gaps")
    print("------------------------")
    print(f"Accessions missing a FASTA sequence: {total_missing_fasta}")

    for key, details in gaps.items():
        print(f"  Species '{key}':")
        for field, label in [
            ("missing_from_metadata", "missing from metadata"),
            ("missing_from_fasta", "missing from FASTA"),
            ("missing_from_composition", "missing from cluster composition"),
        ]:
            accessions = details[field]
            if not accessions:
                continue
            preview = ", ".join(accessions[:10])
            extra = f" ... (+{len(accessions) - 10} more)" if len(accessions) > 10 else ""
            print(f"    {label} ({len(accessions)}): {preview}{extra}")


def compute_assembly_gaps(metadata_rows, cluster_composition_rows):
    """
    Find cluster members whose metadata is missing an assembly.

    Every sequence must be part of an assembly, because assembly is the
    association point between clusters and genomes. Sequences with a
    missing assembly are dropped from cluster_assembly by the ETL, so
    they would silently disappear from assembly-filtered results.

    Returns
    -------
    dict
        {
            species_key: {
                "missing_assembly": [accession, ...]
            }
        }
    """
    assembly_by_accession = {}

    for row in metadata_rows:
        assembly = row.get("assembly")

        if assembly is None:
            continue

        assembly = str(assembly).strip()

        if assembly:
            assembly_by_accession[row["accession"]] = assembly

    gaps = defaultdict(list)

    for row in cluster_composition_rows:
        accession = row["accession"]

        if accession in assembly_by_accession:
            continue

        gaps[species_key_for_row({"species": row.get("organism")})].append(
            accession
        )

    return {
        key: {"missing_assembly": sorted(set(accessions))}
        for key, accessions in sorted(gaps.items())
    }


def report_assembly_gaps(gaps):
    """
    Print a human-readable summary of missing-assembly sequences.
    """
    total = sum(
        len(details["missing_assembly"])
        for details in gaps.values()
    )

    if not gaps:
        print("Assembly completeness: every cluster member belongs to an assembly.")
        return

    print("Assembly completeness gaps")
    print("--------------------------")
    print(f"Sequences missing an assembly: {total}")

    for key, details in gaps.items():
        accessions = details["missing_assembly"]
        preview = ", ".join(accessions[:10])
        extra = f" ... (+{len(accessions) - 10} more)" if len(accessions) > 10 else ""
        print(f"  Species '{key}':")
        print(f"    missing assembly ({len(accessions)}): {preview}{extra}")


def missing_fasta_accessions(gaps):
    """
    Flatten the set of accessions missing a FASTA sequence across species.

    This is the set the NCBI fetch step needs to recover.
    """
    missing = set()

    for details in gaps.values():
        missing.update(details["missing_from_fasta"])

    return sorted(missing)


def validate_species_record_completeness(
    metadata_rows,
    sequence_rows,
    cluster_composition_rows,
    organism=None,
    mode="warn",
):
    """
    Validate that every sequence is present in metadata, cluster composition
    and the FASTA file, grouped by species.

    Parameters
    ----------
    mode : str
        'warn' prints the report and returns; 'fail' raises ValueError if
        any record is incomplete.
    """
    gaps = compute_record_gaps(
        metadata_rows=metadata_rows,
        sequence_rows=sequence_rows,
        cluster_composition_rows=cluster_composition_rows,
        organism=organism,
    )

    report_record_gaps(gaps)

    if mode == "fail" and gaps:
        first = next(iter(gaps.values()))
        raise ValueError(
            "Record completeness check failed: some sequences are not present "
            "in all tables of the same species. "
            "Use --fetch-missing to download missing sequences from NCBI, or "
            "run with --record-completeness warn to continue. "
            f"First missing FASTA examples: {first['missing_from_fasta'][:10]}"
        )

    assembly_gaps = compute_assembly_gaps(
        metadata_rows=metadata_rows,
        cluster_composition_rows=cluster_composition_rows,
    )

    report_assembly_gaps(assembly_gaps)

    if mode == "fail" and assembly_gaps:
        first = next(iter(assembly_gaps.values()))
        raise ValueError(
            "Assembly completeness check failed: some sequences do not belong "
            "to an assembly, so they are excluded from cluster_assembly. "
            f"First missing assembly examples: {first['missing_assembly'][:10]}"
        )

    return gaps

############### SQLite ########################
def connect_sqlite(db_path):
    """
    Open a SQLite connection and enable foreign key checking.

    Args:
        db_path (str):
            Path to the SQLite database file.

    Returns:
        sqlite3.Connection:
            Open SQLite connection.
    """
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON;")
    return conn

def initialise_database(db_path, schema_path, recreate=False):
    """
    Initialise the SQLite database from a schema.sql file.

    Args:
        db_path (str or Path):
            Path to the SQLite database file.

        schema_path (str or Path):
            Path to the schema.sql file.

        recreate (bool):
            If True, delete the existing SQLite database before creating it.

    Returns:
        None
    """
    db_path = Path(db_path)
    schema_path = Path(schema_path)

    db_path.parent.mkdir(parents=True, exist_ok=True)

    if recreate and db_path.exists():
        db_path.unlink()

    if not schema_path.exists():
        raise FileNotFoundError(f"Schema file not found: {schema_path}")

    conn = connect_sqlite(db_path)

    try:
        schema_sql = schema_path.read_text(encoding="utf-8")
        conn.executescript(schema_sql)
        conn.commit()

    finally:
        conn.close()

###### VERSOES BÁSICAS, COM COMMIT ############
def load_metadata(conn, metadata_rows):
    """
    Load rows into the metadata table.

    Expected row structure
    ----------------------
    {
        "accession": ...,
        "organism_name": ...,
        "genbank_refseq": ...,
        "assembly": ...,
        "release_date": ...,
        "isolate": ...,
        "species": ...,
        "length": ...,
        "nuc_completeness": ...,
        "genotype": ...,
        "segment": ...,
        "country": ...,
        "host": ...,
        "collection_date": ...,
        "comment": ...
    }

    Args:
        conn (sqlite3.Connection):
            Open SQLite connection.

        metadata_rows (list[dict]):
            Metadata rows ready for SQLite insertion.

    Returns:
        int:
            Number of rows submitted.
    """
    if not metadata_rows:
        return 0

    sql = """
    INSERT INTO metadata (
        accession,
        organism_name,
        genbank_refseq,
        assembly,
        release_date,
        isolate,
        species,
        length,
        nuc_completeness,
        genotype,
        ha_subtype,
        na_subtype,
        segment,
        country,
        host,
        collection_date,
        comment
    )
    VALUES (
        :accession,
        :organism_name,
        :genbank_refseq,
        :assembly,
        :release_date,
        :isolate,
        :species,
        :length,
        :nuc_completeness,
        :genotype,
        :ha_subtype,
        :na_subtype,
        :segment,
        :country,
        :host,
        :collection_date,
        :comment
    )
    ON CONFLICT(accession) DO UPDATE SET
        organism_name = excluded.organism_name,
        genbank_refseq = excluded.genbank_refseq,
        assembly = excluded.assembly,
        release_date = excluded.release_date,
        isolate = excluded.isolate,
        species = excluded.species,
        length = excluded.length,
        nuc_completeness = excluded.nuc_completeness,
        genotype = excluded.genotype,
        ha_subtype = excluded.ha_subtype,
        na_subtype = excluded.na_subtype,
        segment = excluded.segment,
        country = excluded.country,
        host = excluded.host,
        collection_date = excluded.collection_date,
        comment = excluded.comment;
    """

    conn.executemany(sql, metadata_rows)
    conn.commit()

    return len(metadata_rows)

def load_sequences(conn, sequence_rows):
    """
    Load rows into the sequences table.

    Expected row structure
    ----------------------
    {
        "accession": ...,
        "description": ...,
        "sequence": ...
    }

    Args:
        conn (sqlite3.Connection):
            Open SQLite connection.

        sequence_rows (list[dict]):
            Sequence rows ready for SQLite insertion.

    Returns:
        int:
            Number of rows submitted.
    """
    if not sequence_rows:
        return 0

    sql = """
    INSERT INTO sequences (
        accession,
        description,
        sequence,
        n_N,
        n_degenerate,
        degenerate_breakdown
    )
    VALUES (
        :accession,
        :description,
        :sequence,
        :n_N,
        :n_degenerate,
        :degenerate_breakdown
    )
    ON CONFLICT(accession) DO UPDATE SET
        description = excluded.description,
        sequence = excluded.sequence,
        n_N = excluded.n_N,
        n_degenerate = excluded.n_degenerate,
        degenerate_breakdown = excluded.degenerate_breakdown;
    """

    conn.executemany(sql, sequence_rows)
    conn.commit()

    return len(sequence_rows)

def load_clusters(conn, cluster_rows):
    """
    Load rows into the clusters table.

    Expected row structure
    ----------------------
    {
        "organism": ...,
        "cluster_number": ...,
        "centroid": ...
    }

    Args:
        conn (sqlite3.Connection):
            Open SQLite connection.

        cluster_rows (list[dict]):
            Cluster rows ready for SQLite insertion.

    Returns:
        int:
            Number of rows submitted.
    """
    if not cluster_rows:
        return 0

    sql = """
    INSERT INTO clusters (
        organism,
        cluster_number,
        centroid
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid
    )
    ON CONFLICT(organism, cluster_number) DO UPDATE SET
        centroid = excluded.centroid;
    """

    conn.executemany(sql, cluster_rows)
    conn.commit()

    return len(cluster_rows)

def load_cluster_composition(conn, cluster_composition_rows):
    """
    Load rows into the cluster_composition table.

    Expected row structure
    ----------------------
    {
        "organism": ...,
        "cluster_number": ...,
        "accession": ...,
        "identity_to_centroid": ...
    }

    Args:
        conn (sqlite3.Connection):
            Open SQLite connection.

        cluster_composition_rows (list[dict]):
            Cluster composition rows ready for SQLite insertion.

    Returns:
        int:
            Number of rows submitted.
    """
    if not cluster_composition_rows:
        return 0

    sql = """
    INSERT INTO cluster_composition (
        organism,
        cluster_number,
        accession,
        identity_to_centroid
    )
    VALUES (
        :organism,
        :cluster_number,
        :accession,
        :identity_to_centroid
    )
    ON CONFLICT(accession) DO UPDATE SET
        organism = excluded.organism,
        cluster_number = excluded.cluster_number,
        identity_to_centroid = excluded.identity_to_centroid;
    """

    conn.executemany(sql, cluster_composition_rows)
    conn.commit()

    return len(cluster_composition_rows)

def load_cluster_summary(conn, cluster_summary_rows):
    """
    Load rows into cluster_summary.
    """
    if not cluster_summary_rows:
        return 0

    sql = """
    INSERT INTO cluster_summary (
        organism,
        cluster_number,
        centroid,
        segment,
        n_segments,
        n_sequences,
        n_hosts,
        n_countries,
        n_genotypes,
        length_min,
        length_q1,
        length_mean,
        length_median,
        length_q3,
        length_max,
        length_std
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :segment,
        :n_segments,
        :n_sequences,
        :n_hosts,
        :n_countries,
        :n_genotypes,
        :length_min,
        :length_q1,
        :length_mean,
        :length_median,
        :length_q3,
        :length_max,
        :length_std
    )
    ON CONFLICT(organism, cluster_number) DO UPDATE SET
        centroid = excluded.centroid,
        segment = excluded.segment,
        n_segments = excluded.n_segments,
        n_sequences = excluded.n_sequences,
        n_hosts = excluded.n_hosts,
        n_countries = excluded.n_countries,
        n_genotypes = excluded.n_genotypes,
        length_min = excluded.length_min,
        length_q1 = excluded.length_q1,
        length_mean = excluded.length_mean,
        length_median = excluded.length_median,
        length_q3 = excluded.length_q3,
        length_max = excluded.length_max,
        length_std = excluded.length_std;
    """

    conn.executemany(sql, cluster_summary_rows)
    conn.commit()

    return len(cluster_summary_rows)


def load_cluster_background_counts(conn, cluster_background_count_rows):
    """
    Load rows into cluster_background_counts.
    """
    if not cluster_background_count_rows:
        return 0

    sql = """
    INSERT INTO cluster_background_counts (
        organism,
        cluster_number,
        centroid,
        segment,
        background_type,
        background_value,
        n_sequences
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :segment,
        :background_type,
        :background_value,
        :n_sequences
    )
    ON CONFLICT(
        organism,
        cluster_number,
        segment,
        background_type,
        background_value
    ) DO UPDATE SET
        centroid = excluded.centroid,
        n_sequences = excluded.n_sequences;
    """

    conn.executemany(sql, cluster_background_count_rows)
    conn.commit()

    return len(cluster_background_count_rows)


def load_cluster_descriptors(conn, cluster_descriptor_rows):
    """
    Load rows into cluster_descriptors.
    """
    if not cluster_descriptor_rows:
        return 0

    sql = """
    INSERT INTO cluster_descriptors (
        organism,
        cluster_number,
        representative,
        n_sequences,
        genotypes_json,
        segments_json,
        hosts_json,
        countries_json,
        dominant_genotype,
        dominant_segment,
        dominant_host,
        dominant_country,
        n_genotypes,
        n_segments,
        n_hosts,
        n_countries,
        is_singleton,
        is_mixed_genotype,
        is_mixed_segment,
        is_multi_host,
        is_multi_country,
        length_min,
        length_mean,
        length_median,
        length_max,
        length_std
    )
    VALUES (
        :organism,
        :cluster_number,
        :representative,
        :n_sequences,
        :genotypes_json,
        :segments_json,
        :hosts_json,
        :countries_json,
        :dominant_genotype,
        :dominant_segment,
        :dominant_host,
        :dominant_country,
        :n_genotypes,
        :n_segments,
        :n_hosts,
        :n_countries,
        :is_singleton,
        :is_mixed_genotype,
        :is_mixed_segment,
        :is_multi_host,
        :is_multi_country,
        :length_min,
        :length_mean,
        :length_median,
        :length_max,
        :length_std
    )
    ON CONFLICT(organism, cluster_number) DO UPDATE SET
        representative = excluded.representative,
        n_sequences = excluded.n_sequences,
        genotypes_json = excluded.genotypes_json,
        segments_json = excluded.segments_json,
        hosts_json = excluded.hosts_json,
        countries_json = excluded.countries_json,
        dominant_genotype = excluded.dominant_genotype,
        dominant_segment = excluded.dominant_segment,
        dominant_host = excluded.dominant_host,
        dominant_country = excluded.dominant_country,
        n_genotypes = excluded.n_genotypes,
        n_segments = excluded.n_segments,
        n_hosts = excluded.n_hosts,
        n_countries = excluded.n_countries,
        is_singleton = excluded.is_singleton,
        is_mixed_genotype = excluded.is_mixed_genotype,
        is_mixed_segment = excluded.is_mixed_segment,
        is_multi_host = excluded.is_multi_host,
        is_multi_country = excluded.is_multi_country,
        length_min = excluded.length_min,
        length_mean = excluded.length_mean,
        length_median = excluded.length_median,
        length_max = excluded.length_max,
        length_std = excluded.length_std;
    """

    conn.executemany(sql, cluster_descriptor_rows)
    conn.commit()

    return len(cluster_descriptor_rows)


def load_cluster_filter_counts(conn, cluster_filter_count_rows):
    """
    Load rows into cluster_filter_counts.
    """
    if not cluster_filter_count_rows:
        return 0

    sql = """
    INSERT INTO cluster_filter_counts (
        organism,
        cluster_number,
        centroid,
        segment,
        genotype,
        ha_subtype,
        na_subtype,
        host,
        country,
        n_sequences
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :segment,
        :genotype,
        :ha_subtype,
        :na_subtype,
        :host,
        :country,
        :n_sequences
    )
    ON CONFLICT(
        organism,
        cluster_number,
        segment,
        genotype,
        ha_subtype,
        na_subtype,
        host,
        country
    ) DO UPDATE SET
        centroid = excluded.centroid,
        n_sequences = excluded.n_sequences;
    """

    conn.executemany(sql, cluster_filter_count_rows)
    conn.commit()

    return len(cluster_filter_count_rows)


def load_cluster_assembly(conn, cluster_assembly_rows):
    """
    Load rows into cluster_assembly.
    """
    if not cluster_assembly_rows:
        return 0

    sql = """
    INSERT INTO cluster_assembly (
        organism,
        cluster_number,
        centroid,
        assembly,
        segment,
        n_sequences,
        is_centroid_assembly
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :assembly,
        :segment,
        :n_sequences,
        :is_centroid_assembly
    )
    ON CONFLICT(
        organism,
        cluster_number,
        assembly,
        segment
    ) DO UPDATE SET
        centroid = excluded.centroid,
        n_sequences = excluded.n_sequences,
        is_centroid_assembly = excluded.is_centroid_assembly;
    """

    conn.executemany(sql, cluster_assembly_rows)
    conn.commit()

    return len(cluster_assembly_rows)


def load_cluster_assembly_links(conn, cluster_assembly_link_rows):
    """
    Load rows into cluster_assembly_links.
    """
    if not cluster_assembly_link_rows:
        return 0

    sql = """
    INSERT INTO cluster_assembly_links (
        organism,
        cluster_number,
        linked_organism,
        linked_cluster_number,
        linked_cluster_segment,
        n_connections
    )
    VALUES (
        :organism,
        :cluster_number,
        :linked_organism,
        :linked_cluster_number,
        :linked_cluster_segment,
        :n_connections
    )
    ON CONFLICT(
        organism,
        cluster_number,
        linked_organism,
        linked_cluster_number
    ) DO UPDATE SET
        linked_cluster_segment = excluded.linked_cluster_segment,
        n_connections = excluded.n_connections;
    """

    conn.executemany(sql, cluster_assembly_link_rows)
    conn.commit()

    return len(cluster_assembly_link_rows)


def load_all_tables(
    db_path,
    metadata_rows,
    sequence_rows,
    cluster_rows,
    cluster_composition_rows,
    cluster_summary_rows=None,
    cluster_background_count_rows=None,
    cluster_descriptor_rows=None,
    cluster_filter_count_rows=None,
    cluster_assembly_rows=None,
    cluster_assembly_link_rows=None
):
    """
    Load all ETL-ready rows into the SQLite database.

    Loading order
    -------------
    1. metadata
    2. sequences
    3. clusters
    4. cluster_composition
    5. cluster_summary
    6. cluster_background_counts
    7. cluster_descriptors
    8. cluster_filter_counts
    9. cluster_assembly
    10. cluster_assembly_links
    """
    cluster_summary_rows = cluster_summary_rows or []
    cluster_background_count_rows = cluster_background_count_rows or []
    cluster_descriptor_rows = cluster_descriptor_rows or []
    cluster_filter_count_rows = cluster_filter_count_rows or []
    cluster_assembly_rows = cluster_assembly_rows or []
    cluster_assembly_link_rows = cluster_assembly_link_rows or []

    conn = connect_sqlite(db_path)

    try:
        n_metadata = load_metadata(conn, metadata_rows)
        n_sequences = load_sequences(conn, sequence_rows)
        n_clusters = load_clusters(conn, cluster_rows)
        n_cluster_composition = load_cluster_composition(
            conn,
            cluster_composition_rows
        )

        n_cluster_summary = load_cluster_summary(
            conn,
            cluster_summary_rows
        )

        n_cluster_background_counts = load_cluster_background_counts(
            conn,
            cluster_background_count_rows
        )

        n_cluster_descriptors = load_cluster_descriptors(
            conn,
            cluster_descriptor_rows
        )

        n_cluster_filter_counts = load_cluster_filter_counts(
            conn,
            cluster_filter_count_rows
        )

        n_cluster_assembly = load_cluster_assembly(
            conn,
            cluster_assembly_rows
        )

        n_cluster_assembly_links = load_cluster_assembly_links(
            conn,
            cluster_assembly_link_rows
        )

        summary = {
            "metadata": n_metadata,
            "sequences": n_sequences,
            "clusters": n_clusters,
            "cluster_composition": n_cluster_composition,
            "cluster_summary": n_cluster_summary,
            "cluster_background_counts": n_cluster_background_counts,
            "cluster_descriptors": n_cluster_descriptors,
            "cluster_filter_counts": n_cluster_filter_counts,
            "cluster_assembly": n_cluster_assembly,
            "cluster_assembly_links": n_cluster_assembly_links,
        }

        return summary

    finally:
        conn.close()

###### VERSOES SEGURAS, SEM COMMIT ############
def load_metadata_no_commit(conn, metadata_rows):
    if not metadata_rows:
        return 0

    sql = """
    INSERT INTO metadata (
        accession,
        organism_name,
        genbank_refseq,
        assembly,
        release_date,
        isolate,
        species,
        length,
        nuc_completeness,
        genotype,
        ha_subtype,
        na_subtype,
        segment,
        country,
        host,
        collection_date,
        comment
    )
    VALUES (
        :accession,
        :organism_name,
        :genbank_refseq,
        :assembly,
        :release_date,
        :isolate,
        :species,
        :length,
        :nuc_completeness,
        :genotype,
        :ha_subtype,
        :na_subtype,
        :segment,
        :country,
        :host,
        :collection_date,
        :comment
    )
    ON CONFLICT(accession) DO UPDATE SET
        organism_name = excluded.organism_name,
        genbank_refseq = excluded.genbank_refseq,
        assembly = excluded.assembly,
        release_date = excluded.release_date,
        isolate = excluded.isolate,
        species = excluded.species,
        length = excluded.length,
        nuc_completeness = excluded.nuc_completeness,
        genotype = excluded.genotype,
        ha_subtype = excluded.ha_subtype,
        na_subtype = excluded.na_subtype,
        segment = excluded.segment,
        country = excluded.country,
        host = excluded.host,
        collection_date = excluded.collection_date,
        comment = excluded.comment;
    """

    conn.executemany(sql, metadata_rows)
    return len(metadata_rows)


def load_sequences_no_commit(conn, sequence_rows):
    if not sequence_rows:
        return 0

    sql = """
    INSERT INTO sequences (
        accession,
        description,
        sequence,
        n_N,
        n_degenerate,
        degenerate_breakdown
    )
    VALUES (
        :accession,
        :description,
        :sequence,
        :n_N,
        :n_degenerate,
        :degenerate_breakdown
    )
    ON CONFLICT(accession) DO UPDATE SET
        description = excluded.description,
        sequence = excluded.sequence,
        n_N = excluded.n_N,
        n_degenerate = excluded.n_degenerate,
        degenerate_breakdown = excluded.degenerate_breakdown;
    """

    conn.executemany(sql, sequence_rows)
    return len(sequence_rows)


def load_clusters_no_commit(conn, cluster_rows):
    if not cluster_rows:
        return 0

    sql = """
    INSERT INTO clusters (
        organism,
        cluster_number,
        centroid
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid
    )
    ON CONFLICT(organism, cluster_number) DO UPDATE SET
        centroid = excluded.centroid;
    """

    conn.executemany(sql, cluster_rows)
    return len(cluster_rows)


def load_cluster_composition_no_commit(conn, cluster_composition_rows):
    if not cluster_composition_rows:
        return 0

    sql = """
    INSERT INTO cluster_composition (
        organism,
        cluster_number,
        accession,
        identity_to_centroid
    )
    VALUES (
        :organism,
        :cluster_number,
        :accession,
        :identity_to_centroid
    )
    ON CONFLICT(accession) DO UPDATE SET
        organism = excluded.organism,
        cluster_number = excluded.cluster_number,
        identity_to_centroid = excluded.identity_to_centroid;
    """

    conn.executemany(sql, cluster_composition_rows)
    return len(cluster_composition_rows)

def load_cluster_summary_no_commit(conn, cluster_summary_rows):
    if not cluster_summary_rows:
        return 0

    sql = """
    INSERT INTO cluster_summary (
        organism,
        cluster_number,
        centroid,
        segment,
        n_segments,
        n_sequences,
        n_hosts,
        n_countries,
        n_genotypes,
        length_min,
        length_q1,
        length_mean,
        length_median,
        length_q3,
        length_max,
        length_std
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :segment,
        :n_segments,
        :n_sequences,
        :n_hosts,
        :n_countries,
        :n_genotypes,
        :length_min,
        :length_q1,
        :length_mean,
        :length_median,
        :length_q3,
        :length_max,
        :length_std
    )
    ON CONFLICT(organism, cluster_number) DO UPDATE SET
        centroid = excluded.centroid,
        segment = excluded.segment,
        n_segments = excluded.n_segments,
        n_sequences = excluded.n_sequences,
        n_hosts = excluded.n_hosts,
        n_countries = excluded.n_countries,
        n_genotypes = excluded.n_genotypes,
        length_min = excluded.length_min,
        length_q1 = excluded.length_q1,
        length_mean = excluded.length_mean,
        length_median = excluded.length_median,
        length_q3 = excluded.length_q3,
        length_max = excluded.length_max,
        length_std = excluded.length_std;
    """

    conn.executemany(sql, cluster_summary_rows)
    return len(cluster_summary_rows)


def load_cluster_background_counts_no_commit(conn, cluster_background_count_rows):
    if not cluster_background_count_rows:
        return 0

    sql = """
    INSERT INTO cluster_background_counts (
        organism,
        cluster_number,
        centroid,
        segment,
        background_type,
        background_value,
        n_sequences
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :segment,
        :background_type,
        :background_value,
        :n_sequences
    )
    ON CONFLICT(
        organism,
        cluster_number,
        segment,
        background_type,
        background_value
    ) DO UPDATE SET
        centroid = excluded.centroid,
        n_sequences = excluded.n_sequences;
    """

    conn.executemany(sql, cluster_background_count_rows)
    return len(cluster_background_count_rows)


def load_cluster_descriptors_no_commit(conn, cluster_descriptor_rows):
    if not cluster_descriptor_rows:
        return 0

    sql = """
    INSERT INTO cluster_descriptors (
        organism,
        cluster_number,
        representative,
        n_sequences,
        genotypes_json,
        segments_json,
        hosts_json,
        countries_json,
        dominant_genotype,
        dominant_segment,
        dominant_host,
        dominant_country,
        n_genotypes,
        n_segments,
        n_hosts,
        n_countries,
        is_singleton,
        is_mixed_genotype,
        is_mixed_segment,
        is_multi_host,
        is_multi_country,
        length_min,
        length_mean,
        length_median,
        length_max,
        length_std
    )
    VALUES (
        :organism,
        :cluster_number,
        :representative,
        :n_sequences,
        :genotypes_json,
        :segments_json,
        :hosts_json,
        :countries_json,
        :dominant_genotype,
        :dominant_segment,
        :dominant_host,
        :dominant_country,
        :n_genotypes,
        :n_segments,
        :n_hosts,
        :n_countries,
        :is_singleton,
        :is_mixed_genotype,
        :is_mixed_segment,
        :is_multi_host,
        :is_multi_country,
        :length_min,
        :length_mean,
        :length_median,
        :length_max,
        :length_std
    )
    ON CONFLICT(organism, cluster_number) DO UPDATE SET
        representative = excluded.representative,
        n_sequences = excluded.n_sequences,
        genotypes_json = excluded.genotypes_json,
        segments_json = excluded.segments_json,
        hosts_json = excluded.hosts_json,
        countries_json = excluded.countries_json,
        dominant_genotype = excluded.dominant_genotype,
        dominant_segment = excluded.dominant_segment,
        dominant_host = excluded.dominant_host,
        dominant_country = excluded.dominant_country,
        n_genotypes = excluded.n_genotypes,
        n_segments = excluded.n_segments,
        n_hosts = excluded.n_hosts,
        n_countries = excluded.n_countries,
        is_singleton = excluded.is_singleton,
        is_mixed_genotype = excluded.is_mixed_genotype,
        is_mixed_segment = excluded.is_mixed_segment,
        is_multi_host = excluded.is_multi_host,
        is_multi_country = excluded.is_multi_country,
        length_min = excluded.length_min,
        length_mean = excluded.length_mean,
        length_median = excluded.length_median,
        length_max = excluded.length_max,
        length_std = excluded.length_std;
    """

    conn.executemany(sql, cluster_descriptor_rows)
    return len(cluster_descriptor_rows)


def load_cluster_filter_counts_no_commit(conn, cluster_filter_count_rows):
    if not cluster_filter_count_rows:
        return 0

    sql = """
    INSERT INTO cluster_filter_counts (
        organism,
        cluster_number,
        centroid,
        segment,
        genotype,
        ha_subtype,
        na_subtype,
        host,
        country,
        n_sequences
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :segment,
        :genotype,
        :ha_subtype,
        :na_subtype,
        :host,
        :country,
        :n_sequences
    )
    ON CONFLICT(
        organism,
        cluster_number,
        segment,
        genotype,
        ha_subtype,
        na_subtype,
        host,
        country
    ) DO UPDATE SET
        centroid = excluded.centroid,
        n_sequences = excluded.n_sequences;
    """

    conn.executemany(sql, cluster_filter_count_rows)
    return len(cluster_filter_count_rows)


def load_cluster_assembly_no_commit(conn, cluster_assembly_rows):
    """
    Load rows into cluster_assembly without committing.
    """
    if not cluster_assembly_rows:
        return 0

    sql = """
    INSERT INTO cluster_assembly (
        organism,
        cluster_number,
        centroid,
        assembly,
        segment,
        n_sequences,
        is_centroid_assembly
    )
    VALUES (
        :organism,
        :cluster_number,
        :centroid,
        :assembly,
        :segment,
        :n_sequences,
        :is_centroid_assembly
    )
    ON CONFLICT(
        organism,
        cluster_number,
        assembly,
        segment
    ) DO UPDATE SET
        centroid = excluded.centroid,
        n_sequences = excluded.n_sequences,
        is_centroid_assembly = excluded.is_centroid_assembly;
    """

    conn.executemany(sql, cluster_assembly_rows)
    return len(cluster_assembly_rows)


def load_cluster_assembly_links_no_commit(conn, cluster_assembly_link_rows):
    """
    Load rows into cluster_assembly_links without committing.
    """
    if not cluster_assembly_link_rows:
        return 0

    sql = """
    INSERT INTO cluster_assembly_links (
        organism,
        cluster_number,
        linked_organism,
        linked_cluster_number,
        linked_cluster_segment,
        n_connections
    )
    VALUES (
        :organism,
        :cluster_number,
        :linked_organism,
        :linked_cluster_number,
        :linked_cluster_segment,
        :n_connections
    )
    ON CONFLICT(
        organism,
        cluster_number,
        linked_organism,
        linked_cluster_number
    ) DO UPDATE SET
        linked_cluster_segment = excluded.linked_cluster_segment,
        n_connections = excluded.n_connections;
    """

    conn.executemany(sql, cluster_assembly_link_rows)
    return len(cluster_assembly_link_rows)


def load_all_tables_transactional(
    db_path,
    metadata_rows,
    sequence_rows,
    cluster_rows,
    cluster_composition_rows,
    cluster_summary_rows=None,
    cluster_background_count_rows=None,
    cluster_descriptor_rows=None,
    cluster_filter_count_rows=None,
    cluster_assembly_rows=None,
    cluster_assembly_link_rows=None
):
    """
    Load all ETL-ready rows into SQLite inside a single transaction.

    If any table fails, the whole import is rolled back.
    """
    cluster_summary_rows = cluster_summary_rows or []
    cluster_background_count_rows = cluster_background_count_rows or []
    cluster_descriptor_rows = cluster_descriptor_rows or []
    cluster_filter_count_rows = cluster_filter_count_rows or []
    cluster_assembly_rows = cluster_assembly_rows or []
    cluster_assembly_link_rows = cluster_assembly_link_rows or []

    conn = connect_sqlite(db_path)

    try:
        conn.execute("BEGIN;")

        n_metadata = load_metadata_no_commit(conn, metadata_rows)
        n_sequences = load_sequences_no_commit(conn, sequence_rows)
        n_clusters = load_clusters_no_commit(conn, cluster_rows)
        n_cluster_composition = load_cluster_composition_no_commit(
            conn,
            cluster_composition_rows
        )

        n_cluster_summary = load_cluster_summary_no_commit(
            conn,
            cluster_summary_rows
        )

        n_cluster_background_counts = load_cluster_background_counts_no_commit(
            conn,
            cluster_background_count_rows
        )

        n_cluster_descriptors = load_cluster_descriptors_no_commit(
            conn,
            cluster_descriptor_rows
        )

        n_cluster_filter_counts = load_cluster_filter_counts_no_commit(
            conn,
            cluster_filter_count_rows
        )

        n_cluster_assembly = load_cluster_assembly_no_commit(
            conn,
            cluster_assembly_rows
        )

        n_cluster_assembly_links = load_cluster_assembly_links_no_commit(
            conn,
            cluster_assembly_link_rows
        )

        conn.commit()

        return {
            "metadata": n_metadata,
            "sequences": n_sequences,
            "clusters": n_clusters,
            "cluster_composition": n_cluster_composition,
            "cluster_summary": n_cluster_summary,
            "cluster_background_counts": n_cluster_background_counts,
            "cluster_descriptors": n_cluster_descriptors,
            "cluster_filter_counts": n_cluster_filter_counts,
            "cluster_assembly": n_cluster_assembly,
            "cluster_assembly_links": n_cluster_assembly_links,
        }

    except Exception:
        conn.rollback()
        raise

    finally:
        conn.close()

def parse_args():
    """
    Parse command-line arguments for the ETL pipeline.

    Required inputs
    ---------------
    --metadata
        Path to the metadata CSV/semicolon-separated file.

    --fasta
        Path to the FASTA file containing sequences.

    --clstr
        Path to the CD-HIT .clstr file.

    --schema
        Path to the SQLite schema.sql file.

    --db
        Path to the SQLite database file to create/use.

    --load-mode
        Loading mode:
        - basic: commits table by table
        - transactional: all-or-nothing transaction.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Import metadata, sequence and cluster data into a SQLite database."
        )
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to the metadata CSV file."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to the FASTA file containing sequences."
    )

    parser.add_argument(
        "--clstr",
        required=True,
        help="Path to the CD-HIT .clstr file."
    )

    parser.add_argument(
        "--schema",
        required=True,
        help="Path to the SQLite schema.sql file."
    )

    parser.add_argument(
        "--db",
        required=True,
        help="Path to the SQLite database file."
    )

    parser.add_argument(
        "--load-mode",
        choices=["basic", "transactional"],
        default="transactional",
        help=(
            "SQLite loading mode. "
            "'basic' commits table by table; "
            "'transactional' loads all tables in one transaction. "
            "Default: transactional."
        )
    )

    parser.add_argument(
        "--recreate",
        action="store_true",
        help=(
            "Delete the existing SQLite database before applying the schema "
            "and loading data."
        )
    )

    parser.add_argument(
        "--metadata-delimiter",
        default=";",
        help="Delimiter used in the metadata file. Default: ';'."
    )

    parser.add_argument(
        "--record-completeness",
        choices=["fail", "warn"],
        default="warn",
        help=(
            "Behaviour of the per-species record completeness check: every "
            "sequence must have metadata, a cluster composition entry and a "
            "FASTA sequence. 'fail' aborts the load; 'warn' prints the report "
            "and continues. Default: warn."
        )
    )

    parser.add_argument(
        "--fetch-missing",
        action="store_true",
        help=(
            "Download sequences missing from the FASTA file from NCBI using "
            "Bio.Entrez before loading. The download is resumable: see "
            "--fetch-output and --fetch-progress."
        )
    )

    parser.add_argument(
        "--entrez-email",
        default=None,
        help=(
            "Email address required by NCBI Entrez. Required when "
            "--fetch-missing is used."
        )
    )

    parser.add_argument(
        "--entrez-api-key",
        default=None,
        help="NCBI Entrez API key (optional, raises rate limits)."
    )

    parser.add_argument(
        "--fetch-output",
        default="data/fetched_sequences.fasta",
        help=(
            "FASTA file that fetched sequences are appended to. "
            "Default: data/fetched_sequences.fasta."
        )
    )

    parser.add_argument(
        "--fetch-progress",
        default="data/fetch_progress.json",
        help=(
            "JSON progress file used to resume an interrupted NCBI download. "
            "Default: data/fetch_progress.json."
        )
    )

    parser.add_argument(
        "--fetch-batch-size",
        type=int,
        default=200,
        help="Accessions per NCBI efetch request. Default: 200."
    )

    parser.add_argument(
        "--fetch-sleep",
        type=float,
        default=None,
        help=(
            "Seconds to sleep between NCBI requests. Defaults to 0.34 without "
            "an API key, 0.1 with one."
        )
    )

    return parser.parse_args()

def main():
    '''
    Main function of ETL pipeline, accepts 3 files, 
    extracts and transforms the data, then loads to a sqlite table
    '''
    args=parse_args()
    print("ETL configuration")
    print("-----------------")
    print(f"Metadata file: {args.metadata}")
    print(f"Metadata delimiter: {args.metadata_delimiter!r}")
    print(f"FASTA file: {args.fasta}")
    print(f"CLSTR file: {args.clstr}")
    print(f"Schema file: {args.schema}")
    print(f"SQLite database: {args.db}")
    print(f"Load mode: {args.load_mode}")
    print(f"Recreate database: {args.recreate}")

    metadata_rows=parse_metadata_csv(metadata_csv_path=args.metadata,
                                     delimiter=args.metadata_delimiter)
    organism=infer_single_organism_from_metadata(metadata_rows)
    fasta_dict=seq_get(args.fasta)
    acc_header=decompose_fasta_headers(fasta_dict)
    acc_seq=accession_to_seq_dict(fasta_dict)
    fasta_dict=join_fasta_dicts(acc_header,acc_seq)
    sequence_rows=transform_seqdict(fasta_dict)
    clusters=mine_clstr_table(args.clstr,organism)
    cluster_rows=transform_clustdict(clusters)
    comp_dict=mine_clstr_elements(args.clstr,organism)
    comp_rows=transform_compdict(comp_dict)

    sequence_rows = deduplicate_rows(sequence_rows, label="sequence rows")
    comp_rows = deduplicate_rows(comp_rows, label="cluster composition rows")

    derived_cluster_rows = build_derived_cluster_rows(
        metadata_rows=metadata_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=comp_rows
    )
    cluster_summary_rows = derived_cluster_rows["cluster_summary"]
    cluster_background_count_rows = derived_cluster_rows["cluster_background_counts"]
    cluster_descriptor_rows = derived_cluster_rows["cluster_descriptors"]
    cluster_filter_count_rows = derived_cluster_rows["cluster_filter_counts"]
    cluster_assembly_rows = derived_cluster_rows["cluster_assembly"]
    cluster_assembly_link_rows = derived_cluster_rows["cluster_assembly_links"]

    print("Transformation summary")
    print("----------------------")
    print(f"Metadata rows: {len(metadata_rows)}")
    print(f"Sequence rows: {len(sequence_rows)}")
    print(f"Cluster rows: {len(cluster_rows)}")
    print(f"Cluster composition rows: {len(comp_rows)}")
    print(f"Inferred organism: {organism}")
    print(f"Cluster summary rows: {len(cluster_summary_rows)}")
    print(f"Cluster background count rows: {len(cluster_background_count_rows)}")
    print(f"Cluster descriptor rows: {len(cluster_descriptor_rows)}")
    print(f"Cluster filter count rows: {len(cluster_filter_count_rows)}")
    print(f"Cluster assembly rows: {len(cluster_assembly_rows)}")
    print(f"Cluster assembly link rows: {len(cluster_assembly_link_rows)}")

    validate_cross_references(
        metadata_rows=metadata_rows,
        sequence_rows=sequence_rows,
        cluster_rows=cluster_rows,
        cluster_composition_rows=comp_rows
    )

    # -- Per-species record completeness + optional NCBI backfill ---------

    initial_gaps = compute_record_gaps(
        metadata_rows=metadata_rows,
        sequence_rows=sequence_rows,
        cluster_composition_rows=comp_rows,
        organism=organism,
    )

    report_record_gaps(initial_gaps)

    missing_fasta = missing_fasta_accessions(initial_gaps)

    if missing_fasta and args.fetch_missing:
        if not args.entrez_email:
            raise ValueError(
                "--fetch-missing requires --entrez-email (NCBI usage policy)."
            )

        print("Fetching missing sequences from NCBI")
        print("------------------------------------")

        fetched_rows = ncbi_fetch.fetch_and_build_sequence_rows(
            accessions=missing_fasta,
            output_fasta=args.fetch_output,
            progress_file=args.fetch_progress,
            email=args.entrez_email,
            api_key=args.entrez_api_key,
            batch_size=args.fetch_batch_size,
            sleep_seconds=args.fetch_sleep,
        )

        if fetched_rows:
            print(f"Fetched {len(fetched_rows)} sequences from NCBI.")
            sequence_rows = deduplicate_rows(
                sequence_rows + fetched_rows,
                label="sequence rows after fetch",
            )

            remaining_gaps = compute_record_gaps(
                metadata_rows=metadata_rows,
                sequence_rows=sequence_rows,
                cluster_composition_rows=comp_rows,
                organism=organism,
            )

            report_record_gaps(remaining_gaps)

            if missing_fasta_accessions(remaining_gaps):
                print(
                    "NOTE: some sequences are still missing a FASTA sequence. "
                    "Re-run db/ncbi_fetch.py (or this command) to resume the "
                    "download where it left off."
                )

    if args.record_completeness == "fail":
        validate_species_record_completeness(
            metadata_rows=metadata_rows,
            sequence_rows=sequence_rows,
            cluster_composition_rows=comp_rows,
            organism=organism,
            mode="fail",
        )
    else:
        report_assembly_gaps(
            compute_assembly_gaps(
                metadata_rows=metadata_rows,
                cluster_composition_rows=comp_rows,
            )
        )

    if missing_fasta and not args.fetch_missing:
        print(
            "NOTE: some sequences are missing from the FASTA file. "
            "Run db/ncbi_fetch.py or re-run with --fetch-missing to download "
            "them from NCBI."
        )

    initialise_database(
        db_path=args.db,
        schema_path=args.schema,
        recreate=args.recreate
    )

    if args.load_mode == "basic":
        load_summary = load_all_tables(
            db_path=args.db,
            metadata_rows=metadata_rows,
            sequence_rows=sequence_rows,
            cluster_rows=cluster_rows,
            cluster_composition_rows=comp_rows,
            cluster_summary_rows=cluster_summary_rows,
            cluster_background_count_rows=cluster_background_count_rows,
            cluster_descriptor_rows=cluster_descriptor_rows,
            cluster_filter_count_rows=cluster_filter_count_rows,
            cluster_assembly_rows=cluster_assembly_rows,
            cluster_assembly_link_rows=cluster_assembly_link_rows
        )

    elif args.load_mode == "transactional":
        load_summary = load_all_tables_transactional(
            db_path=args.db,
            metadata_rows=metadata_rows,
            sequence_rows=sequence_rows,
            cluster_rows=cluster_rows,
            cluster_composition_rows=comp_rows,
            cluster_summary_rows=cluster_summary_rows,
            cluster_background_count_rows=cluster_background_count_rows,
            cluster_descriptor_rows=cluster_descriptor_rows,
            cluster_filter_count_rows=cluster_filter_count_rows,
            cluster_assembly_rows=cluster_assembly_rows,
            cluster_assembly_link_rows=cluster_assembly_link_rows
        )

    else:
        raise ValueError(f"Unknown load mode: {args.load_mode}")

    print("Load summary")
    print("------------")
    for table_name, n_rows in load_summary.items():
        print(f"{table_name}: {n_rows}")

if __name__ == "__main__":
    main()