import pandas as pd
import numpy as np
import requests
import time
import os
import re
import itertools
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis

############################################
# FILES
############################################

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIOGRID_FILE = os.path.join(BASE_DIR, "data", "raw", "BIOGRID-PROJECT-covid19_coronavirus_project-5.0.251", "BIOGRID-PROJECT-covid19_coronavirus_project-INTERACTIONS-5.0.251.tab3.txt")
HUMAN_FASTA_FILE = os.path.join(BASE_DIR, "data", "processed", "human_sequences_clean.csv")
OUT_DATASET = os.path.join(BASE_DIR, "data", "processed", "final_ppi_dataset.csv")

MIN_SEQ_LEN = 30
NEG_RATIO = 1

############################################
# AMINO ACIDS
############################################

AA = "ACDEFGHIKLMNPQRSTVWY"

############################################
# FEATURE FUNCTIONS
############################################

def clean_seq(s):
    """Remove non-standard amino acid characters."""
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(s).upper())

def aac(seq):
    c = Counter(seq)
    return [c[a] / len(seq) for a in AA]

def dpc(seq):
    pairs = [seq[i:i+2] for i in range(len(seq)-1)]
    c = Counter(pairs)
    total = len(pairs)
    return [c[a+b] / total for a in AA for b in AA]

def physchem(seq):
    """Compute physicochemical properties using Biopython ProteinAnalysis."""
    if len(seq) < 5:
        return [len(seq), 0.0, 0.0, 0.0, 0.0]
    try:
        pa = ProteinAnalysis(seq)
        return [len(seq), pa.molecular_weight(), pa.gravy(), pa.aromaticity(), pa.instability_index()]
    except Exception:
        return [len(seq), 0.0, 0.0, 0.0, 0.0]

############################################
# SAFE VIRAL FASTA FETCH (ONLY ~19)
############################################

def fetch_fasta(uniprot):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200 and r.text.startswith(">"):
            return "".join(r.text.splitlines()[1:])
    except:
        pass
    return None

############################################
# LOAD BIOGRID
############################################

print("\n=== LOADING BIOGRID ===")

bg = pd.read_csv(BIOGRID_FILE, sep="\t", low_memory=False)

bg = bg[
    (bg["Organism Name Interactor A"] == "Severe acute respiratory syndrome coronavirus 2") &
    (bg["Organism Name Interactor B"] == "Homo sapiens")
]

bg = bg[[
    "SWISS-PROT Accessions Interactor A",
    "SWISS-PROT Accessions Interactor B"
]].dropna()

bg.columns = ["viral_uniprot", "human_uniprot"]
bg = bg[~bg["viral_uniprot"].str.contains("\\|")]
bg = bg[~bg["human_uniprot"].str.contains("\\|")]

print("Positive interactions:", len(bg))
print("Unique viral proteins:", bg["viral_uniprot"].nunique())
print("Unique human proteins:", bg["human_uniprot"].nunique())

############################################
# LOAD HUMAN SEQUENCES (LOCAL – FAST)
############################################

print("\n=== LOADING HUMAN FASTA (LOCAL) ===")

human_df = pd.read_csv(HUMAN_FASTA_FILE)
human_seq = dict(zip(human_df["uniprot"], human_df["sequence"]))

############################################
# FETCH VIRAL SEQUENCES (ONLY 19)
############################################

print("\n=== FETCHING VIRAL FASTA ===")

viral_seq = {}
for v in bg["viral_uniprot"].unique():
    s = fetch_fasta(v)
    if s and len(s) >= MIN_SEQ_LEN:
        viral_seq[v] = s

############################################
# FILTER VALID PAIRS
############################################

bg = bg[
    bg["viral_uniprot"].isin(viral_seq) &
    bg["human_uniprot"].isin(human_seq)
]

print("Usable positives:", len(bg))
print("Final viral proteins:", bg["viral_uniprot"].nunique())

############################################
# GENERATE NEGATIVES
############################################

all_pairs = set(itertools.product(viral_seq.keys(), human_seq.keys()))
positive_pairs = set(zip(bg["viral_uniprot"], bg["human_uniprot"]))

negatives = list(all_pairs - positive_pairs)
np.random.shuffle(negatives)
negatives = negatives[:len(bg) * NEG_RATIO]

pos_df = bg.copy()
pos_df["label"] = 1

neg_df = pd.DataFrame(negatives, columns=["viral_uniprot", "human_uniprot"])
neg_df["label"] = 0

data = pd.concat([pos_df, neg_df], ignore_index=True)

############################################
# FEATURE EXTRACTION
############################################

print("\n=== FEATURE EXTRACTION ===")

rows = []

for _, r in data.iterrows():
    v = clean_seq(viral_seq[r["viral_uniprot"]])
    h = clean_seq(human_seq[r["human_uniprot"]])

    if len(v) < 5 or len(h) < 5:
        continue

    feat = (
        aac(v) + dpc(v) + physchem(v) +
        aac(h) + dpc(h) + physchem(h)
    )

    rows.append([r["viral_uniprot"], r["human_uniprot"], r["label"]] + feat)

############################################
# SAVE FINAL DATASET
############################################

columns = (
    ["viral_uniprot", "human_uniprot", "label"] +
    [f"v_aac_{a}" for a in AA] +
    [f"v_dpc_{a}{b}" for a in AA for b in AA] +
    ["v_len","v_mw","v_gravy","v_arom","v_instab"] +
    [f"h_aac_{a}" for a in AA] +
    [f"h_dpc_{a}{b}" for a in AA for b in AA] +
    ["h_len","h_mw","h_gravy","h_arom","h_instab"]
)

final_df = pd.DataFrame(rows, columns=columns)
final_df.to_csv(OUT_DATASET, index=False)

print("\n=== DONE ===")
print("Final samples:", len(final_df))
print("Saved →", OUT_DATASET)
