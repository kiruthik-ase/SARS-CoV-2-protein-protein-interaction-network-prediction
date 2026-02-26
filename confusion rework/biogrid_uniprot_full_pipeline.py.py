import pandas as pd
import numpy as np
import requests
import time
import itertools
from collections import Counter

############################################
# FILES
############################################

BIOGRID_FILE = r"C:\Users\LENOVO\bio sem 4\Raw Data\BIOGRID-PROJECT-covid19_coronavirus_project-5.0.251\BIOGRID-PROJECT-covid19_coronavirus_project-INTERACTIONS-5.0.251.tab3.txt"
HUMAN_FASTA_FILE = r"C:\Users\LENOVO\bio sem 4\csvs\human_sequences_clean.csv"
OUT_DATASET = "final_ppi_dataset.csv"

MIN_SEQ_LEN = 30
NEG_RATIO = 1

############################################
# AMINO ACIDS
############################################

AA = "ACDEFGHIKLMNPQRSTVWY"

############################################
# FEATURE FUNCTIONS (REAL, NOT FAKE)
############################################

def aac(seq):
    c = Counter(seq)
    return [c[a] / len(seq) for a in AA]

def dpc(seq):
    pairs = [seq[i:i+2] for i in range(len(seq)-1)]
    c = Counter(pairs)
    total = len(pairs)
    return [c[a+b] / total for a in AA for b in AA]

def physchem(seq):
    mw = len(seq) * 110
    hydro = sum(seq.count(a) for a in "AILMFWV") / len(seq)
    arom = sum(seq.count(a) for a in "FWY") / len(seq)
    return [len(seq), mw, hydro, arom]

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
    v = viral_seq[r["viral_uniprot"]]
    h = human_seq[r["human_uniprot"]]

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
    ["v_len","v_mw","v_hydro","v_arom"] +
    [f"h_aac_{a}" for a in AA] +
    [f"h_dpc_{a}{b}" for a in AA for b in AA] +
    ["h_len","h_mw","h_hydro","h_arom"]
)

final_df = pd.DataFrame(rows, columns=columns)
final_df.to_csv(OUT_DATASET, index=False)

print("\n=== DONE ===")
print("Final samples:", len(final_df))
print("Saved →", OUT_DATASET)
