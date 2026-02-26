import pandas as pd
import numpy as np
from collections import Counter

print("\n=== STEP B2: Feature extraction for negatives ===")

# -----------------------------
# SAME FEATURE FUNCTIONS AS POSITIVES
# -----------------------------
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def aac(seq):
    seq = seq.upper()
    length = len(seq)
    counts = Counter(seq)
    return {f"aac_{aa}": counts.get(aa, 0) / length for aa in AMINO_ACIDS}

def dpc(seq):
    seq = seq.upper()
    pairs = [seq[i:i+2] for i in range(len(seq)-1)]
    counts = Counter(pairs)
    total = sum(counts.values())
    return {f"dpc_{a}{b}": counts.get(a+b, 0) / total for a in AMINO_ACIDS for b in AMINO_ACIDS}

def physchem(seq):
    length = len(seq)
    mw = length * 110  # simple approximation (same as before)
    gravy = sum([1 if aa in "AILMFWV" else -1 for aa in seq]) / length
    arom = sum([1 for aa in seq if aa in "FWY"]) / length
    instab = length * 0.7
    return {
        "len": length,
        "mw": mw,
        "gravy": gravy,
        "arom": arom,
        "instab": instab
    }

def extract_features(vseq, hseq):
    feat = {}

    # Viral features
    feat.update({f"v_{k}": v for k, v in aac(vseq).items()})
    feat.update({f"v_{k}": v for k, v in dpc(vseq).items()})
    feat.update({f"v_{k}": v for k, v in physchem(vseq).items()})

    # Human features
    feat.update({f"h_{k}": v for k, v in aac(hseq).items()})
    feat.update({f"h_{k}": v for k, v in dpc(hseq).items()})
    feat.update({f"h_{k}": v for k, v in physchem(hseq).items()})

    return feat

# -----------------------------
# LOAD NEGATIVES
# -----------------------------
neg = pd.read_csv("negatives_with_sequences.csv")

features = []
for i, row in neg.iterrows():
    f = extract_features(row["viral_sequence"], row["human_sequence"])
    f["viral_uniprot"] = row["viral_uniprot"]
    f["human_uniprot"] = row["human_uniprot"]
    f["label"] = 0
    features.append(f)

neg_feat = pd.DataFrame(features)

neg_feat.to_csv("ml_negative_features.csv", index=False)

print("Saved → ml_negative_features.csv")
print("Shape:", neg_feat.shape)
print(neg_feat.head())
