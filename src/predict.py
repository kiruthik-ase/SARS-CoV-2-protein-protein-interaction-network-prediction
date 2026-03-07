"""
Prediction script for SARS-CoV-2 × Human PPI prediction.
Given a viral UniProt ID and a human UniProt ID (or raw sequences),
predict the probability of interaction using the trained XGBoost model.
"""

import pandas as pd
import numpy as np
import os
import re
import sys
import json
import joblib
import argparse
import requests
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from xgboost import XGBClassifier

# -------------------------
# Paths
# -------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(BASE_DIR, "models", "ppi_xgboost_model.json")
FEATURE_COLS_PATH = os.path.join(BASE_DIR, "models", "feature_columns.pkl")

# -------------------------
# Feature extraction (same as pipeline.py)
# -------------------------
AA = "ACDEFGHIKLMNPQRSTVWY"

def clean_seq(s):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(s).upper())

def aac(seq):
    c = Counter(seq)
    return [c[a] / max(1, len(seq)) for a in AA]

def dpc(seq):
    pairs = [seq[i:i+2] for i in range(len(seq)-1)]
    c = Counter(pairs)
    total = max(1, len(pairs))
    return [c[a+b] / total for a in AA for b in AA]

def physchem(seq):
    if len(seq) < 5:
        return [len(seq), 0.0, 0.0, 0.0, 0.0]
    try:
        pa = ProteinAnalysis(seq)
        return [len(seq), pa.molecular_weight(), pa.gravy(), pa.aromaticity(), pa.instability_index()]
    except Exception:
        return [len(seq), 0.0, 0.0, 0.0, 0.0]

def extract_pair_features(viral_seq, human_seq):
    """Extract feature vector for a viral-human protein pair."""
    v = clean_seq(viral_seq)
    h = clean_seq(human_seq)
    if len(v) < 5 or len(h) < 5:
        raise ValueError(f"Sequences too short after cleaning (viral: {len(v)}, human: {len(h)})")
    feat = aac(v) + dpc(v) + physchem(v) + aac(h) + dpc(h) + physchem(h)
    return np.array(feat)

# -------------------------
# Sequence fetching
# -------------------------
def fetch_sequence(uniprot_id):
    """Fetch protein sequence from UniProt API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200 and r.text.startswith(">"):
            return "".join(r.text.splitlines()[1:])
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
    return None

# -------------------------
# Load model
# -------------------------
def load_model():
    model = XGBClassifier()
    model.load_model(MODEL_PATH)
    feature_cols = joblib.load(FEATURE_COLS_PATH)
    return model, feature_cols

# -------------------------
# Predict
# -------------------------
def predict_interaction(viral_id_or_seq, human_id_or_seq, model=None, feature_cols=None):
    """
    Predict interaction probability between a viral and human protein.
    
    Args:
        viral_id_or_seq: UniProt ID (e.g., "P0DTC2") or raw amino acid sequence
        human_id_or_seq: UniProt ID (e.g., "Q9BYF1") or raw amino acid sequence
    
    Returns:
        dict with probability and prediction
    """
    if model is None:
        model, feature_cols = load_model()

    # Determine if input is UniProt ID or raw sequence
    if len(viral_id_or_seq) < 30 and not any(c in viral_id_or_seq for c in "acdefghiklmnpqrstvwy"):
        print(f"Fetching viral sequence for {viral_id_or_seq}...")
        viral_seq = fetch_sequence(viral_id_or_seq)
        if viral_seq is None:
            raise ValueError(f"Could not fetch sequence for {viral_id_or_seq}")
        viral_label = viral_id_or_seq
    else:
        viral_seq = viral_id_or_seq
        viral_label = "custom_viral"

    if len(human_id_or_seq) < 30 and not any(c in human_id_or_seq for c in "acdefghiklmnpqrstvwy"):
        print(f"Fetching human sequence for {human_id_or_seq}...")
        human_seq = fetch_sequence(human_id_or_seq)
        if human_seq is None:
            raise ValueError(f"Could not fetch sequence for {human_id_or_seq}")
        human_label = human_id_or_seq
    else:
        human_seq = human_id_or_seq
        human_label = "custom_human"

    # Extract features
    feat_vector = extract_pair_features(viral_seq, human_seq)
    X = pd.DataFrame([feat_vector], columns=feature_cols)

    # Predict
    prob = model.predict_proba(X)[0, 1]
    pred = int(prob >= 0.5)

    return {
        "viral_protein": viral_label,
        "human_protein": human_label,
        "interaction_probability": round(float(prob), 4),
        "prediction": "INTERACTING" if pred else "NON-INTERACTING",
        "viral_seq_length": len(clean_seq(viral_seq)),
        "human_seq_length": len(clean_seq(human_seq))
    }

# -------------------------
# Batch predict
# -------------------------
def batch_predict(pairs_csv, output_csv=None):
    """
    Batch predict from a CSV with columns: viral_uniprot, human_uniprot
    """
    model, feature_cols = load_model()
    pairs = pd.read_csv(pairs_csv)
    
    results = []
    for _, row in pairs.iterrows():
        try:
            result = predict_interaction(
                row["viral_uniprot"], row["human_uniprot"],
                model=model, feature_cols=feature_cols
            )
            results.append(result)
        except Exception as e:
            results.append({
                "viral_protein": row["viral_uniprot"],
                "human_protein": row["human_uniprot"],
                "interaction_probability": None,
                "prediction": f"ERROR: {e}",
                "viral_seq_length": 0,
                "human_seq_length": 0
            })
    
    results_df = pd.DataFrame(results)
    if output_csv:
        results_df.to_csv(output_csv, index=False)
        print(f"Saved predictions → {output_csv}")
    return results_df

# -------------------------
# CLI
# -------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict SARS-CoV-2 × Human PPI")
    parser.add_argument("--viral", required=True, help="Viral UniProt ID or sequence")
    parser.add_argument("--human", required=True, help="Human UniProt ID or sequence")
    parser.add_argument("--batch", help="CSV file with viral_uniprot, human_uniprot columns")
    parser.add_argument("--output", help="Output CSV for batch predictions")
    args = parser.parse_args()

    if args.batch:
        results = batch_predict(args.batch, args.output)
        print(results)
    else:
        result = predict_interaction(args.viral, args.human)
        print("\n=== PREDICTION RESULT ===")
        for k, v in result.items():
            print(f"  {k}: {v}")
