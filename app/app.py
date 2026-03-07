"""
Flask web app for SARS-CoV-2 PPI prediction.
Provides interactive prediction, network visualization, and model inferences.
"""

from flask import Flask, render_template, request, jsonify
import pandas as pd
import numpy as np
import os
import re
import sys
import json
import joblib
import requests
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from xgboost import XGBClassifier

app = Flask(__name__)

# -------------------------
# Paths
# -------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(BASE_DIR, "models", "ppi_xgboost_model.json")
FEATURE_COLS_PATH = os.path.join(BASE_DIR, "models", "feature_columns.pkl")
DATASET_PATH = os.path.join(BASE_DIR, "data", "processed", "final_ppi_dataset.csv")
RESULTS_PATH = os.path.join(BASE_DIR, "results", "viral_wise_results.csv")
IMPORTANCE_PATH = os.path.join(BASE_DIR, "results", "feature_importance_full.csv")

# -------------------------
# Feature extraction (same as pipeline)
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
        return [len(seq), 0.0, 0.0, 0.0]
    try:
        pa = ProteinAnalysis(seq)
        return [len(seq), pa.molecular_weight(), pa.gravy(), pa.aromaticity()]
    except Exception:
        return [len(seq), 0.0, 0.0, 0.0]

def fetch_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200 and r.text.startswith(">"):
            seq = "".join(r.text.splitlines()[1:])
            # Also get protein name from header
            header = r.text.splitlines()[0]
            name = header.split("|")[-1].split(" OS=")[0].strip() if "|" in header else uniprot_id
            return seq, name
    except:
        pass
    return None, None

# -------------------------
# Load model + data on startup
# -------------------------
print("Loading model...")
model = XGBClassifier()
model.load_model(MODEL_PATH)
feature_cols = joblib.load(FEATURE_COLS_PATH)

# Load network data for visualization
print("Loading network data...")
network_data = {"nodes": [], "edges": []}
viral_results = {}

try:
    df = pd.read_csv(DATASET_PATH, usecols=["viral_uniprot", "human_uniprot", "label"])
    positives = df[df["label"] == 1]

    # Build network: viral proteins and their top human targets
    viral_proteins = positives["viral_uniprot"].unique().tolist()
    
    # Viral protein name mapping
    VIRAL_NAMES = {
        "P0DTC1": "ORF1ab (pp1a)", "P0DTC2": "Spike (S)", "P0DTC3": "ORF3a",
        "P0DTC4": "Envelope (E)", "P0DTC5": "Membrane (M)", "P0DTC6": "ORF6",
        "P0DTC7": "ORF7a", "P0DTC8": "ORF8", "P0DTC9": "Nucleocapsid (N)",
        "P0DTD1": "ORF1ab (replicase)", "P0DTD2": "ORF9b", "P0DTD3": "ORF14",
        "P0DTD8": "ORF7b", "A0A663DJA2": "ORF10"
    }
    
    nodes = []
    edges = []
    
    for vp in viral_proteins:
        nodes.append({"id": vp, "label": VIRAL_NAMES.get(vp, vp), "type": "viral"})
        
        # Get top human targets for this viral protein
        targets = positives[positives["viral_uniprot"] == vp]["human_uniprot"].value_counts().head(8).index.tolist()
        for hp in targets:
            if not any(n["id"] == hp for n in nodes):
                nodes.append({"id": hp, "label": hp, "type": "human"})
            edges.append({"source": vp, "target": hp})
    
    network_data = {"nodes": nodes, "edges": edges}
    
    # Compute inferences
    interaction_counts = positives.groupby("viral_uniprot")["human_uniprot"].nunique().to_dict()
    total_interactions = len(positives)
    total_human = positives["human_uniprot"].nunique()
    
except Exception as e:
    print(f"Warning: Could not load dataset: {e}")
    interaction_counts = {}
    total_interactions = 0
    total_human = 0

# Load results
try:
    results_df = pd.read_csv(RESULTS_PATH)
    viral_results = results_df.to_dict('records')
except:
    viral_results = []

# Load feature importance
top_features = []
try:
    imp_df = pd.read_csv(IMPORTANCE_PATH).head(15)
    top_features = imp_df.to_dict('records')
except:
    pass

print("App ready!")

# -------------------------
# Routes
# -------------------------
@app.route("/")
def index():
    return render_template("index.html",
                         network_data=json.dumps(network_data),
                         viral_results=viral_results,
                         top_features=top_features,
                         interaction_counts=interaction_counts,
                         total_interactions=total_interactions,
                         total_human=total_human,
                         viral_names=json.dumps(VIRAL_NAMES))

@app.route("/predict", methods=["POST"])
def predict():
    data = request.json
    viral_id = data.get("viral_id", "").strip()
    human_id = data.get("human_id", "").strip()
    
    if not viral_id or not human_id:
        return jsonify({"error": "Both viral and human protein IDs are required"}), 400
    
    try:
        # Fetch sequences
        viral_seq, viral_name = fetch_sequence(viral_id)
        if not viral_seq:
            return jsonify({"error": f"Could not fetch sequence for viral protein: {viral_id}"}), 400
        
        human_seq, human_name = fetch_sequence(human_id)
        if not human_seq:
            return jsonify({"error": f"Could not fetch sequence for human protein: {human_id}"}), 400
        
        # Clean sequences
        v_clean = clean_seq(viral_seq)
        h_clean = clean_seq(human_seq)
        
        if len(v_clean) < 5 or len(h_clean) < 5:
            return jsonify({"error": "Sequences too short after cleaning"}), 400
        
        # Extract features
        feat = aac(v_clean) + dpc(v_clean) + physchem(v_clean) + \
               aac(h_clean) + dpc(h_clean) + physchem(h_clean)
        
        X = pd.DataFrame([feat], columns=feature_cols)
        
        # Predict
        prob = float(model.predict_proba(X)[0, 1])
        pred = "INTERACTING" if prob >= 0.5 else "NON-INTERACTING"
        
        # Physicochemical properties for display
        v_pa = ProteinAnalysis(v_clean) if len(v_clean) >= 5 else None
        h_pa = ProteinAnalysis(h_clean) if len(h_clean) >= 5 else None
        
        result = {
            "prediction": pred,
            "probability": round(prob, 4),
            "confidence": round(abs(prob - 0.5) * 200, 1),
            "viral": {
                "id": viral_id,
                "name": viral_name or viral_id,
                "seq_length": len(v_clean),
                "mw": round(v_pa.molecular_weight(), 1) if v_pa else 0,
                "gravy": round(v_pa.gravy(), 3) if v_pa else 0,
            },
            "human": {
                "id": human_id,
                "name": human_name or human_id,
                "seq_length": len(h_clean),
                "mw": round(h_pa.molecular_weight(), 1) if h_pa else 0,
                "gravy": round(h_pa.gravy(), 3) if h_pa else 0,
            }
        }
        
        return jsonify(result)
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, port=5000)
