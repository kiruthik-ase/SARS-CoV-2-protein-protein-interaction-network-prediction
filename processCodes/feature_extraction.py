# feature_extraction.py
import pandas as pd
import numpy as np
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
import argparse
import os

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def clean_seq(s):
    if pd.isna(s): return ""
    s = str(s).upper()
    return re.sub(r"[^A-Z]", "", s)

def aa_composition(seq):
    L = max(1, len(seq))
    c = Counter(seq)
    return [c.get(aa,0)/L for aa in AA_LIST]

def dipeptide_composition(seq):
    pairs = [a+b for a in AA_LIST for b in AA_LIST]
    L = len(seq)
    denom = max(1, L-1)
    c = Counter([seq[i:i+2] for i in range(len(seq)-1)])
    return [c.get(p,0)/denom for p in pairs]

def physchem(seq):
    if len(seq) < 5:
        return [len(seq), 0.0, 0.0, 0.0, 0.0]
    try:
        pa = ProteinAnalysis(seq)
        return [len(seq), pa.molecular_weight(), pa.gravy(), pa.aromaticity(), pa.instability_index()]
    except Exception:
        return [len(seq), 0.0, 0.0, 0.0, 0.0]

def seq_features(seq):
    seq = clean_seq(seq)
    aac = aa_composition(seq)
    dpc = dipeptide_composition(seq)
    phys = physchem(seq)
    return np.array(aac + dpc + phys)

def build_features(infile, out_prefix, max_rows=None):
    print("Loading interactions:", infile)
    df = pd.read_csv(infile)
    if max_rows:
        df = df.head(max_rows).copy()
    # ensure columns present
    assert all(c in df.columns for c in ["viral_uniprot","viral_sequence","human_uniprot","human_sequence"])
    # drop rows missing sequences
    before = df.shape[0]
    df = df[df['viral_sequence'].notna() & df['human_sequence'].notna()].copy()
    print("Rows before:", before, "after dropping missing seq:", df.shape[0])

    # unique sequences to compute features once
    viral_seqs = df[['viral_uniprot','viral_sequence']].drop_duplicates().set_index('viral_uniprot')['viral_sequence'].to_dict()
    human_seqs = df[['human_uniprot','human_sequence']].drop_duplicates().set_index('human_uniprot')['human_sequence'].to_dict()

    print("Unique viral seqs:", len(viral_seqs), "unique human seqs:", len(human_seqs))

    # compute caches
    viral_feat = {}
    for k,s in viral_seqs.items():
        viral_feat[k] = seq_features(s)

    human_feat = {}
    for k,s in human_seqs.items():
        human_feat[k] = seq_features(s)

    # build feature rows
    rows = []
    for idx, r in df.iterrows():
        v = r['viral_uniprot']; h = r['human_uniprot']
        vf = viral_feat.get(v)
        hf = human_feat.get(h)
        if vf is None or hf is None:
            continue
        feat = np.concatenate([vf, hf])
        rows.append({
            'viral_uniprot': v,
            'human_uniprot': h,
            'label': 1,   # positives only here; negatives will be added later
            'features': feat
        })
    # Create DataFrame with columns
    feat_dim = len(next(iter(viral_feat.values()))) + len(next(iter(human_feat.values())))
    cols = []
    # viral names
    cols += [f"v_aac_{aa}" for aa in AA_LIST]
    cols += [f"v_dpc_{a}{b}" for a in AA_LIST for b in AA_LIST]
    cols += ["v_len","v_mw","v_gravy","v_arom","v_instab"]
    # human
    cols += [f"h_aac_{aa}" for aa in AA_LIST]
    cols += [f"h_dpc_{a}{b}" for a in AA_LIST for b in AA_LIST]
    cols += ["h_len","h_mw","h_gravy","h_arom","h_instab"]

    # expand rows into DataFrame
    X = np.vstack([r['features'] for r in rows])
    meta = pd.DataFrame([{'viral_uniprot': r['viral_uniprot'], 'human_uniprot': r['human_uniprot'], 'label': r['label']} for r in rows])
    feat_df = pd.DataFrame(X, columns=cols)
    out_df = pd.concat([meta.reset_index(drop=True), feat_df.reset_index(drop=True)], axis=1)

    out_path = f"{out_prefix}_classical_features.csv"
    out_df.to_csv(out_path, index=False)
    print("Saved:", out_path, "shape:", out_df.shape)
    return out_path

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--infile", default="interactions_full_sequences_fixed.csv")
    p.add_argument("--out_prefix", default="ml")
    p.add_argument("--max_rows", type=int, default=None)
    args = p.parse_args()
    build_features(args.infile, args.out_prefix, args.max_rows)
