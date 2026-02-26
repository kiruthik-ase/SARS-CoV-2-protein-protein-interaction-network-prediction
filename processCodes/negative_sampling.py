# negative_sampling.py
import pandas as pd
import random
import argparse
from collections import defaultdict

def random_negatives(pos_file, human_pool_file, neg_ratio=1, seed=42, out_file="neg_random.csv"):
    pos = pd.read_csv(pos_file)
    # get unique viral IDs from pos
    viral_ids = pos['viral_uniprot'].unique().tolist()
    # human pool: list of human_uniprot available (full sequences)
    human_df = pd.read_csv(human_pool_file)
    if 'uniprot' in human_df.columns:
        human_ids = human_df['uniprot'].unique().tolist()
    elif 'human_uniprot' in human_df.columns:
        human_ids = human_df['human_uniprot'].unique().tolist()
    else:
        raise ValueError("human pool file must have 'uniprot' or 'human_uniprot' column")

    existing = set(zip(pos['viral_uniprot'], pos['human_uniprot']))
    n_pos = pos.shape[0]
    n_neg = int(n_pos * neg_ratio)

    neg = set()
    random.seed(seed)
    attempts = 0
    while len(neg) < n_neg and attempts < n_neg * 20:
        v = random.choice(viral_ids)
        h = random.choice(human_ids)
        if (v,h) in existing or (v,h) in neg:
            attempts += 1
            continue
        neg.add((v,h))
        attempts += 1

    neg_df = pd.DataFrame(list(neg), columns=['viral_uniprot','human_uniprot'])
    neg_df['label'] = 0
    neg_df.to_csv(out_file, index=False)
    print("Saved negatives:", out_file, "count:", neg_df.shape[0])
    return out_file

def degree_preserving_negatives(pos_file, human_pool_file, neg_ratio=1, seed=42, out_file="neg_degree.csv"):
    pos = pd.read_csv(pos_file)
    human_df = pd.read_csv(human_pool_file)
    if 'uniprot' in human_df.columns:
        human_ids = human_df['uniprot'].unique().tolist()
    elif 'human_uniprot' in human_df.columns:
        human_ids = human_df['human_uniprot'].unique().tolist()
    else:
        raise ValueError("human pool file must have 'uniprot' or 'human_uniprot' column")

    # compute per-viral degree (#positives)
    viral_groups = pos.groupby('viral_uniprot')['human_uniprot'].nunique().to_dict()
    neg_rows = []
    random.seed(seed)
    for v, deg in viral_groups.items():
        want = int(deg * neg_ratio)
        sampled = 0
        attempts = 0
        while sampled < want and attempts < want * 20:
            h = random.choice(human_ids)
            if (v,h) in set(zip(pos['viral_uniprot'], pos['human_uniprot'])):
                attempts += 1
                continue
            neg_rows.append((v,h))
            sampled += 1
            attempts += 1
    neg_df = pd.DataFrame(neg_rows, columns=['viral_uniprot','human_uniprot'])
    neg_df['label'] = 0
    neg_df.to_csv(out_file, index=False)
    print("Saved degree-preserving negatives:", out_file, "count:", neg_df.shape[0])
    return out_file

if __name__=="__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pos_file", default="interactions_full_sequences_fixed.csv")
    p.add_argument("--human_pool", default="human_sequences_updated.csv")
    p.add_argument("--neg_ratio", type=float, default=1.0)
    p.add_argument("--method", choices=['random','degree'], default='random')
    p.add_argument("--out", default=None)
    args = p.parse_args()
    out = args.out or ("neg_"+args.method+".csv")
    if args.method == 'random':
        random_negatives(args.pos_file, args.human_pool, args.neg_ratio, out_file=out)
    else:
        degree_preserving_negatives(args.pos_file, args.human_pool, args.neg_ratio, out_file=out)
