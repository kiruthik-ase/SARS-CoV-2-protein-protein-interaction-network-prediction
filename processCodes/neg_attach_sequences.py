import pandas as pd

print("\n=== STEP B1: Attach sequences to negatives ===")

neg = pd.read_csv("neg_random.csv")
viral = pd.read_csv("viral_sequences_mapped.csv")
human = pd.read_csv("human_sequences_clean.csv")

viral_seq = dict(zip(viral["viral_uniprot"], viral["viral_sequence"]))
human_seq = dict(zip(human["uniprot"], human["sequence"]))

neg["viral_sequence"] = neg["viral_uniprot"].map(viral_seq)
neg["human_sequence"] = neg["human_uniprot"].map(human_seq)

before = len(neg)
neg = neg.dropna(subset=["viral_sequence", "human_sequence"])
after = len(neg)

print(f"Negatives before: {before}")
print(f"Negatives after dropping missing seq: {after}")

neg.to_csv("negatives_with_sequences.csv", index=False)
print("Saved → negatives_with_sequences.csv")
print(neg.head())
