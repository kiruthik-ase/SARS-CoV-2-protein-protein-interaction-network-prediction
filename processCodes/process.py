import pandas as pd

print("\n=== FIXING HUMAN SEQUENCES ===")

human = pd.read_csv("human_sequences_updated.csv")

def fasta_to_seq(fasta):
    if not isinstance(fasta, str):
        return None
    lines = fasta.splitlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq

human["sequence"] = human["fasta"].apply(fasta_to_seq)

# Drop rows with empty sequences
human = human.dropna(subset=["sequence"])

human = human[["uniprot", "sequence"]]

human.to_csv("human_sequences_clean.csv", index=False)

print("Saved → human_sequences_clean.csv")
print("Rows:", len(human))
print(human.head())
