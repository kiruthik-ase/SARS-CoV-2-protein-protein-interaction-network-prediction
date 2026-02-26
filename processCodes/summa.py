import pandas as pd

# Load data
inter = pd.read_csv(r"C:\Users\LENOVO\bio sem 4\clean_interactions_fixed.csv")    # viral_uniprot, human_uniprot, ...
virus = pd.read_csv(r"C:\Users\LENOVO\bio sem 4\viral_sequences_mapped.csv")      # viral_uniprot, viral_sequence
human = pd.read_csv(r"C:\Users\LENOVO\bio sem 4\human_sequences_updated.csv")     # uniprot, fasta (cleaned seq)

# Normalize column names
human = human.rename(columns={"uniprot": "human_uniprot", "fasta": "human_sequence"})
virus = virus.rename(columns={"viral_uniprot": "viral_uniprot", "viral_sequence": "viral_sequence"})

print("INPUT SHAPES")
print("interactions:", inter.shape)
print("viral seq table:", virus.shape)
print("human seq table:", human.shape)

# Remove rows with human_uniprot == '-' from interactions
num_before = inter.shape[0]
inter = inter[inter["human_uniprot"].astype(str) != "-"].copy()
num_after = inter.shape[0]
removed_dash = num_before - num_after
print(f"\nRemoved rows with human_uniprot == '-': {removed_dash}")

# Deduplicate sequence tables: keep first occurrence per UniProt
human_before = human.shape[0]
human = human.drop_duplicates(subset=["human_uniprot"], keep="first").copy()
human_after = human.shape[0]
print(f"Human sequences deduped: {human_before} -> {human_after} (dropped {human_before-human_after})")

virus_before = virus.shape[0]
virus = virus.drop_duplicates(subset=["viral_uniprot"], keep="first").copy()
virus_after = virus.shape[0]
print(f"Viral sequences deduped: {virus_before} -> {virus_after} (dropped {virus_before-virus_after})")

# Check which human IDs are still missing
inter_human_ids = set(inter["human_uniprot"].astype(str).unique())
human_ids_present = set(human["human_uniprot"].astype(str).unique())
missing_after = sorted(list(inter_human_ids - human_ids_present))
print("\nDistinct human UniProt IDs referenced in interactions:", len(inter_human_ids))
print("Distinct human UniProt IDs present after dedupe:", len(human_ids_present))
print("Missing human IDs (will cause dropped rows):", missing_after)
print("Count missing:", len(missing_after))

# Drop interaction rows referencing missing human IDs and record how many
rows_before = inter.shape[0]
if missing_after:
    inter = inter[~inter["human_uniprot"].isin(missing_after)].copy()
rows_after = inter.shape[0]
dropped_for_missing = rows_before - rows_after
print(f"Dropped {dropped_for_missing} interaction rows due to missing human sequences (after dedupe).")

# Now safe merges
m1 = inter.merge(virus, on="viral_uniprot", how="left", validate="m:1")
print("\nAfter merging viral seq - shape:", m1.shape, "missing viral_sequence:", m1['viral_sequence'].isna().sum())

m2 = m1.merge(human, on="human_uniprot", how="left", validate="m:1")
print("After merging human seq - shape:", m2.shape, "missing human_sequence:", m2['human_sequence'].isna().sum())

# Keep only rows with both sequences
final = m2[m2["viral_sequence"].notna() & m2["human_sequence"].notna()].copy()
print("\nFinal rows with both sequences:", final.shape)

# Save final merged file
final.to_csv("interactions_full_sequences_fixed.csv", index=False)

# Save diagnostics
with open("merge_diagnostics.txt", "w", encoding="utf-8") as f:
    f.write("Merge diagnostics\n")
    f.write(f"Initial interactions: {inter.shape[0] + dropped_for_missing}\n")
    f.write(f"Removed '-' rows: {removed_dash}\n")
    f.write(f"Rows dropped for missing human sequences after dedupe: {dropped_for_missing}\n")
    f.write(f"Final rows with both sequences: {final.shape[0]}\n")
    f.write("Missing human IDs:\n")
    for mid in missing_after:
        f.write(str(mid) + "\n")

print("\n--- MERGE SUMMARY ---")
print("Initial interactions (after '-' removal):", rows_before)
print("Rows dropped due to missing human seq:", dropped_for_missing)
print("Final interactions with both sequences:", final.shape[0])
print("Distinct missing human IDs (listed in merge_diagnostics.txt):", len(missing_after))
print("Missing IDs sample:", missing_after[:20])
print("\nSaved: interactions_full_sequences_fixed.csv")
print("Saved: merge_diagnostics.txt")
