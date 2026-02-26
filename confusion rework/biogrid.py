import pandas as pd

print("\n=== BIOGRID VIRAL PROTEIN AUDIT (GROUND TRUTH) ===\n")

# =========================
# 1. Load BioGRID file
# =========================
# ⚠️ Change filename if needed
BIOGRID_FILE = r"C:\Users\LENOVO\bio sem 4\Raw Data\BIOGRID-PROJECT-covid19_coronavirus_project-5.0.251\BIOGRID-PROJECT-covid19_coronavirus_project-INTERACTIONS-5.0.251.tab3.txt"

biogrid = pd.read_csv(BIOGRID_FILE, sep="\t", low_memory=False)

print("Total BioGRID rows:", len(biogrid))

# =========================
# 2. Filter SARS-CoV-2 interactions
# =========================
# SARS-CoV-2 taxonomy ID = 2697049
viral_rows = biogrid[biogrid["Organism ID Interactor A"] == 2697049]

print("Total SARS-CoV-2 interaction rows:", len(viral_rows))

# =========================
# 3. Extract viral protein identifiers
# =========================
viral_uniprot = (
    viral_rows["SWISS-PROT Accessions Interactor A"]
    .dropna()
    .str.split("|")
    .str[0]
)

viral_refseq = (
    viral_rows["REFSEQ Accessions Interactor A"]
    .dropna()
    .str.split("|")
    .str[0]
)

viral_systematic = viral_rows["Systematic Name Interactor A"].dropna()

print("\n--- UNIQUE VIRAL PROTEINS IN BIOGRID ---")
print("Unique UniProt IDs  :", viral_uniprot.nunique())
print("Unique RefSeq IDs   :", viral_refseq.nunique())
print("Unique Systematic   :", viral_systematic.nunique())

# Save ALL viral UniProt IDs from BioGRID
biogrid_viral_ids = pd.DataFrame({
    "viral_uniprot": viral_uniprot.unique()
})

biogrid_viral_ids.to_csv("biogrid_all_viral_uniprot_ids.csv", index=False)
print("\nSaved → biogrid_all_viral_uniprot_ids.csv")

# =========================
# 4. Load viral sequence file
# =========================
# ⚠️ Change filename if needed
VIRAL_SEQ_FILE = r"C:\Users\LENOVO\bio sem 4\csvs\viral_sequences.csv"

viral_seq = pd.read_csv(VIRAL_SEQ_FILE)

if "viral_uniprot" not in viral_seq.columns:
    raise ValueError("viral_sequences.csv must contain 'viral_uniprot' column")

print("\nTotal viral sequences available:",
      viral_seq["viral_uniprot"].nunique())

# =========================
# 5. Intersection analysis
# =========================
biogrid_set = set(biogrid_viral_ids["viral_uniprot"])
sequence_set = set(viral_seq["viral_uniprot"])

intersection = biogrid_set & sequence_set
missing_seq = biogrid_set - sequence_set
extra_seq = sequence_set - biogrid_set

print("\n--- FINAL TRUTH ---")
print("BioGRID viral proteins       :", len(biogrid_set))
print("Sequence viral proteins      :", len(sequence_set))
print("USABLE (intersection)        :", len(intersection))
print("Missing sequences            :", len(missing_seq))
print("Sequences not in BioGRID     :", len(extra_seq))

# =========================
# 6. Save missing & extra lists
# =========================
pd.DataFrame({"missing_uniprot": list(missing_seq)}).to_csv(
    "missing_viral_sequences.csv", index=False
)

pd.DataFrame({"extra_uniprot": list(extra_seq)}).to_csv(
    "extra_viral_sequences.csv", index=False
)

print("\nSaved:")
print(" → missing_viral_sequences.csv")
print(" → extra_viral_sequences.csv")

print("\n=== AUDIT COMPLETE ===")
