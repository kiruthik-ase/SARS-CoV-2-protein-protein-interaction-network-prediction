import pandas as pd
import os

print("\n=== BIOGRID VIRAL PROTEIN AUDIT (GROUND TRUTH) ===\n")

# =========================
# Paths
# =========================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIOGRID_FILE = os.path.join(BASE_DIR, "data", "raw", "BIOGRID-PROJECT-covid19_coronavirus_project-5.0.251", "BIOGRID-PROJECT-covid19_coronavirus_project-INTERACTIONS-5.0.251.tab3.txt")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# =========================
# 1. Load BioGRID file
# =========================
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

biogrid_viral_ids.to_csv(os.path.join(RESULTS_DIR, "biogrid_all_viral_uniprot_ids.csv"), index=False)
print("\nSaved → results/biogrid_all_viral_uniprot_ids.csv")

# =========================
# 4. Load final dataset to check viral proteins used
# =========================
DATASET_FILE = os.path.join(BASE_DIR, "data", "processed", "final_ppi_dataset.csv")

if os.path.exists(DATASET_FILE):
    dataset = pd.read_csv(DATASET_FILE, usecols=["viral_uniprot"])
    dataset_viral = set(dataset["viral_uniprot"].unique())
    biogrid_set = set(biogrid_viral_ids["viral_uniprot"])

    intersection = biogrid_set & dataset_viral
    missing = biogrid_set - dataset_viral
    extra = dataset_viral - biogrid_set

    print("\n--- FINAL TRUTH ---")
    print("BioGRID viral proteins       :", len(biogrid_set))
    print("Dataset viral proteins       :", len(dataset_viral))
    print("USABLE (intersection)        :", len(intersection))
    print("Missing from dataset         :", len(missing))
    print("In dataset but not BioGRID   :", len(extra))
else:
    print("\n⚠ final_ppi_dataset.csv not found — skipping intersection analysis")

print("\n=== AUDIT COMPLETE ===")
