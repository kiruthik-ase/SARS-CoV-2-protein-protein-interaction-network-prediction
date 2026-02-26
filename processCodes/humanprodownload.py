import pandas as pd
import requests
from time import sleep

# Load your interaction file
df = pd.read_csv("clean_interactions.csv")

# Extract unique human UniProt IDs
human_ids = df["SWISS-PROT Accessions Interactor B"].unique()
print("Total unique human proteins =", len(human_ids))

sequences = []

for uid in human_ids:
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    try:
        r = requests.get(url)
        if r.status_code == 200 and r.text.startswith(">"):
            sequences.append({
                "uniprot": uid,
                "fasta": r.text
            })
            print("Fetched:", uid)
        else:
            print("Failed:", uid)
    except Exception as e:
        print("Error fetching", uid, ":", str(e))

    sleep(0.1)  # to avoid rate limits

# Save results
pd.DataFrame(sequences).to_csv("human_sequences.csv", index=False)

print("Saved: human_sequences.csv")
