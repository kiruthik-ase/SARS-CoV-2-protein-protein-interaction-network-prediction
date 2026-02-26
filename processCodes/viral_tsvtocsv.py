import pandas as pd
from Bio import SeqIO

# Full RefSeq → UniProt mapping including fragments
mapping = {
    # NSPs (all from ORF1ab → P0DTC1)
    "YP_009742608.1": "P0DTC1",  # nsp1
    "YP_009742609.1": "P0DTC1",  # nsp2
    "YP_009742610.1": "P0DTC1",  # nsp3
    "YP_009742611.1": "P0DTC1",  # nsp4
    "YP_009742612.1": "P0DTC1",  # nsp5
    "YP_009742613.1": "P0DTC1",  # nsp6
    "YP_009742614.1": "P0DTC1",  # nsp7
    "YP_009742615.1": "P0DTC1",  # nsp8
    "YP_009742616.1": "P0DTC1",  # nsp9
    "YP_009742617.1": "P0DTC1",  # nsp10
    "YP_009725311.1": "P0DTC1",  # nsp11
    "YP_009725312.1": "P0DTC1",  # nsp12
    "YP_009725313.1": "P0DTC1",  # nsp13
    "YP_009725314.1": "P0DTC1",  # nsp14
    "YP_009725315.1": "P0DTC1",  # nsp15
    "YP_009725316.1": "P0DTC1",  # nsp16

    # Structural proteins
    "YP_009724390.1": "P0DTC2",  # Spike
    "YP_009724391.1": "P0DTC3",  # ORF3a
    "YP_009724392.1": "P0DTC4",  # Envelope
    "YP_009724393.1": "P0DTC5",  # Membrane
    "YP_009724397.2": "P0DTC9",  # Nucleocapsid

    # Accessory proteins
    "YP_009724389.1": "Q7TLC7",     # ORF3b
    "YP_009724394.1": "P0DTC6",     # ORF6
    "YP_009724395.1": "P0DTC7",     # ORF7a
    "YP_009725318.1": "P0DTD8",     # ORF7b
    "YP_009724396.1": "P0DTC8",     # ORF8
    "YP_009725255.1": "A0A663DJA2", # ORF10

    # ORF1ab fragments → give stable dummy UniProt-like IDs
    "YP_009725295.1": "ORF1AB_F1",
    "YP_009725297.1": "ORF1AB_F2",
    "YP_009725298.1": "ORF1AB_F3",
    "YP_009725299.1": "ORF1AB_F4",
    "YP_009725300.1": "ORF1AB_F5",
    "YP_009725301.1": "ORF1AB_F6",
    "YP_009725302.1": "ORF1AB_F7",
    "YP_009725303.1": "ORF1AB_F8",
    "YP_009725304.1": "ORF1AB_F9",
    "YP_009725305.1": "ORF1AB_F10",
    "YP_009725306.1": "ORF1AB_F11",
    "YP_009725307.1": "ORF1AB_F12",
    "YP_009725308.1": "ORF1AB_F13",
    "YP_009725309.1": "ORF1AB_F14",
    "YP_009725310.1": "ORF1AB_F15",
}

virus_records = list(SeqIO.parse(r"C:\Users\LENOVO\bio sem 4\data\38ViralSequences.fasta", "fasta"))

rows = []
missing = []

for rec in virus_records:
    ref = rec.id
    seq = str(rec.seq)

    if ref not in mapping:
        missing.append(ref)
        continue

    rows.append({
        "viral_refseq": ref,
        "viral_uniprot": mapping[ref],
        "viral_sequence": seq
    })

df = pd.DataFrame(rows)
df.to_csv("viral_sequences_mapped.csv", index=False)

print(df)
print("\nTotal mapped:", df.shape[0])
print("Missing mappings:", missing)
