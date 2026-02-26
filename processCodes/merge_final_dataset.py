import pandas as pd

print("\n=== STEP C: Merge positives and negatives ===")

pos = pd.read_csv("ml_classical_features.csv")
neg = pd.read_csv("ml_negative_features.csv")

print("Positive shape:", pos.shape)
print("Negative shape:", neg.shape)

# Concatenate
full = pd.concat([pos, neg], ignore_index=True)

# Shuffle
full = full.sample(frac=1, random_state=42).reset_index(drop=True)

full.to_csv("ml_final_dataset.csv", index=False)

print("Saved → ml_final_dataset.csv")
print("Final shape:", full.shape)
print(full["label"].value_counts())
