import pandas as pd

# Load features
pos = pd.read_csv(r"C:\Users\LENOVO\bio sem 4\csvs\ml_classical_features.csv")     # label = 1
neg = pd.read_csv(r"C:\Users\LENOVO\bio sem 4\csvs\ml_negative_features.csv")      # label = 0

# Combine
full = pd.concat([pos, neg], ignore_index=True)

print("Combined label distribution:")
print(full["label"].value_counts())

# Load human-wise split lists
train_humans = set(pd.read_csv("train_set.csv")["human_uniprot"])
test_humans  = set(pd.read_csv("test_set.csv")["human_uniprot"])

# Apply split
train = full[full["human_uniprot"].isin(train_humans)]
test  = full[full["human_uniprot"].isin(test_humans)]

print("\nFinal Train labels:")
print(train["label"].value_counts())

print("\nFinal Test labels:")
print(test["label"].value_counts())

# Save
train.to_csv("train_final.csv", index=False)
test.to_csv("test_final.csv", index=False)

print("\nSaved:")
print(" → train_final.csv")
print(" → test_final.csv")
