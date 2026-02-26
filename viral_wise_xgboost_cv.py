import pandas as pd
import numpy as np

from xgboost import XGBClassifier
from sklearn.metrics import roc_auc_score, average_precision_score

print("\n=== VIRAL-WISE XGBOOST CROSS-VALIDATION ===\n")

# -------------------------
# Load dataset
# -------------------------
df = pd.read_csv(r"C:\Users\LENOVO\bio sem 4\final_ppi_dataset.csv")

print("Total samples:", df.shape[0])
print("Total viral proteins:", df["viral_uniprot"].nunique())

# -------------------------
# Columns
# -------------------------
id_cols = ["viral_uniprot", "human_uniprot", "label"]
feature_cols = [c for c in df.columns if c not in id_cols]

results = []

# -------------------------
# Leave-One-Viral-Out loop
# -------------------------
for virus in sorted(df["viral_uniprot"].unique()):
    print(f"\n--- Testing on viral protein: {virus} ---")

    train_df = df[df["viral_uniprot"] != virus]
    test_df  = df[df["viral_uniprot"] == virus]

    # Safety check
    if test_df["label"].nunique() < 2:
        print("⚠ Skipping (only one class present)")
        continue

    X_train = train_df[feature_cols]
    y_train = train_df["label"]

    X_test  = test_df[feature_cols]
    y_test  = test_df["label"]

    # -------------------------
    # Model
    # -------------------------
    model = XGBClassifier(
        n_estimators=500,
        max_depth=8,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        objective="binary:logistic",
        eval_metric="auc",
        tree_method="hist",
        random_state=42,
        n_jobs=-1
    )

    model.fit(X_train, y_train)

    # -------------------------
    # Evaluation
    # -------------------------
    y_prob = model.predict_proba(X_test)[:, 1]

    roc_auc = roc_auc_score(y_test, y_prob)
    pr_auc  = average_precision_score(y_test, y_prob)

    print(f"ROC-AUC: {roc_auc:.4f} | PR-AUC: {pr_auc:.4f}")

    results.append({
        "viral_uniprot": virus,
        "n_test_samples": len(test_df),
        "roc_auc": roc_auc,
        "pr_auc": pr_auc
    })

# -------------------------
# Final summary
# -------------------------
results_df = pd.DataFrame(results)
results_df.to_csv("viral_wise_results.csv", index=False)

print("\n=== FINAL SUMMARY ===")
print(results_df)

print("\nMean ROC-AUC:", results_df["roc_auc"].mean())
print("Mean PR-AUC :", results_df["pr_auc"].mean())

print("\nSaved → viral_wise_results.csv")
print("=== DONE ===")
