import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt

from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
    ConfusionMatrixDisplay,
    RocCurveDisplay,
    PrecisionRecallDisplay
)

print("\n=== XGBOOST TRAINING WITH FULL EVALUATION ===\n")

# -------------------------
# Load dataset
# -------------------------
df = pd.read_csv("final_ppi_dataset.csv")

id_cols = ["viral_uniprot", "human_uniprot", "label"]
feature_cols = [c for c in df.columns if c not in id_cols]

X = df[feature_cols]
y = df["label"]

print("Total samples :", X.shape[0])
print("Total features:", X.shape[1])
print("Positive ratio:", y.mean())

# -------------------------
# Train / Test split
# -------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    stratify=y,
    random_state=42
)

print("\nTrain samples:", X_train.shape[0])
print("Test samples :", X_test.shape[0])

# -------------------------
# Model
# -------------------------
model = XGBClassifier(
    n_estimators=600,
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

# -------------------------
# Training
# -------------------------
print("\nTraining XGBoost...")
model.fit(X_train, y_train)

# -------------------------
# Predictions
# -------------------------
y_prob = model.predict_proba(X_test)[:, 1]
y_pred = (y_prob >= 0.5).astype(int)

# -------------------------
# Metrics
# -------------------------
roc_auc = roc_auc_score(y_test, y_prob)
pr_auc  = average_precision_score(y_test, y_prob)
acc     = accuracy_score(y_test, y_pred)
prec    = precision_score(y_test, y_pred)
rec     = recall_score(y_test, y_pred)
f1      = f1_score(y_test, y_pred)

print("\n=== EVALUATION METRICS ===")
print(f"ROC-AUC      : {roc_auc:.4f}")
print(f"PR-AUC       : {pr_auc:.4f}")
print(f"Accuracy     : {acc:.4f}")
print(f"Precision    : {prec:.4f}")
print(f"Recall       : {rec:.4f}")
print(f"F1-score     : {f1:.4f}")

# -------------------------
# Confusion Matrix
# -------------------------
cm = confusion_matrix(y_test, y_pred)
disp = ConfusionMatrixDisplay(cm)
disp.plot(cmap="Blues")
plt.title("Confusion Matrix")
plt.show()

# -------------------------
# ROC Curve
# -------------------------
RocCurveDisplay.from_predictions(y_test, y_prob)
plt.title("ROC Curve")
plt.show()

# -------------------------
# Precision-Recall Curve
# -------------------------
PrecisionRecallDisplay.from_predictions(y_test, y_prob)
plt.title("Precision-Recall Curve")
plt.show()

# -------------------------
# Feature Importance (Top 20)
# -------------------------
importances = model.feature_importances_
imp_df = pd.DataFrame({
    "feature": feature_cols,
    "importance": importances
}).sort_values(by="importance", ascending=False)

top20 = imp_df.head(20)

plt.figure(figsize=(8, 6))
plt.barh(top20["feature"][::-1], top20["importance"][::-1])
plt.title("Top 20 Feature Importances")
plt.xlabel("Importance")
plt.tight_layout()
plt.show()

# -------------------------
# Save model & features
# -------------------------
model.get_booster().save_model("ppi_xgboost_model.json")
joblib.dump(feature_cols, "feature_columns.pkl")
imp_df.to_csv("feature_importance_full.csv", index=False)

print("\nSaved files:")
print(" → ppi_xgboost_model.json")
print(" → feature_columns.pkl")
print(" → feature_importance_full.csv")
print("\n=== DONE ===")
