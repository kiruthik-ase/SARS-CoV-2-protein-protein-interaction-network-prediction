import pandas as pd
import numpy as np
import os

from xgboost import XGBClassifier

from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    classification_report,
    confusion_matrix
)

from sklearn.model_selection import StratifiedKFold

print("\n=== XGBOOST TRAINING: VIRAL–HUMAN PPI ===\n")

# -------------------------
# Paths
# -------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRAIN_FILE = os.path.join(BASE_DIR, "data", "processed", "train_final.csv")
TEST_FILE = os.path.join(BASE_DIR, "data", "processed", "test_final.csv")
MODEL_DIR = os.path.join(BASE_DIR, "models")
RESULTS_DIR = os.path.join(BASE_DIR, "results")

# -------------------------
# Load data
# -------------------------
train_df = pd.read_csv(TRAIN_FILE)
test_df  = pd.read_csv(TEST_FILE)

print("Train shape:", train_df.shape)
print("Test shape :", test_df.shape)

# -------------------------
# Separate features & labels
# -------------------------
drop_cols = ["viral_uniprot", "human_uniprot", "label"]

X_train = train_df.drop(columns=drop_cols)
y_train = train_df["label"]

X_test  = test_df.drop(columns=drop_cols)
y_test  = test_df["label"]

print("Feature count:", X_train.shape[1])

# -------------------------
# XGBoost model
# -------------------------
model = XGBClassifier(
    n_estimators=500,
    max_depth=8,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    objective="binary:logistic",
    eval_metric="auc",
    tree_method="hist",     # CPU optimized
    predictor="auto",
    random_state=42,
    n_jobs=-1
)

# -------------------------
# Train
# -------------------------
print("\nTraining XGBoost...")
model.fit(X_train, y_train)

# -------------------------
# Predict
# -------------------------
y_pred_prob = model.predict_proba(X_test)[:, 1]
y_pred = (y_pred_prob >= 0.5).astype(int)

# -------------------------
# Evaluation
# -------------------------
roc_auc = roc_auc_score(y_test, y_pred_prob)
pr_auc  = average_precision_score(y_test, y_pred_prob)
cm = confusion_matrix(y_test, y_pred)

print("\n=== EVALUATION RESULTS ===")
print(f"ROC-AUC : {roc_auc:.4f}")
print(f"PR-AUC  : {pr_auc:.4f}")
print("\nConfusion Matrix:")
print(cm)

print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# -------------------------
# Feature importance
# -------------------------
importance = model.feature_importances_
feat_names = X_train.columns

imp_df = pd.DataFrame({
    "feature": feat_names,
    "importance": importance
}).sort_values(by="importance", ascending=False)

imp_df.to_csv(os.path.join(RESULTS_DIR, "xgboost_feature_importance.csv"), index=False)

# -------------------------
# Save model
# -------------------------
model.save_model(os.path.join(MODEL_DIR, "xgboost_ppi_model.json"))

print("\nSaved:")
print(" → models/xgboost_ppi_model.json")
print(" → results/xgboost_feature_importance.csv")

print("\n=== TRAINING COMPLETE ===")
