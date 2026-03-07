# SARS-CoV-2 Protein–Protein Interaction Network Prediction

Machine learning pipeline to predict **protein-protein interactions (PPIs)** between **SARS-CoV-2 viral proteins** and **human host proteins** using sequence-based features and XGBoost.

## Overview

This project uses experimentally validated interaction data from [BioGRID](https://thebiogrid.org/) to train a binary classifier that predicts whether a given SARS-CoV-2 protein interacts with a human protein based solely on their amino acid sequences.

### Key Results
- **14 SARS-CoV-2 proteins** analyzed (Spike, Nucleocapsid, Envelope, Membrane, ORFs, NSPs)
- **~8,100 human proteins** in the interaction network
- **850 sequence-based features** per protein pair (AAC + DPC + physicochemical)
- **XGBoost classifier** with leave-one-viral-out cross-validation
- **Best ROC-AUC: 0.87** (ORF7a), **Mean ROC-AUC: 0.76**

## Project Structure

```
├── data/
│   ├── raw/                    # BioGRID interaction data, viral FASTA sequences
│   └── processed/              # Cleaned datasets ready for ML
├── src/
│   ├── pipeline.py             # End-to-end data pipeline (BioGRID → features → dataset)
│   ├── train.py                # XGBoost training with full evaluation
│   ├── train_human_split.py    # Training with human-wise split (no data leakage)
│   ├── eval_viral_cv.py        # Leave-one-viral-out cross-validation
│   ├── predict.py              # Prediction/inference script
│   └── audit.py                # BioGRID data audit utility
├── models/                     # Saved XGBoost model + feature definitions
├── results/                    # Feature importance, viral-wise evaluation results
└── website/                    # Interactive showcase website
```

## Setup

```bash
pip install -r requirements.txt
```

## Usage

### Run the full pipeline
```bash
python src/pipeline.py
```

### Train the model
```bash
python src/train.py
```

### Predict interactions
```bash
# Single prediction (by UniProt ID)
python src/predict.py --viral P0DTC2 --human Q9BYF1

# Batch prediction from CSV
python src/predict.py --batch pairs.csv --output predictions.csv
```

### Leave-one-viral-out evaluation
```bash
python src/eval_viral_cv.py
```

## Features

Each protein pair is represented by 850 features:

| Feature Type | Count (per protein) | Description |
|---|---|---|
| AAC | 20 | Amino acid composition (frequency of each of 20 amino acids) |
| DPC | 400 | Dipeptide composition (frequency of 400 possible dipeptide pairs) |
| Physicochemical | 5 | Length, molecular weight, GRAVY, aromaticity, instability index |

Features are computed for both the viral and human protein → **850 total** per pair.

## Data Sources

- **Interactions**: [BioGRID COVID-19 Project](https://thebiogrid.org/project/covid19) (v5.0.251)
- **Viral sequences**: NCBI RefSeq + UniProt mapping
- **Human sequences**: UniProt REST API

## Technology Stack

- Python 3.13
- XGBoost (gradient boosted trees)
- Biopython (sequence analysis)
- scikit-learn (evaluation metrics, splitting)
- pandas / numpy (data processing)
