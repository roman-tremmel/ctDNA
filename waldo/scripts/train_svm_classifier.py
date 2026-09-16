#!/usr/bin/env python3
"""Train an SVM to classify samples as euploid/aneuploid from their
per-chromosome-arm Z-scores, mirroring the genome-wide classification
layer described in Douville et al. 2018 (used there to catch samples
with many individually-nonsignificant arm perturbations, e.g. low
neoplastic-fraction liquid biopsies).

The paper reports testing several ML algorithms and finding SVM optimal,
trained on tens of thousands of synthetic samples; this script provides
the same mechanism (scikit-learn SVC) so you can retrain it on your own
labeled data (synthetic and/or real) - it does not ship a pretrained
model, since the paper's is not published.

Usage:
    train_svm_classifier.py \
        --features feature_matrix.tsv \
        --labels labels.tsv \
        --out model.joblib
"""
import argparse

import joblib
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--features", required=True,
                     help="samples x arms Z-score matrix (build_feature_matrix.py output)")
    ap.add_argument("--labels", required=True,
                     help="TSV with columns: sample, label (0=euploid, 1=aneuploid)")
    ap.add_argument("--kernel", default="rbf", choices=["linear", "rbf", "poly"])
    ap.add_argument("--cv-folds", type=int, default=5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    X = pd.read_csv(args.features, sep="\t").set_index("sample")
    labels = pd.read_csv(args.labels, sep="\t").set_index("sample")["label"]
    labels = labels.loc[X.index]

    model = make_pipeline(
        StandardScaler(),
        SVC(kernel=args.kernel, probability=True, class_weight="balanced"),
    )

    if len(X) >= args.cv_folds:
        scores = cross_val_score(model, X.fillna(0.0), labels, cv=args.cv_folds)
        print(f"{args.cv_folds}-fold CV accuracy: {scores.mean():.3f} +/- {scores.std():.3f}")
    else:
        print(f"WARNING: only {len(X)} training samples (< {args.cv_folds} CV folds); "
              "skipping cross-validation. Train on more samples before trusting this model.")

    model.fit(X.fillna(0.0), labels)
    joblib.dump({"model": model, "arm_columns": list(X.columns)}, args.out)
    print(f"Wrote trained model to {args.out}")


if __name__ == "__main__":
    main()
