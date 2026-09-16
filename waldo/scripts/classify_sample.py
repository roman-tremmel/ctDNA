#!/usr/bin/env python3
"""Apply a trained SVM (train_svm_classifier.py) to classify a new sample
as euploid/aneuploid from its per-arm Z-scores.

Usage:
    classify_sample.py --scores sample.arm_scores.tsv --model model.joblib
"""
import argparse

import joblib
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", required=True, help="call_aneuploidy.py output TSV")
    ap.add_argument("--model", required=True)
    args = ap.parse_args()

    bundle = joblib.load(args.model)
    model, arm_columns = bundle["model"], bundle["arm_columns"]

    df = pd.read_csv(args.scores, sep="\t").set_index("arm")
    x = df["z_score"].reindex(arm_columns).fillna(0.0).to_frame().T

    label = model.predict(x)[0]
    proba = model.predict_proba(x)[0]
    classes = list(model.classes_)
    aneuploid_proba = proba[classes.index(1)] if 1 in classes else float("nan")

    print(f"Predicted label: {'aneuploid' if label == 1 else 'euploid'}")
    print(f"P(aneuploid) = {aneuploid_proba:.3f}")


if __name__ == "__main__":
    main()
