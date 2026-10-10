"""Predict labels and class probabilities for peptide/protein-pair CSV files."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from feature_extraction import extract_features, load_model


def predict_csv(input_path: Path, output_path: Path | None = None) -> Path:
    if not input_path.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    model = load_model()
    if not hasattr(model, "predict_proba"):
        raise TypeError("The loaded model does not provide predict_proba().")

    frame = pd.read_csv(input_path, keep_default_na=False)
    features = extract_features(frame, model)
    predictions = model.predict(features)
    probabilities = model.predict_proba(features)

    classes = np.asarray(model.classes_)
    if len(classes) != 2 or 0 not in classes or 1 not in classes:
        raise ValueError(f"Expected binary model classes 0 and 1; found {classes.tolist()}")
    class_0_index = int(np.flatnonzero(classes == 0)[0])
    class_1_index = int(np.flatnonzero(classes == 1)[0])
    predicted_indices = np.argmax(probabilities, axis=1)

    result = frame.copy()
    result["Y_pred"] = predictions
    result["probability_class_0"] = probabilities[:, class_0_index]
    result["probability_class_1"] = probabilities[:, class_1_index]
    result["Y_pred_probability"] = probabilities[
        np.arange(len(predictions)), predicted_indices
    ]

    destination = output_path or input_path.with_name(
        f"{input_path.stem}_predictions.csv"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(destination, index=False)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract model features and predict labels/probabilities from a CSV."
    )
    parser.add_argument("input_csv", type=Path, help="CSV with amp_sequence and protein_sequence")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output CSV path (default: <input>_predictions.csv)"
    )
    args = parser.parse_args()
    saved_path = predict_csv(args.input_csv, args.output)
    print(f"Predictions saved to: {saved_path.resolve()}")


if __name__ == "__main__":
    main()
