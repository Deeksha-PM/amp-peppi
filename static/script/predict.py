"""Predict labels and class probabilities for peptide/protein-pair CSV files."""

from pathlib import Path

import numpy as np
import pandas as pd

from feature_extraction import extract_features, load_model
from validation import validate_frame

RESULT_COLUMNS = ["Y_pred", "probability_class_0", "probability_class_1", "Y_pred_probability"]


def predict_frame(frame, model, chunk_size=500, progress_callback=None) -> pd.DataFrame:
    """Predict for a frame with amp_sequence / protein_sequence columns.

    Works in chunks so callers can report progress via progress_callback(done, total).
    Returns a frame (RangeIndex) with the RESULT_COLUMNS, row-aligned with `frame`.
    """
    if not hasattr(model, "predict_proba"):
        raise TypeError("The loaded model does not provide predict_proba().")

    classes = np.asarray(model.classes_)
    if len(classes) != 2 or 0 not in classes or 1 not in classes:
        raise ValueError(f"Expected binary model classes 0 and 1; found {classes.tolist()}")
    class_0_index = int(np.flatnonzero(classes == 0)[0])
    class_1_index = int(np.flatnonzero(classes == 1)[0])

    frame = frame.reset_index(drop=True)
    total = len(frame)
    if total == 0:
        raise ValueError("There are no sequence pairs to predict.")

    parts = []
    for start in range(0, total, chunk_size):
        chunk = frame.iloc[start:start + chunk_size].reset_index(drop=True)
        features = extract_features(chunk, model)
        predictions = model.predict(features)
        probabilities = model.predict_proba(features)
        predicted_indices = np.argmax(probabilities, axis=1)
        parts.append(pd.DataFrame({
            "Y_pred": predictions,
            "probability_class_0": probabilities[:, class_0_index],
            "probability_class_1": probabilities[:, class_1_index],
            "Y_pred_probability": probabilities[np.arange(len(predictions)), predicted_indices],
        }))
        if progress_callback:
            progress_callback(min(start + chunk_size, total), total)
    return pd.concat(parts, ignore_index=True)
