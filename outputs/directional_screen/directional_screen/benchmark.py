"""
directional_screen.benchmark
============================

Loader + scorer for the packaged directional target-nomination benchmark.

The benchmark is a set of immune drug targets with a *directional* gold label
(driver_antagonize vs brake_agonize) derived from the mechanism of the approved
drug against each target. It is external ground truth INDEPENDENT of any
perturbation screen, so it can score whether a directional screen-scoring method
calls the therapeutically correct direction.

It ships WITH an honest negative: on the packaged set the reference framework's
directional accuracy is below the majority-class baseline. The point of releasing
it is to let others measure directional resolution ON THEIR OWN screen, not to
claim the method solves direction.
"""
from __future__ import annotations
import json
import os
import numpy as np
import pandas as pd

_HERE = os.path.dirname(__file__)
_DATA = os.path.join(_HERE, "data")


def load_benchmark(version: str = "v1") -> pd.DataFrame:
    """Load the packaged directional benchmark as a DataFrame."""
    path = os.path.join(_DATA, f"directional_benchmark_{version}.csv")
    return pd.read_csv(path)


def load_dictionary(version: str = "v1") -> dict:
    """Load the benchmark's data dictionary + honest-baseline metadata."""
    path = os.path.join(_DATA, f"directional_benchmark_{version}.dictionary.json")
    with open(path) as f:
        return json.load(f)


def evaluate_directions(
    calls: pd.Series,
    version: str = "v1",
) -> dict:
    """Score a set of directional calls against the benchmark.

    Parameters
    ----------
    calls : Series gene -> {'driver_antagonize','brake_agonize'}
        Your method's directional call per gene (any subset of benchmark genes).
    version : str

    Returns
    -------
    dict with raw accuracy, per-class accuracy, balanced accuracy, and a
    permutation p-value against the sign-shuffled null. Balanced accuracy and
    the permutation test are the honest metrics — raw accuracy is inflated by
    class imbalance (most drug targets are antagonize).
    """
    bench = load_benchmark(version).set_index("gene")
    shared = [g for g in calls.index if g in bench.index]
    if not shared:
        raise ValueError("No overlap between calls and benchmark genes.")
    truth = bench.loc[shared, "true_direction"]
    pred = calls.loc[shared]
    correct = (truth.to_numpy() == pred.to_numpy())
    # per-class
    accs = {}
    for cls in ["driver_antagonize", "brake_agonize"]:
        m = truth == cls
        accs[cls] = float(correct[m.to_numpy()].mean()) if m.any() else np.nan
    bal = float(np.nanmean(list(accs.values())))
    # permutation null on balanced accuracy
    rng = np.random.default_rng(0)
    labels = np.array(["driver_antagonize", "brake_agonize"])
    null = []
    tvals = truth.to_numpy()
    for _ in range(2000):
        shuf = rng.choice(labels, size=len(shared))
        c = (tvals == shuf)
        a = [c[(tvals == cls)].mean() for cls in labels if (tvals == cls).any()]
        null.append(np.nanmean(a))
    null = np.array(null)
    p = float((null >= bal).mean())
    return {
        "n_evaluated": len(shared),
        "raw_accuracy": float(correct.mean()),
        "per_class_accuracy": accs,
        "balanced_accuracy": bal,
        "permutation_p": p,
        "majority_baseline": float((truth == "driver_antagonize").mean()),
        "note": "Balanced accuracy + permutation p are the honest metrics; raw accuracy is inflated by antagonize-class imbalance.",
    }
