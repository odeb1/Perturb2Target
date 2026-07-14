"""
directional_screen.calibration
==============================

Calibrate a per-perturbation probability that the *directional* call is correct,
using an independent screen as external truth, and measure how well those
probabilities are calibrated (reliability) and how well they separate correct
from incorrect calls (resolution).

The lesson from the Perturb2Target project: a directional call from a single
screen has MODEST resolution. Concordance with an external screen rises with
effect magnitude but the per-gene direction is not sharply predictable. This
module makes that honesty quantitative rather than assumed — it reports the
reliability curve and the Brier skill score so a user of *their own* screen
learns the resolution ceiling on their data instead of trusting the call blindly.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def concordance_labels(
    scored: pd.DataFrame,
    external_signed: pd.Series,
    min_abs_external: float = 1.0,
) -> pd.DataFrame:
    """Label each perturbation by whether its direction matches an external screen.

    Parameters
    ----------
    scored : DataFrame from score_screen (must have 'program_score')
    external_signed : Series gene -> signed external effect (e.g. an independent
        CRISPRi cytokine z-score). Sign convention must match: negative = driver.
    min_abs_external : float
        Only genes whose |external| exceeds this are used as confident truth.

    Returns
    -------
    DataFrame (subset with confident external truth) with columns:
        program_score, external_signed, concordant (bool)
    """
    shared = scored.index.intersection(external_signed.index)
    df = pd.DataFrame({
        "program_score": scored.loc[shared, "program_score"],
        "external_signed": external_signed.loc[shared],
    }).dropna()
    df = df[df["external_signed"].abs() >= min_abs_external]
    df["concordant"] = (df["program_score"] < 0) == (df["external_signed"] < 0)
    return df


def fit_direction_confidence(
    labelled: pd.DataFrame,
    features: pd.DataFrame,
    n_splits: int = 5,
    seed: int = 0,
) -> pd.DataFrame:
    """Out-of-fold logistic calibration of P(direction correct).

    Parameters
    ----------
    labelled : DataFrame from concordance_labels (needs 'concordant')
    features : DataFrame indexed by perturbation, numeric predictors
        (e.g. abs(program_score), abs(emp_z), knockdown depth, n_downstream).
    n_splits : int   CV folds
    seed : int

    Returns
    -------
    DataFrame aligned to `labelled.index` with an added 'p_direction_correct'
    column (out-of-fold predicted probabilities).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler

    idx = labelled.index.intersection(features.index)
    y = labelled.loc[idx, "concordant"].astype(int).to_numpy()
    X = features.loc[idx].fillna(0.0).to_numpy(dtype=float)
    p = np.full(y.shape, np.nan)
    if len(np.unique(y)) < 2:
        # degenerate: all concordant or all discordant
        p[:] = y.mean()
    else:
        skf = StratifiedKFold(n_splits=min(n_splits, y.sum(), (1 - y).sum().astype(int) or 1),
                              shuffle=True, random_state=seed)
        for tr, te in skf.split(X, y):
            sc = StandardScaler().fit(X[tr])
            lr = LogisticRegression(max_iter=1000).fit(sc.transform(X[tr]), y[tr])
            p[te] = lr.predict_proba(sc.transform(X[te]))[:, 1]
    out = labelled.loc[idx].copy()
    out["p_direction_correct"] = p
    return out


def reliability(cal: pd.DataFrame, n_bins: int = 5) -> pd.DataFrame:
    """Reliability table: predicted vs observed concordance by probability quantile."""
    d = cal.dropna(subset=["p_direction_correct"]).copy()
    d["bin"] = pd.qcut(d["p_direction_correct"], q=min(n_bins, d["p_direction_correct"].nunique()),
                       duplicates="drop")
    tab = d.groupby("bin", observed=True).agg(
        predicted=("p_direction_correct", "mean"),
        observed=("concordant", "mean"),
        n=("concordant", "size"),
    ).reset_index(drop=True)
    return tab


def brier_skill(cal: pd.DataFrame) -> float:
    """Brier skill score vs the base-rate reference (0 = no skill, 1 = perfect)."""
    d = cal.dropna(subset=["p_direction_correct"])
    y = d["concordant"].astype(int).to_numpy()
    p = d["p_direction_correct"].to_numpy()
    bs = np.mean((p - y) ** 2)
    base = y.mean()
    bs_ref = np.mean((base - y) ** 2)
    return float(1 - bs / bs_ref) if bs_ref > 0 else 0.0
