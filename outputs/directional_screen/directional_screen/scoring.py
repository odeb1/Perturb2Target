"""
directional_screen.scoring
===========================

Directional program scoring for massively-multiplexed perturbation screens.

The core idea: given a perturbation x readout matrix (z-scores), and a *signed*
gene program (each program gene weighted +1 if pro-program, -1 if anti-program),
score every perturbation by how strongly and in which *direction* it moves the
program. The sign of the score assigns each perturbation a therapeutic direction:

    negative score  ->  DRIVER  (knockdown lowers the program => antagonize the gene)
    positive score  ->  BRAKE   (knockdown raises the program => agonize the gene)

Statistical specificity is assessed against an empirical null of random signed
gene sets matched in composition, giving a per-perturbation empirical z and a
Benjamini-Hochberg FDR.

This module is assay-agnostic. It was developed on CD4+ T-cell CRISPRi
Perturb-seq (Zhu, Dann et al. 2025) but applies to any screen that yields a
perturbation x readout-gene effect matrix: CRISPRi/a Perturb-seq, pooled
cytokine-sort screens, or functional metaviromics screens where the "program"
is a signed host-response signature.

Author: Perturb2Target project (Claude Science hackathon).
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def program_score(
    zmatrix: pd.DataFrame,
    signature: pd.Series,
) -> pd.Series:
    """Signed program score per perturbation.

    Parameters
    ----------
    zmatrix : DataFrame  (perturbations x readout genes)
        Effect sizes (ideally z-scores) of each perturbation on each readout gene.
    signature : Series indexed by gene -> signed weight in {+1, -1}
        The directional program. Only genes present in both `signature.index`
        and `zmatrix.columns` are used.

    Returns
    -------
    Series indexed like `zmatrix.index` : the mean signed-weighted z over the
    program genes. Sign convention: negative = driver, positive = brake.
    """
    genes = [g for g in signature.index if g in zmatrix.columns]
    if not genes:
        raise ValueError("No signature genes found in zmatrix columns.")
    w = signature.loc[genes].to_numpy(dtype=float)
    Z = zmatrix[genes].to_numpy(dtype=float)
    return pd.Series((Z * w).mean(axis=1), index=zmatrix.index, name="program_score")


def empirical_null(
    zmatrix: pd.DataFrame,
    signature: pd.Series,
    n_null: int = 500,
    seed: int = 0,
) -> tuple[pd.Series, pd.Series]:
    """Empirical null of random *composition-matched* signed gene sets.

    Draws `n_null` random signatures with the same number of +1 and -1 genes as
    the real one, sampled from readout genes NOT in the real signature, and
    returns the per-perturbation mean and standard deviation of the resulting
    program scores. These convert a raw score into an empirical z-score:

        emp_z = (observed_score - null_mean) / null_std

    Parameters
    ----------
    zmatrix : DataFrame (perturbations x readout genes)
    signature : Series gene -> {+1,-1}
    n_null : int    number of random signed sets
    seed : int      RNG seed

    Returns
    -------
    (null_mean, null_std) : two Series indexed like zmatrix.index
    """
    rng = np.random.default_rng(seed)
    sig_genes = set(signature.index)
    pool = [g for g in zmatrix.columns if g not in sig_genes]
    n_pos = int((signature > 0).sum())
    n_neg = int((signature < 0).sum())
    draws = np.empty((n_null, zmatrix.shape[0]), dtype=float)
    Zpool = zmatrix[pool].to_numpy(dtype=float)
    pool_idx = np.arange(len(pool))
    for i in range(n_null):
        pick = rng.choice(pool_idx, size=n_pos + n_neg, replace=False)
        w = np.concatenate([np.ones(n_pos), -np.ones(n_neg)])
        draws[i] = (Zpool[:, pick] * w).mean(axis=1)
    return (pd.Series(draws.mean(0), index=zmatrix.index, name="null_mean"),
            pd.Series(draws.std(0), index=zmatrix.index, name="null_std"))


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. Returns q-values aligned to input order."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def score_screen(
    zmatrix: pd.DataFrame,
    signature: pd.Series,
    n_null: int = 500,
    seed: int = 0,
) -> pd.DataFrame:
    """End-to-end directional scoring of one screen (one condition).

    Returns a DataFrame indexed by perturbation with columns:
        program_score : signed mean-weighted z
        emp_z         : empirical z vs composition-matched null
        emp_p         : two-sided p from emp_z (normal approx)
        emp_fdr       : BH-corrected q-value
        direction     : 'driver_antagonize' (score<0) or 'brake_agonize' (score>0)
    """
    from scipy import stats
    s = program_score(zmatrix, signature)
    mu, sd = empirical_null(zmatrix, signature, n_null=n_null, seed=seed)
    emp_z = (s - mu) / sd.replace(0, np.nan)
    emp_p = pd.Series(2 * stats.norm.sf(emp_z.abs().to_numpy()), index=s.index)
    out = pd.DataFrame({
        "program_score": s,
        "emp_z": emp_z,
        "emp_p": emp_p,
        "emp_fdr": bh_fdr(emp_p.fillna(1.0).to_numpy()),
    })
    out["direction"] = np.where(out["program_score"] < 0,
                                "driver_antagonize", "brake_agonize")
    return out
