#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 3 — Directional, context-specific program scoring with a 500-draw empirical null;
BH-FDR per condition; direction + context-specificity; quality gate. -> perturbation_scores*.

Conda environment: python
Produces artifact: perturbation_scores_wide.parquet
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import io, time, requests, h5py, numpy as np, pandas as pd
from scipy import stats

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


# Load pre-computed data
obs = pd.read_parquet("{{artifact:de_obs_parquet}}")
var = pd.read_parquet("{{artifact:de_var_parquet}}")
sig = pd.read_csv("{{artifact:program_signature_csv}}")

Z = np.load("{{artifact:zscore_f32_npy}}", mmap_mode="r")

var_names = var['gene_name'].astype(str).values
gene2col = {g: i for i, g in enumerate(var_names)}
sig['in_readout'] = sig['gene'].isin(gene2col)
perts = set(obs['target_gene'].unique())
sig['is_perturbed'] = sig['gene'].isin(perts)

sig_meas = sig[sig['in_readout']].copy()
sig_meas['col'] = sig_meas['gene'].map(gene2col)

cols = sig_meas['col'].values
w = sig_meas['weight'].values.astype(np.float32)
genes_sig = sig_meas['gene'].values
arm_sig = sig_meas['arm'].values

Zsig = np.asarray(Z[:, cols])   # (33983, 48)

prog_score = (Zsig * w).mean(axis=1)
obs_s = obs.copy()
obs_s['prog_score'] = prog_score
pro_mask = (w > 0)
reg_mask = (w < 0)
obs_s['pro_mean_z'] = Zsig[:, pro_mask].mean(axis=1)
obs_s['reg_mean_z'] = Zsig[:, reg_mask].mean(axis=1)

WZ = Zsig * w
n_g = WZ.shape[1]
mean_wz = WZ.mean(axis=1)
std_wz = WZ.std(axis=1, ddof=1)
prog_t = mean_wz / (std_wz / np.sqrt(n_g) + 1e-9)
strong = np.abs(Zsig) > 2
concord = ((WZ > 0) & strong)
n_strong = strong.sum(axis=1)
frac_up = np.where(n_strong > 0, concord.sum(axis=1) / np.maximum(n_strong, 1), 0)

sc = obs_s[['target_gene', 'condition', 'row', 'prog_score', 'pro_mean_z', 'reg_mean_z',
            'n_downstream', 'n_cells_target', 'ontarget_effect_size']].copy()
sc['prog_t'] = prog_t
sc['n_sig_strong'] = n_strong
sc['kd_ok'] = sc['ontarget_effect_size'] < -1.0

wide = sc.pivot_table(index='target_gene', columns='condition',
                      values=['prog_score', 'prog_t', 'n_downstream', 'ontarget_effect_size'])
wide.columns = [f"{a}__{b}" for a, b in wide.columns]
ps = wide[[f'prog_score__{c}' for c in ['Rest', 'Stim8hr', 'Stim48hr']]].values
absps = np.abs(ps)
peak = absps.max(axis=1)
meanall = absps.mean(axis=1)
spec_index = np.where(peak > 1e-6, 1 - (meanall * 3 - peak) / (2 * peak + 1e-9), 0)
spec_index = np.clip(spec_index, 0, 1)
peak_cond = np.array(['Rest', 'Stim8hr', 'Stim48hr'])[absps.argmax(axis=1)]
wide['peak_abs_score'] = peak
wide['peak_condition'] = peak_cond
wide['context_specificity'] = spec_index
peak_signed = ps[np.arange(len(ps)), absps.argmax(axis=1)]
wide['peak_signed_score'] = peak_signed
wide['direction'] = np.where(peak_signed < 0, 'driver_antagonize', 'brake_agonize')

rng = np.random.default_rng(0)
B = 500
n_pos, n_neg = int(pro_mask.sum()), int(reg_mask.sum())
all_cols = np.arange(Z.shape[1])
noncols = np.setdiff1d(all_cols, cols)
null_scores = np.empty((Z.shape[0], B), dtype=np.float32)
for b in range(B):
    pos = rng.choice(noncols, n_pos, replace=False)
    neg = rng.choice(noncols, n_neg, replace=False)
    sc_b = (np.asarray(Z[:, pos]).sum(1) - np.asarray(Z[:, neg]).sum(1)) / n_g
    null_scores[:, b] = sc_b
    if (b + 1) % 100 == 0:
        print(b + 1, "draws")
null_mean = null_scores.mean(1)
null_std = null_scores.std(1, ddof=1)
emp_z = (obs_s['prog_score'].values - null_mean) / (null_std + 1e-9)

sc['emp_z'] = emp_z
sc['emp_p'] = 2 * stats.norm.sf(np.abs(sc['emp_z']))


def bh(p):
    p = np.asarray(p)
    n = len(p)
    order = np.argsort(p)
    ranked = np.empty(n)
    ranked[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(ranked, 0, 1)


sc['emp_fdr'] = sc.groupby('condition')['emp_p'].transform(lambda x: bh(x.values))

ez = sc.pivot_table(index='target_gene', columns='condition', values='emp_z')
ez.columns = [f'emp_z__{c}' for c in ez.columns]
fdr = sc.pivot_table(index='target_gene', columns='condition', values='emp_fdr')
fdr.columns = [f'emp_fdr__{c}' for c in fdr.columns]
kd = sc.pivot_table(index='target_gene', columns='condition', values='ontarget_effect_size')
kd.columns = [f'kd__{c}' for c in kd.columns]
W = wide.join(ez).join(fdr).join(kd)

peak_ez = W[[f'emp_z__{c}' for c in ['Rest', 'Stim8hr', 'Stim48hr']]].values
absps2 = np.abs(W[[f'prog_score__{c}' for c in ['Rest', 'Stim8hr', 'Stim48hr']]].values)
pk = absps2.argmax(1)
W['peak_emp_z'] = peak_ez[np.arange(len(W)), pk]
W['peak_emp_fdr'] = W[[f'emp_fdr__{c}' for c in ['Rest', 'Stim8hr', 'Stim48hr']]].values[np.arange(len(W)), pk]
W['peak_kd'] = W[[f'kd__{c}' for c in ['Rest', 'Stim8hr', 'Stim48hr']]].values[np.arange(len(W)), pk]
W['peak_ndown'] = W[[f'n_downstream__{c}' for c in ['Rest', 'Stim8hr', 'Stim48hr']]].values[np.arange(len(W)), pk]

W.to_parquet("perturbation_scores_wide.parquet")
