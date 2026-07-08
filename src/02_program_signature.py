#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2 — Define the pathogenic pro-inflammatory program signature (36 up / 12 down)
and render program_signature_heatmap.png across positive-control perturbations.

Conda environment: python
Produces artifact: program_signature_heatmap.png
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

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


apply_figure_style()

obs = pd.read_parquet("inputs/de_obs.parquet")
var = pd.read_parquet("inputs/de_var.parquet")
Z = np.load("inputs/zscore_f32.npy", mmap_mode="r")
sig = pd.read_csv("inputs/program_signature.csv")

var_names = var['gene_name'].astype(str).values
gene2col = {g: i for i, g in enumerate(var_names)}
sig['in_readout'] = sig['gene'].isin(gene2col)

sig_meas = sig[sig['in_readout']].copy()
sig_meas['col'] = sig_meas['gene'].map(gene2col)

cols = sig_meas['col'].values
w = sig_meas['weight'].values.astype(np.float32)
genes_sig = sig_meas['gene'].values
arm_sig = sig_meas['arm'].values

Zsig = np.asarray(Z[:, cols])

prog_score = (Zsig * w).mean(axis=1)
obs_s = obs.copy()
obs_s['prog_score'] = prog_score
pro_mask = (w > 0)
reg_mask = (w < 0)
obs_s['pro_mean_z'] = Zsig[:, pro_mask].mean(axis=1)
obs_s['reg_mean_z'] = Zsig[:, reg_mask].mean(axis=1)

drivers = ['TBX21', 'BATF', 'IRF4', 'STAT3', 'RORC', 'STAT4', 'NFKB1', 'REL', 'IL2RA']
brakes  = ['STAT6', 'IL4R', 'INPP5D', 'FOXP3', 'CTLA4', 'PTPN2', 'SOCS1', 'SOCS3', 'CBLB']
qc_genes = drivers + brakes

cond = 'Stim48hr'
sub = obs_s[obs_s['condition'] == cond].set_index('target_gene')
qc_present = [g for g in qc_genes if g in sub.index]

order = np.argsort(~pro_mask)
col_order = order
genes_ord = genes_sig[col_order]
arm_ord = arm_sig[col_order]

rowidx = [sub.loc[g, 'row'] if np.ndim(sub.loc[g, 'row']) == 0 else sub.loc[g, 'row'].iloc[0] for g in qc_present]
M = np.asarray(Z[np.array(rowidx)][:, cols[col_order]])

fig, ax = plt.subplots(figsize=(13, 5.6))
M_clip = np.clip(M, -8, 8)
norm = TwoSlopeNorm(vmin=-8, vcenter=0, vmax=8)
im = ax.imshow(M_clip, aspect='auto', cmap='RdBu_r', norm=norm)
ax.set_xticks(range(len(genes_ord)))
ax.set_xticklabels([f"$\\it{{{g}}}$" for g in genes_ord], rotation=90, fontsize=6)
ax.set_yticks(range(len(qc_present)))
ax.set_yticklabels([f"$\\it{{{g}}}$" for g in qc_present], fontsize=7.5)
n_pro = int(pro_mask.sum())
ax.axvline(n_pro - 0.5, color='black', lw=1.4)
ax.axhline(len(drivers) - 0.5, color='black', lw=1.4)
ax.text((n_pro - 1) / 2, -1.1, 'pro-inflammatory program genes', ha='center', va='bottom', fontsize=8, fontweight='bold')
ax.text(n_pro + (len(genes_ord) - n_pro - 1) / 2, -1.1, 'regulatory genes', ha='center', va='bottom', fontsize=8, fontweight='bold')
ax.text(-6.5, (len(drivers) - 1) / 2, 'known\ndrivers', ha='right', va='center', fontsize=8, fontweight='bold', color='#b2182b')
ax.text(-6.5, len(drivers) + (len(qc_present) - len(drivers) - 1) / 2, 'known\nbrakes', ha='right', va='center', fontsize=8, fontweight='bold', color='#2166ac')
cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, extend='both')
cbar.set_label('perturbation effect (z-score)\nblue = gene down on KD', fontsize=7)
fig.suptitle('Signature QC: knockdown of canonical regulators moves the program coherently (Stim 48 hr)',
             fontsize=9.5, x=0.5, y=1.02)
fig.tight_layout()
fig.savefig("program_signature_heatmap.png", dpi=200, bbox_inches='tight')
