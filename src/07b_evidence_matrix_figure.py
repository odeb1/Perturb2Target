#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 7 (figure) — top_targets_evidence_matrix.png: evidence heatmap x composite x modality/disease.

Conda environment: python
Produces artifact: top_targets_evidence_matrix.png
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mp
from matplotlib.colors import LinearSegmentedColormap

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


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


apply_figure_style()

full2 = pd.read_parquet("inputs/candidates_step6.parquet")

# Reconstruct df with integrated score
df = full2.copy()

df['causal_component'] = (df['abs_emp_z'].rank(pct=True)).clip(0, 1)

g = df['genetics_autoimmune'].fillna(0.0).clip(0, 1)
lof_bonus = df['lof_constrained'].fillna(False).astype(float) * 0.15
df['genetics_component'] = (g + lof_bonus).clip(0, 1)

trac_map = {'clinical_precedent': 1.0, 'discovery_precedent': 0.75, 'predicted_tractable': 0.55, 'difficult': 0.25}
df['drug_component'] = df['tractability_tier'].map(trac_map).fillna(0.15)
df.loc[df['suggested_modality'].str.contains('hard', na=False), 'drug_component'] = df['drug_component'].clip(upper=0.30)

nov_map = {'novel_undrugged': 1.0, 'preclinical_or_unknown': 0.7, 'not_assessed': 0.5, 'clinical_stage_drug': 0.35, 'approved_drug_exists': 0.15}
df['novelty_component'] = df['novelty_class'].map(nov_map).fillna(0.5)

W = dict(causal=0.34, genetics=0.30, drug=0.22, novelty=0.14)
df['integrated_score'] = (W['causal'] * df['causal_component'] + W['genetics'] * df['genetics_component']
                          + W['drug'] * df['drug_component'] + W['novelty'] * df['novelty_component'])

df['evidence_n'] = (df['program_hit'].astype(int) + (df['genetics_tier'].isin(['strong', 'moderate'])).astype(int)
                    + df['is_druggable'].fillna(False).astype(int) + (df['novelty_class'] == 'novel_undrugged').astype(int))
df = df.sort_values('integrated_score', ascending=False)
df['rank'] = np.arange(1, len(df) + 1)

df['novel_priority'] = ((df['novelty_class'].isin(['novel_undrugged', 'preclinical_or_unknown']))
                        & (df['genetics_tier'].isin(['strong', 'moderate']))
                        & (df['is_druggable'].fillna(False)))

top_overall = df.head(14).index.tolist()
top_novel = df[df['novel_priority']].head(16).index.tolist()
sel = list(dict.fromkeys(top_overall + top_novel))[:26]
M = df.loc[sel]
comp_cols = ['causal_component', 'genetics_component', 'drug_component', 'novelty_component']
comp_lab = ['causal\nstrength', 'human\ngenetics', 'drugg-\nability', 'clinical\nnovelty']
Z = M[comp_cols].values
N = len(sel)

fig = plt.figure(figsize=(13.6, 8.6))
gs = fig.add_gridspec(1, 3, width_ratios=[3.0, 1.0, 3.0], wspace=0.05)

axA = fig.add_subplot(gs[0, 0])
cmap = LinearSegmentedColormap.from_list('ev', ['#f7f7f7', '#c6dbef', '#4292c6', '#08306b'])
im = axA.imshow(Z, aspect='auto', cmap=cmap, vmin=0, vmax=1)
axA.set_xticks(range(4)); axA.set_xticklabels(comp_lab, fontsize=7.5); axA.xaxis.tick_top()
axA.set_yticks(range(N))
ylabs = []
for g in sel:
    r = M.loc[g]; mark = ''
    if r['novel_priority']: mark += ' ★'
    if r['in_paper_noms']: mark += '°'
    ylabs.append(f"$\\it{{{g}}}$" + mark)
axA.set_yticklabels(ylabs, fontsize=7)
for i in range(N):
    for j in range(4):
        v = Z[i, j]; axA.text(j, i, f"{v:.2f}", ha='center', va='center', fontsize=5.6, color='white' if v > 0.55 else '#333')
axA.set_title('Evidence components', fontsize=9.5, pad=22, loc='left')
panel_letter(axA, 'a')

axB = fig.add_subplot(gs[0, 1])
dir_col = {'driver_antagonize': '#b2182b', 'brake_agonize': '#2166ac'}
bcols = [dir_col[d] for d in M['direction']]
axB.barh(range(N), M['integrated_score'].values, color=bcols, height=0.72, alpha=0.9)
axB.set_ylim(N - 0.5, -0.5); axB.set_yticks([])
axB.set_xlim(0, 1.0); axB.set_xlabel('integrated\nscore', fontsize=7.5)
for i, v in enumerate(M['integrated_score'].values):
    axB.text(v + 0.02, i, f"{v:.2f}", va='center', fontsize=5.8, color='#333')
axB.set_title('Composite', fontsize=9, pad=22, loc='left')
axB.spines[['right', 'top', 'left']].set_visible(False)
panel_letter(axB, 'b')

axC = fig.add_subplot(gs[0, 2]); axC.axis('off')
axC.set_title('Modality · genetic disease anchor', fontsize=9, pad=22, loc='left')
mod_short = {'small molecule': 'SM', 'antibody/biologic (surface)': 'Ab', 'biologic (secreted)': 'Ab(sec)', 'hard (TF — degrader/PPI)': 'TF*', 'undetermined': '?'}
for i, g in enumerate(sel):
    r = M.loc[g]; mod = mod_short.get(r['suggested_modality'], '?')
    dis = r['ot_autoimmune_disease'] or r['ot_top_disease'] or '—'
    if isinstance(dis, str) and len(dis) > 34: dis = dis[:33] + '…'
    axC.text(0.0, i, mod, fontsize=6.3, va='center', fontweight='bold', color='#444')
    axC.text(0.13, i, dis, fontsize=6.0, va='center', color='#333')
axC.set_xlim(0, 1); axC.set_ylim(N - 0.5, -0.5)
axA.set_ylim(N - 0.5, -0.5)

cbar = fig.colorbar(im, ax=axC, fraction=0.025, pad=0.02, location='right')
cbar.set_label('component score (0–1)', fontsize=6.5); cbar.ax.tick_params(labelsize=6)
handles = [mp.Patch(color='#b2182b', label='antagonize driver'), mp.Patch(color='#2166ac', label='agonize brake')]
axB.legend(handles=handles, frameon=False, fontsize=6, loc='lower center', bbox_to_anchor=(0.5, -0.14))
fig.text(0.012, 0.02, '★ novel-priority nomination   ° paper-nominated   SM small molecule · Ab antibody · TF* transcription factor (hard)', fontsize=6, color='#555')
fig.suptitle('Top drug-target nominations — integrated causal + genetic + druggability + novelty evidence', fontsize=10.5, x=0.5, y=0.99)
fig.savefig("top_targets_evidence_matrix.png", dpi=200, bbox_inches='tight')
