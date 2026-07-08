#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 8 — Deep-dive 6 novel nominations: context-specific program trajectories with empirical-null
significance. -> top_target_vignettes.png, deep_dive_vignettes.csv.

Conda environment: python
Produces artifact: top_target_vignettes.png
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as ml

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


# Load data
full2 = pd.read_parquet("inputs/candidates_step6.parquet")

# Reconstruct druggable class map
def load_list(path):
    return set(pd.read_csv(path, header=None)[0].astype(str).str.strip())

def tb(v):
    return bool(v) if (v is True or v is False) else (False if pd.isna(v) else bool(v))

def num(v):
    return 0 if pd.isna(v) else float(v)

# Reconstruct tractability tier and other derived columns if not present
easy_sm = {'kinase', 'GPCR', 'enzyme', 'ion_channel', 'nuclear_receptor', 'transporter', 'catalytic_receptor'}
surface = {'cytokine_receptor', 'catalytic_receptor', 'GPCR', 'ion_channel', 'transporter'}

def tract_tier(r):
    if tb(r.get('sm_clinical')) or tb(r.get('ab_clinical')) or num(r.get('n_drugs')) > 0:
        return 'clinical_precedent'
    if tb(r.get('sm_ligand')) or tb(r.get('sm_pocket')): return 'sm_tractable'
    if tb(r.get('ab_surface')): return 'ab_tractable'
    pc = r.get('primary_class')
    if pc in easy_sm: return 'class_tractable_sm'
    if pc in ('cytokine', 'cytokine_receptor'): return 'class_tractable_bio'
    if pc == 'transcription_factor': return 'difficult_tf'
    return 'unknown'

sm_classes = {'kinase', 'enzyme', 'GPCR', 'ion_channel', 'nuclear_receptor', 'transporter'}

def modality3(r):
    pc = r.get('primary_class')
    if pc == 'cytokine': return 'biologic (secreted)'
    if pc in ('cytokine_receptor', 'catalytic_receptor'): return 'antibody/biologic (surface)'
    if pc in sm_classes: return 'small molecule'
    if tb(r.get('sm_ligand')) or tb(r.get('sm_pocket')): return 'small molecule'
    if tb(r.get('ab_surface')): return 'antibody/biologic (surface)'
    if pc == 'transcription_factor': return 'hard (TF — degrader/PPI)'
    return 'undetermined'

# Reconstruct stage ranking
stage_rank = {'APPROVAL': 4, 'PHASE_3': 3, 'PHASE_2_3': 3, 'PHASE_2': 2, 'PHASE_1_2': 2, 'PHASE_1': 1, 'UNKNOWN': 0}

def max_stage(s):
    if not isinstance(s, str) or not s: return np.nan
    vals = [stage_rank.get(t.strip()) for t in s.split(';') if t.strip() in stage_rank]
    return max(vals) if vals else np.nan

full2['max_drug_stage'] = full2['drug_stages'].apply(max_stage)

def novelty(r):
    ns = r.get('n_drugs'); ms = r.get('max_drug_stage')
    if pd.isna(ns): return 'not_assessed'
    if ns == 0: return 'novel_undrugged'
    if ms is not None and not pd.isna(ms) and ms >= 4: return 'approved_drug_exists'
    if ms is not None and not pd.isna(ms) and ms >= 1: return 'clinical_stage_drug'
    return 'preclinical_or_unknown'

full2['novelty_class'] = full2.apply(novelty, axis=1)

paper_noms = {'STAT3', 'IRF4', 'BATF', 'LRRC25', 'FAM20B'}
full2['in_paper_noms'] = full2.index.isin(paper_noms)

# Reconstruct tractability and modality
full2['tractability_tier'] = full2.apply(tract_tier, axis=1)
full2['suggested_modality'] = full2.apply(modality3, axis=1)
full2['is_druggable'] = ~full2['tractability_tier'].isin(['difficult_tf', 'unknown'])

# Reconstruct integrated score
full2['causal_component'] = (full2['abs_emp_z'].rank(pct=True)).clip(0, 1)
g = full2['genetics_autoimmune'].fillna(0.0).clip(0, 1)
lof_bonus = full2['lof_constrained'].fillna(False).astype(float) * 0.15
full2['genetics_component'] = (g + lof_bonus).clip(0, 1)

trac_map = {'clinical_precedent': 1.0, 'discovery_precedent': 0.75, 'predicted_tractable': 0.55, 'difficult': 0.25}
full2['drug_component'] = full2['tractability_tier'].map(trac_map).fillna(0.15)
full2.loc[full2['suggested_modality'].str.contains('hard', na=False), 'drug_component'] = full2['drug_component'].clip(upper=0.30)

nov_map = {'novel_undrugged': 1.0, 'preclinical_or_unknown': 0.7, 'not_assessed': 0.5, 'clinical_stage_drug': 0.35, 'approved_drug_exists': 0.15}
full2['novelty_component'] = full2['novelty_class'].map(nov_map).fillna(0.5)

W = dict(causal=0.34, genetics=0.30, drug=0.22, novelty=0.14)
full2['integrated_score'] = (W['causal'] * full2['causal_component'] + W['genetics'] * full2['genetics_component']
                             + W['drug'] * full2['drug_component'] + W['novelty'] * full2['novelty_component'])

full2['evidence_n'] = (full2['program_hit'].astype(int) + (full2['genetics_tier'].isin(['strong', 'moderate'])).astype(int)
                       + full2['is_druggable'].fillna(False).astype(int) + (full2['novelty_class'] == 'novel_undrugged').astype(int))
full2 = full2.sort_values('integrated_score', ascending=False)
full2['rank'] = np.arange(1, len(full2) + 1)

full2['novel_priority'] = ((full2['novelty_class'].isin(['novel_undrugged', 'preclinical_or_unknown']))
                           & (full2['genetics_tier'].isin(['strong', 'moderate']))
                           & (full2['is_druggable'].fillna(False)))

df = full2

# Deep-dive picks
picks = ['STAT5B', 'CTSH', 'NEK6', 'PTPN2', 'SIK2', 'STAT6']

show = ['direction', 'integrated_score', 'peak_condition', 'context_specificity', 'genetics_autoimmune',
        'lof_bin', 'ot_autoimmune_disease', 'ot_top_disease', 'primary_class', 'suggested_modality',
        'tractability_tier', 'prog_score__Rest', 'prog_score__Stim8hr', 'prog_score__Stim48hr',
        'emp_fdr__Rest', 'emp_fdr__Stim8hr', 'emp_fdr__Stim48hr', 'peak_kd', 'ensembl_id']

# Plot
apply_figure_style()
conds = ['Rest', 'Stim8hr', 'Stim48hr']
clab = ['Rest', 'Stim 8h', 'Stim 48h']
dir_col = {'driver_antagonize': '#b2182b', 'brake_agonize': '#2166ac'}

fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4))
axes = axes.ravel()
for k, g in enumerate(picks):
    ax = axes[k]
    r = df.loc[g]
    ps = [r[f'prog_score__{c}'] for c in conds]
    fdr = [r[f'emp_fdr__{c}'] for c in conds]
    col = dir_col[r['direction']]
    x = np.arange(3)
    ax.axhline(0, color='#bbb', lw=0.8, zorder=1)
    ax.plot(x, ps, '-', color=col, lw=1.8, zorder=2)
    for xi, (p, q) in enumerate(zip(ps, fdr)):
        sig = q < 0.05
        ax.plot(xi, p, 'o', ms=9 if sig else 6, color=col if sig else 'white',
                mec=col, mew=1.5, zorder=3)
        if sig:
            star = '***' if q < 1e-3 else ('**' if q < 1e-2 else '*')
            ax.annotate(star, (xi, p), xytext=(0, 8 if p >= 0 else -14), textcoords='offset points',
                        ha='center', fontsize=8, color=col, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(clab, fontsize=7.5)
    ax.set_xlim(-0.35, 2.35)
    ymax = max(1.2, max(abs(v) for v in ps) * 1.35)
    ax.set_ylim(-ymax, ymax)
    ax.set_title(f"$\\it{{{g}}}$", fontsize=11, loc='left', color=col, pad=14)
    mod = {'small molecule': 'small molecule', 'antibody/biologic (surface)': 'antibody'}.get(r['suggested_modality'], r['suggested_modality'])
    dis = r['ot_autoimmune_disease'] or r['ot_top_disease'] or ''
    if len(dis) > 28: dis = dis[:27] + '…'
    role = 'driver (KD lowers program)' if r['direction'] == 'driver_antagonize' else 'brake (KD raises program)'
    txt = (f"{role}\n"
           f"genetics: {r['genetics_tier']} (AI {r['genetics_autoimmune']:.2f})\n"
           f"disease: {dis}\n"
           f"modality: {mod}\n"
           f"novelty: {r['novelty_class'].replace('_', ' ')}")
    ax.text(0.03, 0.03 if ps[np.argmax(np.abs(ps))] > 0 else 0.97, txt, transform=ax.transAxes,
            fontsize=5.9, va='bottom' if ps[np.argmax(np.abs(ps))] > 0 else 'top', ha='left', color='#333',
            bbox=dict(boxstyle='round,pad=0.35', fc='#f7f7f7', ec='#ddd', lw=0.5))
    if k % 3 == 0: ax.set_ylabel('pro-inflammatory\nprogram score')
    ax.spines[['right', 'top']].set_visible(False)
    panel_letter(ax, chr(97 + k))

h = [ml.Line2D([], [], color='#b2182b', lw=1.8, marker='o', label='antagonize driver (KD ↓ program)'),
     ml.Line2D([], [], color='#2166ac', lw=1.8, marker='o', label='agonize brake (KD ↑ program)'),
     ml.Line2D([], [], color='#888', lw=0, marker='o', mfc='#888', label='filled = FDR<0.05'),
     ml.Line2D([], [], color='#888', lw=0, marker='o', mfc='white', mec='#888', label='open = n.s.')]
fig.legend(handles=h, frameon=False, fontsize=7, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.02))
fig.suptitle('Deep-dive: top novel drug-target nominations — context-specific perturbation trajectories', fontsize=10.5, y=1.0)
fig.text(0.5, 0.945, '* FDR<0.05   ** FDR<0.01   *** FDR<0.001   (empirical null, per condition)', ha='center', fontsize=6.5, color='#666')
fig.tight_layout(rect=[0, 0.03, 1, 0.94])
fig.savefig("top_target_vignettes.png", dpi=200, bbox_inches='tight')
