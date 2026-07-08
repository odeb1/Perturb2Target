#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 7 — Integrate 4 evidence components (causal .34 / genetics .30 / drug .22 / novelty .14)
into integrated_score; novel-priority flag. -> target_shortlist.csv, candidates_step7.parquet.

Conda environment: python
Produces artifact: candidates_step7.parquet
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import pandas as pd, numpy as np

full2 = pd.read_parquet("inputs/candidates_step6.parquet")

df = full2.copy()

# ---- 1. CAUSAL component: strength of directional program perturbation (empirical z at peak) ----
# rank-normalize abs_emp_z to [0,1]
df['causal_component'] = (df['abs_emp_z'].rank(pct=True)).clip(0, 1)

# ---- 2. GENETICS component: autoimmune OT association + LoF-constraint bonus ----
g = df['genetics_autoimmune'].fillna(0.0).clip(0, 1)
lof_bonus = df['lof_constrained'].fillna(False).astype(float) * 0.15
df['genetics_component'] = (g + lof_bonus).clip(0, 1)

# ---- 3. DRUGGABILITY component: tractability tier ----
trac_map = {'clinical_precedent': 1.0, 'discovery_precedent': 0.75, 'predicted_tractable': 0.55, 'difficult': 0.25}
df['drug_component'] = df['tractability_tier'].map(trac_map).fillna(0.15)
# hard-TF targets (no modality) get floored low
df.loc[df['suggested_modality'].str.contains('hard', na=False), 'drug_component'] = df['drug_component'].clip(upper=0.30)

# ---- 4. NOVELTY component: undrugged>preclinical>clinical>approved ----
nov_map = {'novel_undrugged': 1.0, 'preclinical_or_unknown': 0.7, 'not_assessed': 0.5, 'clinical_stage_drug': 0.35, 'approved_drug_exists': 0.15}
df['novelty_component'] = df['novelty_class'].map(nov_map).fillna(0.5)

# ---- Composite: weighted geometric-leaning sum. Causal & genetics are the evidence spine. ----
W = dict(causal=0.34, genetics=0.30, drug=0.22, novelty=0.14)
df['integrated_score'] = (W['causal'] * df['causal_component'] + W['genetics'] * df['genetics_component']
                          + W['drug'] * df['drug_component'] + W['novelty'] * df['novelty_component'])

# require: KD achieved & program hit already true for all 1923; add evidence gate for shortlist ranking
df['evidence_n'] = (df['program_hit'].astype(int) + (df['genetics_tier'].isin(['strong', 'moderate'])).astype(int)
                    + df['is_druggable'].fillna(False).astype(int) + (df['novelty_class'] == 'novel_undrugged').astype(int))
df = df.sort_values('integrated_score', ascending=False)
df['rank'] = np.arange(1, len(df) + 1)

# Novel-priority flag: the differentiated nominations (novel/undrugged, decent genetics, druggable)
df['novel_priority'] = ((df['novelty_class'].isin(['novel_undrugged', 'preclinical_or_unknown']))
                        & (df['genetics_tier'].isin(['strong', 'moderate']))
                        & (df['is_druggable'].fillna(False)))

df.to_parquet("candidates_step7.parquet")
print("shortlist rows:", len(df), "| novel_priority:", int(df['novel_priority'].sum()))
