#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 6 — Clinical novelty from OT known-drug/clinical-candidate max stage. -> candidates_step6.parquet,
clinical_novelty.csv, novelty_breakdown.png.

Conda environment: python
Produces artifact: candidates_step6.parquet
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import pandas as pd
import numpy as np

# Load inputs
dcm = pd.read_csv('inputs/druggable_class_map.csv')
cand = pd.read_parquet('inputs/candidates_step4.parquet')
obs = pd.read_parquet('inputs/de_obs.parquet')
ann = pd.read_parquet('inputs/ot_annotation.parquet')

# Reconstruct candidate mask
cand['kd_achieved'] = cand['peak_kd'] < -1.0
cand['program_hit'] = (cand['peak_emp_fdr'] < 0.05) & (cand['peak_abs_score'].abs() > 0.2)
cand['abs_emp_z'] = cand['peak_emp_z'].abs()

# Merge OT annotation onto candidates
full = cand.join(ann, how='left')
full['genetics_autoimmune'] = full[['autoimmune_ot_score', 'ot_autoimmune_score']].max(axis=1)
full['lof_constrained'] = full['lof_bin'] <= 1

def gtier(r):
    s = r['genetics_autoimmune']
    if s >= 0.5: return 'strong'
    if s >= 0.25: return 'moderate'
    if s > 0.1: return 'weak'
    return 'none'

full['genetics_tier'] = full.apply(gtier, axis=1)

# Merge druggable class map
full2 = full.reset_index().merge(dcm, left_on='target_gene', right_on=dcm.columns[0], how='left').set_index('target_gene')

def tb(v):
    return bool(v) if (v is True or v is False) else (False if pd.isna(v) else bool(v))

def num(v):
    return 0 if pd.isna(v) else float(v)

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

full2['tractability_tier'] = full2.apply(tract_tier, axis=1)

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

full2['suggested_modality'] = full2.apply(modality3, axis=1)
full2['is_druggable'] = ~full2['tractability_tier'].isin(['difficult_tf', 'unknown'])

stage_rank = {'APPROVAL': 4, 'PHASE_3': 3, 'PHASE_2_3': 3, 'PHASE_2': 2, 'PHASE_1_2': 2, 'PHASE_1': 1, 'UNKNOWN': 0}

def max_stage(s):
    if not isinstance(s, str) or not s: return np.nan
    vals = [stage_rank.get(t.strip()) for t in s.split(';') if t.strip() in stage_rank]
    return max(vals) if vals else np.nan

full2['max_drug_stage'] = full2['drug_stages'].apply(max_stage)

def novelty(r):
    ns = r.get('n_drugs')
    ms = r.get('max_drug_stage')
    if pd.isna(ns):
        return 'not_assessed'
    if ns == 0: return 'novel_undrugged'
    if ms is not None and not pd.isna(ms) and ms >= 4: return 'approved_drug_exists'
    if ms is not None and not pd.isna(ms) and ms >= 1: return 'clinical_stage_drug'
    return 'preclinical_or_unknown'

full2['novelty_class'] = full2.apply(novelty, axis=1)

paper_noms = {'STAT3', 'IRF4', 'BATF', 'LRRC25', 'FAM20B'}
full2['in_paper_noms'] = full2.index.isin(paper_noms)

full2.to_parquet('candidates_step6.parquet')
