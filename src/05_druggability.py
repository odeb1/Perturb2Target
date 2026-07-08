#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 5 — Druggability/tractability: protein-class -> suggested modality (priority-ordered),
tractability tier. -> candidates_step5.parquet, druggability_landscape.png.

Conda environment: python
Produces artifact: candidates_step5.parquet
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import pandas as pd
import numpy as np

dcm = pd.read_csv('inputs/druggable_class_map.csv')
cand = pd.read_parquet('inputs/candidates_step4.parquet')
ann = pd.read_parquet('inputs/ot_annotation.parquet')

full = cand.join(ann, how='left')
full['genetics_autoimmune'] = full[['autoimmune_ot_score','ot_autoimmune_score']].max(axis=1)
full['lof_constrained'] = full['lof_bin'] <= 1

def gtier(r):
    s = r['genetics_autoimmune']
    if s >= 0.5: return 'strong'
    if s >= 0.25: return 'moderate'
    if s > 0.1: return 'weak'
    return 'none'
full['genetics_tier'] = full.apply(gtier, axis=1)

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

full2.to_parquet('candidates_step5.parquet')
