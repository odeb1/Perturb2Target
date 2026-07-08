#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 4 — Human-genetics support: Open Targets GraphQL (autoimmune assoc, LoF constraint,
tractability, known drugs) + gnomAD; genetics tier. -> ot_annotation.parquet, genetics_support.png.

Conda environment: python
Produces artifact: ot_annotation.parquet
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import json, pandas as pd, numpy as np

out = json.load(open("handoff/ot_deep.json"))

AUTOIMMUNE_TA = {'EFO_0005140'}
autoimmune_kw = ['autoimmune','rheumatoid','lupus','psoriasis','psoriatic','inflammatory bowel','crohn',
                 'ulcerative colitis','multiple sclerosis','type 1 diabetes','type i diabetes','celiac',
                 'ankylosing','vitiligo','graves','hashimoto','sjogren','asthma','atopic','allergic',
                 'juvenile idiopathic','spondyloarthropathy','alopecia','immune']
rows=[]
for sym,t in out.items():
    if not t or 'error' in t or 'errors' in t: 
        rows.append({'target_gene':sym}); continue
    gc={c['constraintType']:c for c in (t.get('geneticConstraint') or [])}
    lof=gc.get('lof',{})
    ad=(t.get('associatedDiseases') or {})
    adr=ad.get('rows',[])
    best_ai=0.0; best_ai_name=None; top_name=None; top_score=0.0
    for r in adr:
        dis=r.get('disease',{}); nm=(dis.get('name') or '').lower(); sc=r.get('score',0)
        tas=[ (ta.get('name') or '').lower() for ta in (dis.get('therapeuticAreas') or []) ]
        if sc>top_score: top_score=sc; top_name=dis.get('name')
        is_ai = any(k in nm for k in autoimmune_kw) or any('immune' in ta or 'inflammat' in ta for ta in tas)
        if is_ai and sc>best_ai: best_ai=sc; best_ai_name=dis.get('name')
    tr=t.get('tractability') or []
    def has(mod,label): 
        return any(x['modality']==mod and x['label']==label and x['value'] for x in tr)
    sm_clin = has('SM','Approved Drug') or has('SM','Advanced Clinical') or has('SM','Phase 1 Clinical')
    sm_lig = has('SM','High-Quality Ligand') or has('SM','Structure with Ligand')
    sm_pocket = has('SM','High-Quality Pocket') or has('SM','Med-Quality Pocket')
    sm_family = has('SM','Druggable Family')
    ab_clin = has('AB','Approved Drug') or has('AB','Advanced Clinical') or has('AB','Phase 1 Clinical')
    ab_surface = has('AB','UniProt loc high conf') or has('AB','GO CC high conf') or has('AB','UniProt SigP or TMHMM')
    dc=t.get('drugAndClinicalCandidates') or {}
    ndrug=dc.get('count',0)
    stages=[ (row.get('drug') or {}).get('maximumClinicalStage') for row in dc.get('rows',[]) ]
    rows.append({'target_gene':sym,'lof_score':lof.get('score'),'lof_bin':lof.get('upperBin'),
                 'mis_z':gc.get('mis',{}).get('score'),
                 'ot_top_disease':top_name,'ot_top_score':round(top_score,3),
                 'ot_autoimmune_disease':best_ai_name,'ot_autoimmune_score':round(best_ai,3),
                 'ot_n_assoc_disease':ad.get('count'),
                 'sm_clinical':sm_clin,'sm_ligand':sm_lig,'sm_pocket':sm_pocket,'sm_druggable_family':sm_family,
                 'ab_clinical':ab_clin,'ab_surface':ab_surface,
                 'n_drugs':ndrug,'drug_stages':';'.join(sorted({s for s in stages if s}))})
ann=pd.DataFrame(rows).set_index('target_gene')
ann.to_parquet("ot_annotation.parquet")
