# Novelty audit & citation-linked evidence dossiers

*Cross-source evidence via the Amass life-science API (PatentCore, TrialCore, DrugCore, BioMedCore) for the 11 top novel target nominations. This layer is **verification and annotation, not a prioritisation axis** — patent/literature volume tracks study bias, so it sits beside the ranking, never inside it.*

**Method.** For each target: (1) patents mentioning the gene vs those explicitly *target-directed* (claiming an inhibitor/modulator/antagonist of the gene), sub-classified by immune indication; (2) interventional trials naming the gene; (3) drugs indexed against the target; (4) top-cited primary literature, and literature linking the target to its disease anchor, with PMIDs/DOIs.

## Summary

| Target | Dir. | Rank | Disease anchor | Patents (mention / target-directed / immune) | Trials | Drugs | Novelty tier |
|---|---|---|---|---|---|---|---|
| **SIK2** | agon | 19 | asthma | 77 / 10 / 4 | 0 | 0 | emerging immune IP |
| **CTSH** | antag | 31 | type 1 diabetes mellitus | 96 / 0 / 0 | 0 | 0 | unencumbered |
| **NEK6** | antag | 45 | Abnormality of the skeletal  | 92 / 2 / 0 | 0 | 0 | patented other indication |
| **DAGLB** | antag | 51 | osteogenesis imperfecta, typ | 100 / 0 / 0 | 0 | 0 | unencumbered |
| **ADAM19** | agon | 57 | coronary artery disorder | 88 / 0 / 0 | 0 | 0 | unencumbered |
| **BRD1** | antag | 70 | synovium disorder | 100 / 0 / 0 | 0 | 0 | unencumbered |
| **PRKX** | antag | 85 | neurodegenerative disease | 94 / 0 / 0 | 0 | 0 | unencumbered |
| **IRAK3** | agon | 91 | mixed connective tissue dise | 87 / 3 / 1 | 1 | 0 | emerging immune IP |
| **PIK3C3** | agon | 119 | neurodegenerative disease | 70 / 0 / 0 | 0 | 0 | unencumbered |
| **NCOA2** | antag | 138 | atrial fibrillation | 100 / 0 / 0 | 0 | 0 | unencumbered |
| **GPR55** | antag | 164 | preeclampsia | 94 / 3 / 2 | 0 | 0 | emerging immune IP |

**Headline:** 7/11 nominations are fully unencumbered — no target-directed patents, no interventional trials, no drugs. All 11 have **zero approved/clinical drugs against the target**. Three (SIK2, IRAK3, GPR55) show emerging target-directed IP, only partly in immune indications — flagged honestly as reduced novelty.

---

## Per-target dossiers

### SIK2 — brake agonize · rank 19 · kinase
*Disease anchor: asthma. Novelty tier: ⚠ Emerging immune IP (early patents/trials in immune space)*

**IP/clinical status:** 77 patents mention SIK2; 10 are target-directed (4 in immune indications); 0 trials; 0 drugs against target.

Key target-directed patents:
- US-2014256704-A1 — Substituted 1H-Pyrrolo [2, 3-b] pyridine and 1H-Pyrazolo [3, 4-b] pyri · ARRIEN PHARMACEUTICALS LLC _(immune)_
- US-2023045929-A1 — Isoquinoline derivatives as sik2 inhibitors · CANCER RESEARCH TECH LTD _(immune)_
- US-2022296583-A1 — Compositions and methods for treatment of ovarian and breast cancer · UNIV TEXAS
- WO-2021084265-A1 — Isoquinoline derivatives as sik2 inhibitors · CANCER RESEARCH TECH LTD _(immune)_
- WO-2026060155-A1 — Sik2 modulators and uses thereof · NIMBUS SALACIA INC
- TW-202547512-A — Sik2 modulators and uses thereof · NIMBUS SALACIA INC

Top-cited primary literature:
- [160 cit] Adipose-specific expression, phosphorylation of Ser794 in insulin recept (2003, *The Journal of biological chemistry*) — https://doi.org/10.1074/jbc.M211770200
- [149 cit] SIK2 is a key regulator for neuronal survival after ischemia via TORC1-C (2011, *Neuron*) — https://doi.org/10.1016/j.neuron.2010.12.004
- [143 cit] SIK2 is a centrosome kinase required for bipolar mitotic spindle formati (2010, *Cancer cell*) — https://doi.org/10.1016/j.ccr.2010.06.018
- [88 cit] SIK2 regulates CRTCs, HDAC4 and glucose uptake in adipocytes (2015, *Journal of cell science*) — https://doi.org/10.1242/jcs.153932

### CTSH — driver antagonize · rank 31 · enzyme
*Disease anchor: type 1 diabetes mellitus. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 96 patents mention CTSH; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [105 cit] CTSH regulates β-cell function and disease progression in newly diagnose (2014, *Proceedings of the National Academy of Sciences of the United States of America*) — https://doi.org/10.1073/pnas.1402571111
- [57 cit] Cathepsin H regulated by the thyroid hormone receptors associate with tu (2011, *Oncogene*) — https://doi.org/10.1038/onc.2010.585
- [33 cit] The human cathepsin H gene encodes two novel minor histocompatibility an (2006, *British journal of haematology*) — https://doi.org/10.1111/j.1365-2141.2006.06205.x
- [31 cit] Effects of chitin and sepia ink hybrid hemostatic sponge on the blood pa (2014, *Marine drugs*) — https://doi.org/10.3390/md12042269

Literature linking CTSH to type 1 diabetes mellitus:
- Integrating multi-omics data to analyze the potential pathogenic mechani (2024) — https://doi.org/10.1093/bfgp/elad052
- Genetic and environmental factors regulate the type 1 diabetes gene CTSH (2021) — https://doi.org/10.1016/j.jbc.2021.100774
- CTSH regulates β-cell function and disease progression in newly diagnose (2014) — https://doi.org/10.1073/pnas.1402571111

### NEK6 — driver antagonize · rank 45 · kinase
*Disease anchor: Abnormality of the skeletal system. Novelty tier: ○ Patented for other indication (no immune IP)*

**IP/clinical status:** 92 patents mention NEK6; 2 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Key target-directed patents:
- US-2021147395-A1 — Nek6 kinase inhibitors useful for the treatment of solid tumors · CONSIGLIO NAZIONALE RICERCHE, MOLIPHARMA S R L
- CN-118236502-A — Application of NEK6 inhibitor combined with CDK4/6 inhibitor in prepar · RUIJIN HOSPITAL SHANGHAI JIAOTONG UNIV SCHOOL MEDICINE

Top-cited primary literature:
- [51 cit] Role of NEK6 in tumor promoter-induced transformation in JB6 C141 mouse  (2010, *The Journal of biological chemistry*) — https://doi.org/10.1074/jbc.M110.137190
- [49 cit] NIMA-related kinases 6, 4, and 5 interact with each other to regulate mi (2011, *The Plant journal : for cell and molecular biology*) — https://doi.org/10.1111/j.1365-313X.2011.04652.x
- [40 cit] NIMA-related kinase NEK6 affects plant growth and stress response in Ara (2011, *The Plant journal : for cell and molecular biology*) — https://doi.org/10.1111/j.1365-313X.2011.04733.x
- [40 cit] Overexpression of NIMA-related kinase 6 (NEK6) contributes to malignant  (2018, *Pathology, research and practice*) — https://doi.org/10.1016/j.prp.2018.07.030

### DAGLB — driver antagonize · rank 51 · enzyme
*Disease anchor: osteogenesis imperfecta, type 21. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 100 patents mention DAGLB; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [212 cit] DAGLβ inhibition perturbs a lipid network involved in macrophage inflamm (2012, *Nature chemical biology*) — https://doi.org/10.1038/nchembio.1105
- [64 cit] Deficiency in endocannabinoid synthase DAGLB contributes to early onset  (2022, *Nature communications*) — https://doi.org/10.1038/s41467-022-31168-9
- [54 cit] AP-4-mediated axonal transport controls endocannabinoid production in ne (2022, *Nature communications*) — https://doi.org/10.1038/s41467-022-28609-w
- [21 cit] Case of 7p22.1 Microduplication Detected by Whole Genome Microarray (REV (2015, *Case reports in genetics*) — https://doi.org/10.1155/2015/212436

Literature linking DAGLB to osteogenesis imperfecta, type :
- Modern classification and molecular-genetic aspects of osteogenesis impe (2020) — https://doi.org/10.18699/VJ20.614
- The Spine in Patients With Osteogenesis Imperfecta (2017) — https://doi.org/10.5435/JAAOS-D-15-00169
- Potential of gene therapy for treating osteogenesis imperfecta (2000) — https://doi.org/10.1097/00003086-200010001-00017

### ADAM19 — brake agonize · rank 57 · enzyme
*Disease anchor: coronary artery disorder. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 88 patents mention ADAM19; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [136 cit] Essential role for ADAM19 in cardiovascular morphogenesis (2004, *Molecular and cellular biology*) — https://doi.org/10.1128/MCB.24.1.96-104.2004
- [126 cit] Catalytic properties of ADAM19 (2003, *The Journal of biological chemistry*) — https://doi.org/10.1074/jbc.M302781200
- [60 cit] ADAM19/adamalysin 19 structure, function, and role as a putative target  (2009, *Current pharmaceutical design*) — https://doi.org/10.2174/138161209788682352
- [45 cit] Xenopus ADAM19 is involved in neural, neural crest and muscle developmen (2008, *Mechanisms of development*) — https://doi.org/10.1016/j.mod.2008.10.010

Literature linking ADAM19 to coronary artery disorder:
- A rare case of congenital atresia of the left main coronary artery (2022) — https://doi.org/10.1016/j.radcr.2021.10.012
- Coronary artery ectasia in a patient with myocardial infarction (2011) — https://doi.org/10.5830/CVJA-2010-033
- Spontaneous coronary artery dissection in a patient with autosomal domin (2016) — https://doi.org/10.1186/s13256-016-0832-8

### BRD1 — driver antagonize · rank 70 · enzyme
*Disease anchor: synovium disorder. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 100 patents mention BRD1; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [290 cit] Isolation and characterization of a rice dwarf mutant with a defect in b (2002, *Plant physiology*) — https://doi.org/10.1104/pp.007179
- [77 cit] Evidence implicating BRD1 with brain development and susceptibility to b (2006, *Molecular psychiatry*) — https://doi.org/10.1038/sj.mp.4001885
- [54 cit] Expression of BrD1, a plant defensin from Brassica rapa, confers resista (2009, *Molecules and cells*) — https://doi.org/10.1007/s10059-009-0117-9
- [52 cit] Support of association between BRD1 and both schizophrenia and bipolar a (2010, *American journal of medical genetics. Part B, Neuropsychiatric genetics : the official publication of the International Society of Psychiatric Genetics*) — https://doi.org/10.1002/ajmg.b.31023

Literature linking BRD1 to synovium disorder:
- CD4+ tissue-resident memory Th17 cells are a major source of IL-17A in S (2025) — https://doi.org/10.1016/j.ard.2025.04.018
- Pathology of the synovium (2000) — https://doi.org/10.1309/LWW3-5XK0-FKG9-HDRK

### PRKX — driver antagonize · rank 85 · kinase
*Disease anchor: neurodegenerative disease. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 94 patents mention PRKX; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [143 cit] Abnormal XY interchange between a novel isolated protein kinase gene, PR (1997, *Human molecular genetics*) — https://doi.org/10.1093/hmg/6.11.1985
- [94 cit] PrKX is a novel catalytic subunit of the cAMP-dependent protein kinase r (1999, *The Journal of biological chemistry*) — https://doi.org/10.1074/jbc.274.9.5370
- [64 cit] Adeno-associated virus Rep78 protein interacts with protein kinase A and (1998, *Journal of virology*) — https://doi.org/10.1128/JVI.72.10.7916-7925.1998
- [58 cit] Transforming growth factor (TGF)-β-activated kinase 1 (TAK1) activation  (2014, *The Journal of biological chemistry*) — https://doi.org/10.1074/jbc.M114.559963

Literature linking PRKX to neurodegenerative disease:
- Neurodegenerative disease (2007) — https://doi.org/10.1007/s11307-007-0099-y

### IRAK3 — brake agonize · rank 91 · kinase
*Disease anchor: mixed connective tissue disease. Novelty tier: ⚠ Emerging immune IP (early patents/trials in immune space)*

**IP/clinical status:** 87 patents mention IRAK3; 3 are target-directed (1 in immune indications); 1 trials; 0 drugs against target.

Key target-directed patents:
- US-2026028343-A1 — Substituted Imidazopyrazine Compounds as Ligand Directed Degraders of  · CELGENE CORP
- EP-4638454-A1 — Substituted imidazo-based compounds as ligand directed degraders of ir · CELGENE CORP
- US-2023374606-A1 — Methods for treating melanoma · UNIV LOUISVILLE RES FOUND INC _(immune)_

Top-cited primary literature:
- [148 cit] Inhibition of interleukin-1 receptor-associated kinase 1 (IRAK1) as a th (2018, *Oncotarget*) — https://doi.org/10.18632/oncotarget.26058
- [57 cit] IRAK3 modulates downstream innate immune signalling through its guanylat (2019, *Scientific reports*) — https://doi.org/10.1038/s41598-019-51913-3
- [34 cit] Dimeric Structure of the Pseudokinase IRAK3 Suggests an Allosteric Mecha (2021, *Structure (London, England : 1993)*) — https://doi.org/10.1016/j.str.2020.11.004
- [28 cit] A systematic review and meta-analyses of interleukin-1 receptor associat (2022, *PloS one*) — https://doi.org/10.1371/journal.pone.0263968

Literature linking IRAK3 to mixed connective tissue diseas:
- Mixed connective-tissue disease (1977) — https://doi.org/10.1056/NEJM197702242960812
- Mixed connective tissue disease (1993) — https://doi.org/10.1093/rheumatology/32.7.645
- Mixed connective tissue disease (1977) — https://pubmed.ncbi.nlm.nih.gov/921429

### PIK3C3 — brake agonize · rank 119 · enzyme
*Disease anchor: neurodegenerative disease. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 70 patents mention PIK3C3; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [92 cit] The mammalian class 3 PI3K (PIK3C3) is required for early embryogenesis  (2011, *PloS one*) — https://doi.org/10.1371/journal.pone.0016358
- [65 cit] Autophagy-related protein PIK3C3/VPS34 controls T cell metabolism and fu (2021, *Autophagy*) — https://doi.org/10.1080/15548627.2020.1752979
- [62 cit] Deficiency in class III PI3-kinase confers postnatal lethality with IBD- (2018, *Nature communications*) — https://doi.org/10.1038/s41467-018-05105-8
- [49 cit] Pik3c3 deletion in pyramidal neurons results in loss of synapses, extens (2011, *Neuroscience*) — https://doi.org/10.1016/j.neuroscience.2010.10.035

Literature linking PIK3C3 to neurodegenerative disease:
- Conditional knockout of pik3c3 causes a murine muscular dystrophy (2014) — https://doi.org/10.1016/j.ajpath.2014.02.012
- Development of in vitro PIK3C3/VPS34 complex protein assay for autophagy (2015) — https://doi.org/10.1016/j.ab.2015.04.004

### NCOA2 — driver antagonize · rank 138 · enzyme
*Disease anchor: atrial fibrillation. Novelty tier: ✓ Unencumbered (no target-directed patents, trials, or drugs)*

**IP/clinical status:** 100 patents mention NCOA2; 0 are target-directed (0 in immune indications); 0 trials; 0 drugs against target.

Top-cited primary literature:
- [149 cit] Fusion of the AHRR and NCOA2 genes through a recurrent translocation t(5 (2012, *Genes, chromosomes & cancer*) — https://doi.org/10.1002/gcc.21939
- [97 cit] Androgen deprivation-induced NCoA2 promotes metastatic and castration-re (2014, *The Journal of clinical investigation*) — https://doi.org/10.1172/JCI76412
- [51 cit] Soft tissue angiofibroma: Clinicopathologic, immunohistochemical and mol (2017, *Genes, chromosomes & cancer*) — https://doi.org/10.1002/gcc.22478
- [43 cit] Disruption of NCOA2 by recurrent fusion with LACTB2 in colorectal cancer (2016, *Oncogene*) — https://doi.org/10.1038/onc.2015.72

Literature linking NCOA2 to atrial fibrillation:
- Atrial fibrillation as a marker of occult cancer (2014) — https://doi.org/10.1371/journal.pone.0102861
- A Discussion of the Contemporary Prediction Models for Atrial Fibrillati (2023) — https://doi.org/10.18103/mra.v11i10.4481
- Atrial fibrillation among patients under investigation for suspected obs (2017) — https://doi.org/10.1371/journal.pone.0171575

### GPR55 — driver antagonize · rank 164 · GPCR
*Disease anchor: preeclampsia. Novelty tier: ⚠ Emerging immune IP (early patents/trials in immune space)*

**IP/clinical status:** 94 patents mention GPR55; 3 are target-directed (2 in immune indications); 0 trials; 0 drugs against target.

Key target-directed patents:
- US-2018271808-A1 — Methods of reducing chemoresistance and treating cancer · MITCHELL WOODS PHARMACEUTICALS INC
- CN-120114599-A — Application of GPR55 inhibitor in preparation of product for enhancing · BIG HEALTH RES INSTITUTE OF HEFEI INTEGRATED SCIENCES CENTER, UNIV SCIENCE & TECHNOLOGY CHINA _(immune)_
- CN-119235846-A — Application of GPR55 antagonist CID16020046 in treatment of diabetic n · UNIV SHANXI MEDICAL _(immune)_

Top-cited primary literature:
- [695 cit] GPR55 is a cannabinoid receptor that increases intracellular calcium and (2008, *Proceedings of the National Academy of Sciences of the United States of America*) — https://doi.org/10.1073/pnas.0711278105
- [236 cit] GPR55: a new member of the cannabinoid receptor clan? (2007, *British journal of pharmacology*) — https://doi.org/10.1038/sj.bjp.0707464
- [176 cit] The orphan G protein-coupled receptor GPR55 promotes cancer cell prolife (2011, *Oncogene*) — https://doi.org/10.1038/onc.2010.402
- [111 cit] GPR55, a G-protein coupled receptor for lysophosphatidylinositol, plays  (2013, *PloS one*) — https://doi.org/10.1371/journal.pone.0060314
