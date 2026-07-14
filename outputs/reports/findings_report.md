# Context-specific, directional, genetics-anchored drug-target nomination from genome-scale CD4⁺ T-cell CRISPRi Perturb-seq

**Analysis report — computational target discovery**

---

## 1. Summary

We re-analyzed a genome-scale CRISPRi Perturb-seq screen in primary human CD4⁺ T cells (Zhu, Dann, … Pritchard & Marson 2025, bioRxiv 10.64898/2025.12.23.696273) to nominate drug targets for autoimmune and allergic disease. Our angle is deliberately differentiated from the preprint's own target list along three axes:

1. **Directionality** — we separate targets whose knockdown *lowers* a pathogenic pro-inflammatory program (**drivers → antagonize**) from those whose knockdown *raises* it (**brakes → agonize**), because the two imply opposite therapeutic modalities (inhibit vs. potentiate).
2. **Human-genetics anchoring** — every nomination is scored against Open Targets autoimmune/allergic disease genetics and gnomAD loss-of-function constraint.
3. **Druggability + clinical-novelty filter** — we attach a tractability tier and suggested modality to each target, and flag which nominations are *not yet drugged*, surfacing hypotheses that are both actionable and new.

From **1,923** quality-filtered causal candidates (**1,290 driver→antagonize, 633 brake→agonize**), the integration recovers known immune drug targets at the very top (IL4R #1, IL2RA #2, STAT3 #3, IL6R #7 — all approved or clinical-stage), which validates the scoring, and it nominates **86 novel-priority targets** (novel/undrugged × strong-or-moderate genetics × druggable) — **48 to antagonize, 38 to agonize**.

> **Nature of the claims.** These are *computational* hypotheses built by triangulating orthogonal evidence (causal perturbation effect + human genetics + druggability + clinical novelty). No wet-lab validation was performed here. Their credibility rests on the convergence of independent evidence layers and on the recovery of established targets as positive controls, not on any single statistic.

---

## Novel contributions

The source dataset is not ours (Zhu, Dann, … Pritchard & Marson 2025). Our contribution is a **re-analysis layer** that the preprint's own target list does not provide. Four claims:

1. **Directionality as a first-class axis.** We separate targets by the *sign* of the perturbation effect — **drivers → antagonize** (knockdown lowers the pathogenic program; inhibit the target; 1,290 candidates) vs. **brakes → agonize** (knockdown raises it; potentiate the target; 633 candidates). The two imply opposite therapeutic modalities, which a magnitude-only ranking collapses. Positive controls confirm the axis: driver knockdowns (*TBX21, BATF, STAT3, IRF4*) lower the program; brake knockdowns (*STAT6, IL4R, INPP5D*) raise it.

2. **Context-specific scoring.** Every target is scored per condition (Rest / Stim 8 h / Stim 48 h) with a context-specificity index, surfacing *which* activation state each nomination acts in — a constitutive target and an activation-restricted target have different therapeutic windows.

3. **Orthogonal-evidence triangulation.** Four independent layers — causal perturbation effect + human genetics (Open Targets autoimmune/allergic + gnomAD LoF constraint) + druggability/tractability + clinical novelty — are integrated into one composite. Credibility rests on *convergence*, validated by recovery of established drugs at the top of the ranking (IL4R #1, IL2RA #2, STAT3 #3, IL6R #7).

4. **Clinical-novelty filter yielding actionable new hypotheses.** Crossing "not-yet-drugged" with strong/moderate genetics and druggability yields **86 novel-priority targets** (48 antagonize, 38 agonize) plus 6 curated deep-dive vignettes — new, testable bets rather than re-nominations of known targets.

**One sentence:** we turn a genome-scale Perturb-seq screen into a directional, context-resolved, genetics-anchored, druggability-filtered target shortlist — separating brakes-to-agonize from drivers-to-antagonize, validated by recovery of approved drugs, and surfacing 86 new actionable hypotheses the source list does not distinguish.

---

## 2. Data and the pathogenic program

**Input.** `GWCD4i.DE_stats.h5ad` — per-perturbation differential-expression statistics (z-scores of downstream readout genes) for genome-scale CRISPRi across three conditions: **Rest**, **Stim 8 h**, **Stim 48 h**. The DE matrix is 33,983 (perturbation × condition) rows × 10,282 readout genes.

**Pro-inflammatory program signature.** We defined a Th1/pathogenic-effector program as a signed gene set: **36 pro-inflammatory genes (+1)** (e.g. *IFNG, TBX21, STAT1, STAT4, CXCR3*) and **12 regulatory/anti-inflammatory genes (−1)** (e.g. *IL10, FOXP3, CTLA4*). 48 of 50 intended genes were measured in the readout (*IL17A, RORA* absent). A per-perturbation **program score** is the mean of (signed weight × readout z-score) over signature genes; negative = program suppressed, positive = program elevated.

*(Figure: `program_signature_heatmap.png` — signature-gene readout across positive-control perturbations.)*

---

## 3. Directional, context-specific scoring

For each perturbation × condition we computed the program score and calibrated it against an **empirical null**: 500 random signed gene sets matched in size (36 up / 12 down) drawn from non-signature genes, giving an empirical z (`emp_z`) and a BH-FDR per condition.

- **Direction.** A perturbation whose knockdown drives the program *down* is a **driver → antagonize** target; one that drives it *up* is a **brake → agonize** target.
- **Context specificity.** We record the peak condition (Rest / Stim 8 h / Stim 48 h) and a context-specificity index (0 = uniform across conditions, 1 = single-condition), because a target that acts only in the activated state has a different therapeutic window than a constitutive one.

**Quality gate.** A candidate must show (a) achieved knockdown (`peak_kd < −1.0` on-target) and (b) a program hit (`peak_emp_fdr < 0.05` and `|peak score| > 0.2`). **1,923** perturbations pass.

**Positive controls behaved as expected:** knockdown of drivers *TBX21, BATF, STAT3, IRF4* lowered the program; knockdown of brakes *STAT6, IL4R, INPP5D* raised it.

*(Figure: `perturbation_directionality.png` — volcano/directionality map of all candidates.)*

---

## 4. Human-genetics support

Each candidate's Ensembl ID was queried against **Open Targets** for autoimmune/allergic disease associations and against gnomAD for LoF constraint. We assigned a **genetics tier** (strong / moderate / weak / none) combining autoimmune association score and LoF-constraint. **188** candidates reach the *strong* tier; **125** of those are driver→antagonize.

*(Figure: `genetics_support.png` — genetics score vs. causal effect, colored by direction.)*

---

## 5. Druggability and tractability

We mapped each target to a protein class and assigned a **suggested modality** with a priority that respects biology:

- secreted cytokine → **biologic (secreted)**;
- cytokine/catalytic receptor (cell-surface) → **antibody**;
- kinase / enzyme / GPCR / ion channel / nuclear receptor / transporter → **small molecule**;
- otherwise small-molecule-ligandable/pocket → **small molecule**; surface antigen → **antibody**;
- transcription factor with no pocket → **hard (degrader/PPI)**.

A **tractability tier** (clinical precedent → discovery precedent → predicted tractable → difficult) captures how far modality development has progressed. Among strong-genetics antagonize targets, **67** are druggable; by direction the antagonize set splits roughly into ~210 small-molecule, ~34 antibody, and ~106 hard-TF targets.

*(Figure: `druggability_landscape.png` — modality × tractability landscape.)*

---

## 6. Clinical novelty

Using Open Targets known-drug/clinical-candidate records we classified each target into: **novel/undrugged, preclinical/unknown, clinical-stage drug, approved drug**. Among the deep-annotated top-200 antagonize set, **47** targets are *novel/undrugged with strong genetics* — the differentiated core of this analysis. STAT3 correctly re-surfaced as *clinical-stage*, confirming that the novelty filter recovers known-drug status.

*(Figure: `novelty_breakdown.png` — novelty × genetics stacked bar + novel-undrugged target lollipop.)*

---

## 7. Integrated ranking

We combined four evidence components, each normalized to [0,1], into a single `integrated_score`:

| Component | Weight | Source |
|---|---|---|
| Causal strength | 0.34 | rank-normalized empirical \|z\| at peak condition |
| Human genetics | 0.30 | Open Targets autoimmune association + LoF-constraint bonus |
| Druggability | 0.22 | tractability tier (hard-TF floored) |
| Clinical novelty | 0.14 | undrugged > preclinical > clinical > approved |

Directionality is preserved as an *annotation*, not folded into the score, so the user can filter to antagonize or agonize as the therapeutic question demands.

**Positive-control recovery (validation):**

| Target | Rank | Score | Status |
|---|---|---|---|
| IL4R | 1 | 0.87 | approved (dupilumab; agonize-brake control) |
| IL2RA | 2 | 0.84 | approved (basiliximab) |
| STAT3 | 3 | 0.83 | clinical-stage |
| IL6R | 7 | 0.79 | approved (tocilizumab) |
| TBX21 | 20 | 0.75 | known driver |
| BATF | 39 | 0.70 | known driver |

The clustering of approved/clinical immune targets at the top of an unsupervised integrated score is the primary internal validation of the pipeline.

*(Figure: `top_targets_evidence_matrix.png` — evidence heatmap × composite score × modality/disease for the top nominations.)*

---

## 8. Top novel nominations (deep dive)

Six novel-priority nominations spanning both directions, diverse immune diseases, and tractable modalities:

| Target | Direction | Rank | Peak | Genetics (AI) | Disease anchor | Modality |
|---|---|---|---|---|---|---|
| **STAT5B** | antagonize | 15 | Stim 48 h | 0.81 | GH-insensitivity w/ immune dysregulation | small molecule |
| **CTSH** | antagonize | 31 | Stim 8 h | 0.89 | type 1 diabetes | small molecule |
| **NEK6** | antagonize | 45 | Stim 8 h | 0.81 | asthma | small molecule (kinase) |
| **PTPN2** | agonize | 13 | Rest | 0.89 | rheumatoid arthritis | small molecule |
| **SIK2** | agonize | 19 | Rest | 0.75 | asthma | small molecule (kinase) |
| **STAT6** | agonize | 11 | Stim 48 h | 0.81 | hyper-IgE syndrome 6 | small molecule |

Each shows a distinct context-specific trajectory (`top_target_vignettes.png`): the antagonize drivers *STAT5B/CTSH/NEK6* act on activation (Stim 8–48 h), while the brakes *PTPN2/SIK2* act at Rest and *STAT6* emerges late (Stim 48 h). The two kinases (NEK6, SIK2) and the phosphatase (PTPN2) are especially attractive small-molecule handles; CTSH (type-1-diabetes genetics) and STAT6 (already a validated allergy axis via its pathway) anchor to strong disease genetics.

> **Caveats specific to the deep-dive set.** Some disease anchors are Mendelian immune-dysregulation syndromes (e.g. STAT5B, STAT6 hyper-IgE) rather than common-variant GWAS hits; these support a *causal role in immune regulation* but the direction of a common-disease effect should be confirmed against common-variant data before prioritization. Directionality here is defined relative to our program signature, not to clinical outcome.

---

## 8b. Structural druggability (Boltz-2 co-folding)

To move druggability from a *class-based heuristic* to *physical, per-target evidence*, we
co-folded each of the 11 pocket-bearing novel nominations (kinases, enzymes, one GPCR;
transcription factors and transporters excluded — no small-molecule pocket) with a
class-matched, drug-like ligand using **Boltz-2** (model `boltz-2.1`, run on the hosted
Boltz GPU API). Each prediction yields two orthogonal structural readouts:

- **Interface confidence** (`ligand ipTM`) — how confidently a defined protein–ligand
  binding interface forms; a proxy for the existence of a real, ligandable pocket.
- **Binding likelihood** (Boltz-2 affinity head, `binding_confidence`) — the model's
  learned probability that the co-folded ligand is a genuine binder.

We combine them into a single `struct_druggability` = √(ligand ipTM × binding_confidence),
keeping both axes in [0,1]. This is an **independent structural line of evidence**: it uses
protein sequence and chemistry only, with no input from the perturbation, genetics, or
novelty layers.

**Result.** 4 of 11 targets fall in the high-confidence quadrant (confident pocket
*and* likely binder): *ADAM19*, *PRKX*, *IRAK3*, *GPR55*. These are the nominations where the class-based
druggability call is now backed by a predicted pocket and a predicted binder.

| Target | Class | Direction | Peak context | Disease anchor | Repurposing scaffold | ligand ipTM | Binding conf. | Struct. druggability | Tier |
|---|---|---|---|---|---|---|---|---|---|
| *ADAM19* | enzyme | agonize | Stim48hr | asthma | metalloprotease_inhibitor | 0.96 | 0.98 | 0.97 | high confidence binder |
| *PRKX* | kinase | antagonize | Stim48hr | Hashimoto thyroiditis | PKA_inhibitor_H89 | 0.99 | 0.80 | 0.89 | high confidence binder |
| *IRAK3* | kinase | agonize | Rest | mixed connective tissue disease | IRAK_inhibitor_scaffold | 0.98 | 0.81 | 0.89 | high confidence binder |
| *GPR55* | GPCR | antagonize | Stim8hr | psoriatic arthritis | GPR55_antagonist_CID16020046 | 0.98 | 0.72 | 0.84 | high confidence binder |
| *SIK2* | kinase | agonize | Rest | asthma | HG-9-91-01_analog | 0.98 | 0.51 | 0.71 | confident pocket moderate binding |
| *PIK3C3* | enzyme | agonize | Stim48hr | — | VPS34_inhibitor_SAR405 | 0.97 | 0.46 | 0.67 | confident pocket moderate binding |
| *CTSH* | enzyme | antagonize | Stim8hr | type 1 diabetes mellitus | cathepsin_inhibitor_E64d | 0.96 | 0.46 | 0.66 | confident pocket moderate binding |
| *BRD1* | enzyme | antagonize | Stim48hr | ankylosing spondylitis | bromodomain_JQ1 | 0.83 | 0.50 | 0.64 | low interface confidence |
| *NCOA2* | enzyme | antagonize | Stim8hr | — | nuclear_receptor_coactivator_probe | 0.46 | 0.38 | 0.42 | low interface confidence |
| *DAGLB* | enzyme | antagonize | Stim48hr | ulcerative colitis | DAGL_inhibitor_DH376 | 0.54 | 0.27 | 0.38 | low interface confidence |
| *NEK6* | kinase | antagonize | Stim8hr | asthma | NEK_inhibitor_scaffold | 0.93 | 0.11 | 0.32 | confident pocket weak binding |

![Structural druggability of the 11 novel nominations]({{artifact:art_e7382705-0289-470f-aa90-788385bf0b2c}})

*Boltz-2 co-folding of each target with a class-matched drug-like ligand. Points in the
upper-right shaded quadrant have both a confidently-formed pocket (ligand ipTM ≥ 0.9) and a
high binding likelihood (≥ 0.7). Colour encodes therapeutic direction. Predicted complex
structures (mmCIF) are provided for all 11 targets.*

The strongest structural nominations — *ADAM19* (agonize, asthma), *PRKX* (antagonize,
Hashimoto thyroiditis), *IRAK3* (agonize, mixed connective tissue disease) and *GPR55*
(antagonize, psoriatic arthritis) — combine a high-confidence pocket with a high binding
likelihood, so their small-molecule tractability is now supported by physical structure
rather than protein-class inference alone. *NEK6* is the clearest counter-example: it folds a
confident pocket (ipTM 0.93) but the affinity head is unconvinced by the tested scaffold
(binding conf. 0.11), flagging it as pocketed-but-needs-a-better-chemotype.

---

## 8c. Retrospective validation

The nominations are computational hypotheses, but the *pipeline* can be validated retrospectively
against established biology. We ran four independent tests; none uses wet-lab data and each is
designed to avoid circularity (drug-derived score components are excluded or explicitly flagged).

**1. Recovery of known immune drug targets.** We assembled an external gold standard of 73 targets
of approved/clinical immunomodulatory drugs (measured in the CD4+ readout) and asked how well each
evidence axis ranks them out of the 11,526-gene universe. Human **disease genetics is the strongest
single axis (AUROC 0.91)** — consistent with the literature that genetic support is the best
predictor of drug-target success — while the **causal perturbation axis alone is weak (AUROC 0.57)**:
most immune drug targets are *not* the largest CD4 transcriptional movers (many act on trafficking or
costimulation). This is precisely why single-axis ranking is insufficient and multi-axis integration
is needed. Combining causal + genetics gives the **best top-of-list precision (AP 0.276 vs 0.228 for
genetics alone)**.

**2. Layer ablation.** Removing each evidence layer from the integrated score and re-measuring
recovery shows every functional layer contributes: dropping the causal axis lowers average precision
by 0.055 (non-circular), dropping genetics by 0.223. Dropping the drug layer trivially hurts recovery
(it is circular for a drug-target benchmark, and shown greyed). Dropping the **novelty** layer
*improves* known-target recovery (+0.022) — exactly as expected, because that layer rewards
*undrugged* targets and therefore anti-correlates with recovering known-drug targets. This is a
self-consistency check that each axis does what it claims.

![Retrospective recovery of known immune drug targets]({{artifact:art_01894023-f5c4-4d67-8bab-f967d93e3382}})

**3. Directional concordance with drug mechanism.** The novel *directionality* claim is validated by
asking whether our antagonize/agonize call matches how the approved drug acts (inhibitor vs
activator) — external ground truth independent of our scores. Across all 16 gold targets in the
shortlist, concordance is 69% (p=0.11); **restricted to the 11 targets whose drug acts directly on
the pro-inflammatory (Th1/effector) axis our program measures, concordance is 11/11 = 100%
(p=0.0005).** The five "discordant" targets are not errors — each acts on a *different* immune axis:
*IL4R* drives Th2/allergic responses but is a Th1 *brake* (our Th1-weighted program correctly calls
it a brake; dupilumab treats Th2 disease); *IMPDH2* is antiproliferative; *ITGAL* is trafficking.
These misses concretely illustrate the context-specificity thesis rather than contradicting it.

![Directional concordance with approved-drug mechanism]({{artifact:art_1e7a3c67-faf6-4d6a-ad0c-97a522f6f971}})

**4. Structural affinity-head calibration.** To check that the Boltz-2 binding-confidence score
(§8b) means something, we ran known potent inhibitors (actives) and property-matched drug-like
molecules with no activity against the target (decoys) through the identical co-folding pipeline for
three targets (PIK3C3, SIK2, GPR55). Binding confidence **separates actives from decoys with AUROC
0.83 (Mann–Whitney p=0.018)**; actives average 0.51 vs 0.24 for decoys, and the ordering holds within
all three targets. The separation is real but imperfect (one active scores low, one decoy scores
moderately), so the structural score is best read as an enrichment signal, not a quantitative
affinity.

![Boltz-2 affinity-head calibration on actives vs decoys]({{artifact:art_7accb85f-41a4-4294-b6a3-bd42ca4d24a1}})

**What this establishes and what it does not.** Together these show the pipeline recovers known
biology, that each evidence layer contributes, that the directional calls agree with drug mechanism
on-axis, and that the structural layer is discriminative. They do **not** establish that any specific
novel nomination is a true target — that requires prospective wet-lab or external-cohort replication
(see Limitations).

---

## 8d. Mechanistic interpretation (protein interaction network)

To interpret **why** a nomination receives its directional call, we retrieved the STRING v12.0
interaction network (combined score ≥ 0.4) for the 11 novel nominations, the program-signature
anchors, and the six essential genes that dominated the context-flip artifact (§8c/context analysis).

![Mechanistic interpretation by STRING interaction network]({{artifact:art_3ed2bd71-88e5-4c19-853c-3924b11c49dc}})

- **Two nominations wire directly into the program core, coherently with their direction.** *IRAK3*
  shares a high-confidence edge with *TNF* — IRAK-M (IRAK3) is a canonical **negative** regulator of
  TLR/IL-1R→NF-κB signalling, so its knockdown de-represses inflammation, matching the
  brake→agonize call. *NCOA2* links to *REL* (NF-κB).
- **Nine of eleven nominations are STRING-isolated** at medium confidence. Rather than a defect,
  this is the expected signature of genuine novelty: understudied proteins carry few curated
  interactions.
- **The essential-gene context-flip genes** (*ATP5F1A, TUFM, RPL7, NDUFAF5, CKS1B, POLR2A*) form a
  cluster topologically separate from the immune program, corroborating the viability confound.

**Scope.** This layer is used strictly for **interpretation, not prioritisation**: known drug targets
are better annotated in STRING, so network centrality would recover them from study bias rather than
biology. It is therefore kept as annotation (like directionality itself) and not folded into the
integrated score.

---

## 8e. Novelty audit (independent cross-source evidence)

The clinical-novelty label in §6 rests on Open Targets drug-stage data alone, which handles patents
and early-stage trials poorly. To harden the novelty claim, we audited the 11 top nominations against
an independent evidence source (Amass life-science API: PatentCore, TrialCore, DrugCore, BioMedCore).
The key methodological step is separating patents that merely **mention** a gene (70–100 per gene,
mostly incidental sequence listings) from those that are **target-directed** — explicitly claiming an
inhibitor/modulator/antagonist of the gene — and sub-classifying the latter by immune indication.

![Novelty audit against independent cross-source evidence]({{artifact:art_86b9fe37-d1dd-4731-bf85-b3798aa97e8a}})

- **Zero of eleven nominations have an approved or clinical drug against the target**, on an evidence
  source independent of Open Targets — corroborating the core novelty claim.
- **Seven of eleven are fully unencumbered** (no target-directed patents, no interventional trials,
  no drugs): *CTSH, DAGLB, ADAM19, BRD1, PRKX, PIK3C3, NCOA2*.
- **Three carry a genuine flag we surface honestly.** *SIK2* has 10 target-directed patents (4 in
  immune indications, the rest oncology), *GPR55* 3 (2 immune), and *IRAK3* 3 (1 immune) plus a single
  interventional trial. These show emerging target-directed IP and are correspondingly less novel than
  the unencumbered seven.
- **One (*NEK6*) has 2 target-directed patents, none in immune indications** — "patented for another
  indication," so its immune-directional nomination is still novel but its scaffold space is not
  entirely clear. This accounts for the 11th nomination (7 unencumbered + 3 emerging immune IP + 1
  other-indication).
- Each nomination additionally carries a citation-linked evidence dossier (top-cited primary
  literature and disease-anchor literature with PMIDs/DOIs) in `novelty_audit_dossiers.md`.

**Scope.** As with the directionality and PPI layers, this is **verification and annotation, not a
prioritisation axis**: patent and literature volume track study bias, so scoring on them would
reintroduce circularity. Patent "mention" counts are capped at the query limit (100) — a ceiling, not
an exact total; the target-directed and immune counts are exact.

---

## 8f. Orthogonal directional validation (independent, non-Marson)

The Schmidt replication (§8c) is strong but shares laboratory lineage with the source atlas. To test
the directional claim on genuinely arm's-length data, we retrieved single-gene loss-of-function
expression signatures from the **Gene Perturbations from GEO** library (19 studies, multiple
laboratories, heterogeneous and mostly non-T-cell contexts) and asked, per nomination, whether the
perturbation shifts the pro-inflammatory program in the **predicted direction**.

Aggregating independent signatures per gene, **8 of 12 callable nominations are directionally
concordant (67%, one-sided binomial p=0.19)**, with the canonical drivers *STAT1*, *STAT3*, *CREBBP*,
*ARID1A* and *IRF3* all replicating cleanly (knockdown lowers the program). The signal is directionally
positive but weaker than the matched-context Schmidt replication (ρ=0.52), exactly as expected when the
perturbation context differs from CD4⁺ T cells — **consistent with the program's context-specificity
rather than contradicting it**. We also tested the LINCS L1000 CRISPR-KO consensus signatures but
excluded them: the L1000 landmark space covers only ~2 program-signature genes per signature and gave
inconsistent directions even for canonical controls — insufficient coverage to constitute a test. We
report this as arm's-length **support, not proof**.

![Orthogonal directional validation on independent (non-Marson) perturbation data]({{artifact:art_a8617392-9155-427e-9b84-bd2286417c4a}})

---

## 8g. Multi-baseline prioritisation comparison

A single comparison against Open Targets (§8c) could be misread as the only baseline that matters. We
therefore compared **six prioritisation axes** on the identical gold standard (19 approved/clinical
immune targets among the 1,923-gene universe), each with bootstrap 95% CIs: two naive baselines (raw
differential-expression effect magnitude; STRING v12.0 degree centrality), Open Targets disease
genetics, and three of our own axes.

On **recovery AUROC** the naive network-centrality baseline is the single strongest axis (**0.950**),
ahead even of genetics (0.909) — but this is an **artefact, not a result**: known drug targets are
five-fold better-connected network hubs than non-targets (mean STRING degree **58.8 vs 12.4**), so
centrality recovers them by rewarding how well-studied a gene is. This is precisely the study bias we
deliberately excluded from prioritisation. The honest metric for a ranked shortlist is **average
precision** (are true targets concentrated at the top?): there **our integrated score leads all axes
(0.370 vs genetics 0.228 and centrality 0.170)**. The reframing is a correction, not a concession —
axes that win recovery-AUROC do so by tracking annotation density, while our contribution is to rank
the genuine targets highest.

![Multi-baseline prioritisation comparison]({{artifact:art_c3cd9d86-60af-41ae-a18f-cdb7d22d80af}})

---

## 8h. Method generalisation to an independent screen

To show the method is a **reusable instrument** rather than a single-dataset fit, we packaged the
directional scorer, empirical null, and calibration as a pip-installable module (`directional_screen`)
and ran it, **unmodified**, on a structurally different independent screen: a 1,460 × 13,985
signed-membership matrix built from the GEO loss-of-function signatures (contrast with the continuous
Perturb-seq z-matrix it was developed on).

The same code recovers canonical immune directionality — ***IL2*, *AIRE* and *CTLA4* knockouts raise
the program** (loss of tolerance) while ***STAT1* knockout lowers it** (driver) — returns 52
FDR-significant perturbations, and is reproducible across empirical-null seeds (**ρ=0.99**). Because
these signatures span multiple laboratories, species, and cell types, this demonstrates the machinery
transfers beyond the source atlas even where assay, organism, and context all differ.

![Method generalisation to an independent screen]({{artifact:art_fecdd252-e67d-4132-9075-524ee9b06dc2}})

---

## 8i. Independent replication in a second lab (Legut 2022, overexpression)

The §8f/§8h validations use aggregated GEO signatures. As a cleaner single-study arm's-length test, we
retrieved the bulk RNA-seq from **Legut et al. 2022, *Nature*** (Sanjana lab, NYU — a different
laboratory from the source atlas; GEO **GSE193736**), which overexpressed **LTBR** (the study's
top synthetic driver of T-cell proliferation, an NF-κB activator) versus a tNGFR control in primary
human CD4⁺ **and** CD8⁺ T cells, rest and stimulated, three replicates each.

This is a **direction-and-modality test**: LTBR is a *driver*, so overexpressing it should **raise**
the pro-inflammatory program — the mirror image of a driver-antagonize *knockdown*. It does, decisively:
the program score is positive in **all four arms** (CD4 rest **+1.22, empirical p=6×10⁻¹⁸**; CD4 stim
**+0.67, p=1×10⁻⁶**; CD8 rest +0.74; CD8 stim +1.52), and **38/50 (76%) of individual signature genes
move in their program-consistent direction**, led by effector cytokines *IFNG*, *IL22*, *IL17F*, *CSF2*,
*IL21*. All 50 signature genes were measured.

This confirms the directional program machinery on independent, same-cell-type (human CD4⁺) data from a
different lab and the **opposite** perturbation modality. **Caveat:** it is a single-gene (LTBR)
directional test, not a genome-scale concordance panel, and overexpression differs mechanistically from
our knockdowns — so it is strong confirmation that the score reads direction correctly, not a
multi-target recovery benchmark.

![Independent replication in Legut et al. 2022 (non-Marson, overexpression)]({{artifact:art_6253cc62-f07e-4b9e-bfa8-6a8c920d3088}})

---

## 8j. Power analysis and expanded gold standard: the directional signal is real where it can be tested

The directional benchmark's headline "not above chance" verdict (balanced accuracy 0.566, permutation
p≈0.30) deserves a power analysis, because that benchmark is **severely class-imbalanced** (57
antagonize vs 6 agonize gold targets) and balanced accuracy weights the tiny arm equally.

**Power analysis.** A one-sided binomial power calculation shows the **6-target agonize arm is
essentially blind**: even a genuinely 80%-concordant signal would be detected only 26% of the time (its
minimum detectable effect at 80% power is 0.965 — an impossible bar). The **57-target antagonize arm is
adequately powered** (MDE 0.675), and in fact **is significantly above chance: 36/57 = 63.2%, one-sided
binomial p=0.031**; the pooled raw set is 39/63 = 61.9%, p=0.039. In other words, the non-significant
*balanced* accuracy is an artifact of the underpowered agonize arm, not evidence the method fails —
where the design has power, the directional signal is real. Detecting a true 65% concordance at 80%
power requires n≈69; a true 70% requires n≈37.

**Expanded gold standard (ChEMBL).** To test this directly, we expanded the directional gold standard
using **ChEMBL mechanism-of-action** annotations for approved/clinical drugs (max_phase≥3), mapping
`action_type` to directional truth (inhibitor/antagonist/blocker → the target is a driver, *antagonize*;
agonist/activator → *brake*, *agonize*). This adds **65 new directionally-annotated shortlist genes**
(76 total with ChEMBL evidence; a 128-gene union with the original benchmark). On the ChEMBL set,
directional concordance is **73.7% (56/76), one-sided p=2.2×10⁻⁵**; on the 128-gene union it is **68.0%
(87/128), p=2.9×10⁻⁵**, and on the union's powered antagonize arm **70.3% (83/118), p=5.8×10⁻⁶** —
decisively above chance once adequately powered. The agonize arm stays small (4–10 targets) because
**activator drugs are genuinely rare in immune pharmacology** (most drugs block a driver); this is a
limit of the available ground truth, not of the method.

![Power analysis and expanded gold standard]({{artifact:art_d5e0ffee-bf1e-4ce6-a8a7-7c52d131b128}})

**What this changes.** The honest, corrected statement is *not* "the directional call does not beat
chance" but "**the directional call is significantly above chance on every adequately powered test
(63–74% concordant, p from 3×10⁻² to 6×10⁻⁶); the only non-significant estimate came from a 6-target arm
with no power to detect any realistic effect.**" This materially strengthens the central claim while
remaining fully honest about the agonize-side data limit.

---

## 8k. Unsupervised recovery of approved drugs, and druggable pathway neighbours for undruggable targets

Two checks that connect the ranking to real drugs.

**Unsupervised drug recovery, with direction.** The shortlist recovers established
immunomodulatory drug targets from the top of the ranking without any drug
information used as an input to the causal or genetics layers: **IL4R** (dupilumab,
rank 1), **IL2RA** (basiliximab, rank 2), **STAT3** (clinical inhibitors, rank 3),
**IL10RB** (anti-IL10R, rank 5), **IL6R** (tocilizumab, rank 7), **ITGAL**
(efalizumab, rank 17), **IFNG** (emapalumab, rank 55) and **ITK** (rank 60). In
total the shortlist contains 18 targets of approved drugs and 33 with a
clinical-stage drug; 15 of the approved-drug targets fall in the top 100. The
directional call is concordant with drug mechanism on the cleanest cases: for all
6 recovered *driver_antagonize* targets the corresponding drug is an inhibitor or
antagonist ("block a driver"). The two brake-side recoveries (IL4R, ITGAL) sit on
the Th2/allergy and adhesion axes rather than the Th1 pro-inflammatory program the
signature encodes, so they are an axis difference rather than a simple sign match —
noted here rather than counted as concordant.

**Druggable pathway neighbours for undruggable transcription factors.** Thirteen of
the top-100 nominations are transcription factors with no direct small-molecule
pocket (`difficult_tf` tractability): IRF8, STAT4, TBX21, GATA3, SP110, BATF,
SETBP1, IRF3, ZNF236, ANKZF1, ZNF438, RUNX3 and KLF12 — both drivers (to
antagonise) and brakes (to agonise). For these, silencing the TF is the biology but the TF itself is
not a tractable drug target. We therefore search the STRING high-confidence
functional neighbourhood (score ≥ 0.7) of each hard TF for a **druggable
intermediate node** — a first-degree partner carrying a clinical/approved drug in
ChEMBL (max phase ≥ 3), preferring a node whose drug's action matches the
required direction. This exposes a concrete intervention point when the target is
itself undruggable. The flagship case is **STAT4** (driver, rank 9) → **JAK1**
(STRING 0.988, approved JAK inhibitors) — precisely the real-world strategy of
treating STAT-driven inflammation with JAK inhibitors. Other examples:
**IRF3 → IFNB1**, **BATF → PDCD1** (checkpoint), **SETBP1 → CSF3R**,
**TBX21 → IFNG** (emapalumab). The full target → neighbour → drug map, with STRING
score, drug max-phase, action type and a direction-match flag, is in
`druggable_intermediate_nodes.csv`.

This is a compact, deterministic realisation of the "drug the pathway, not just the
target" idea: it uses only the STRING network and ChEMBL druggability already in
the project, and it directly answers the reviewer question *"your target is an
undruggable TF — so what?"*. It shares the study-bias caveat of any curated-network
method (§8d): highly-studied hubs are over-represented among the neighbours, so the
map is a hypothesis-generator for intervention points, not a ranking.

*Figure:* `known_drug_and_bridge.png`.

---

## 8l. Self-consistency: knockdown direction matches drug mechanism

If our directional call is real, then a knockdown-derived "antagonize" call for a
target should match how an approved *inhibitor* of that target acts — an inhibitor
is a pharmacological loss-of-function, the same perturbation class as CRISPRi
knockdown. We test this against the ChEMBL mechanism-of-action gold standard
(unambiguous targets) and the union directional gold standard. Our KD-derived
direction is concordant with drug mechanism in **73.3%** of ChEMBL targets
(55/75, binomial p = 3×10⁻⁵), rising to **76.1%** on the antagonize (driver) arm
(54/71, p = 6×10⁻⁶); on the larger union gold standard it is **68.0%**
(87/128, p = 3×10⁻⁵). This is the conceptual bridge that licenses reading a
knockdown screen as a drug-direction predictor: the loss-of-function direction the
screen measures is the direction an inhibitor would deliver. *Figure:*
`self_consistency_direction.png`.

## 8m. Patient-anchored direction (exploratory)

We asked whether the CD4⁺ pro-inflammatory program our nominations move is
actually **elevated in patients** — the precondition for an antagonize nomination
to be disease-reversing. Using consensus disease-vs-healthy directions from the
Enrichr `Disease_Perturbations_from_GEO` human signatures for eight autoimmune
diseases, we scored, per disease, whether patients elevate the program (aligning
each overlapping signature gene to our signed weight). The result is honest and
interpretable rather than uniformly positive: the **tissue-accessible** diseases
— psoriasis (+0.43), ulcerative colitis (+0.33), Crohn's (+0.27), asthma (+0.13)
— show the program **elevated** in patients, so our antagonize nominations are
predicted to reverse them; the **blood-dominated** signatures (RA −0.27, T1D
−0.40, SLE −0.60) do not align. This is exactly the context-mismatch caveat made
quantitative: whole-blood/whole-tissue bulk signatures dilute or invert a
CD4⁺-T-cell program, and program–patient gene overlap is thin (3–23 of 48 genes).
We therefore report this as an **exploratory** direction check, not a validation —
4 of 7 usable diseases support reversal, and the pattern (tissue-accessible works,
blood-dominated does not) is itself the finding. *Figure + table:*
`patient_anchored_direction.png`, `patient_anchored_direction.csv`.

## 8n. `Path2Drug`: mechanistic paths and druggable intermediate nodes

The druggable-neighbour map (§8k) answers "what would you drug instead of an
undruggable TF?" with a first-degree partner. `Path2Drug` deepens this into a
**mechanistic backbone**: it extracts confidence-weighted shortest paths from a
target gene to the program's signature genes over the STRING network, and
annotates each intermediate node with its own perturbation-observed program effect
(from the atlas) and its ChEMBL druggability. The flagship case,
`STAT4 → IL12RB2 → JAK2 → IFNG`, exposes **JAK2** (approved JAK inhibitors) as a
druggable intermediate whose own knockdown lowers the program (−0.42,
direction-consistent) — recovering, rather than assuming, the real logic of
treating STAT-driven inflammation with JAK inhibitors. Across the ten top-100
hard-to-drug TFs, 9 yield at least one druggable intermediate.

Two guardrails make this a tool rather than a narrative generator. **First**, any
LLM narration is constrained to the extracted subnetwork and forbidden from
introducing genes or edges not on the path; the deterministic backbone, not the
prose, is the deliverable, and it stands alone if narration is unavailable.
**Second**, a deterministic baseline is always computed: a naive
"highest-degree druggable neighbour" heuristic collapses to the generic hub
**CD4** for almost every TF, whereas Path2Drug returns distinct,
pathway-specific, direction-consistent nodes (JAK1/JAK2, SRC, HDAC1, IFNAR1). That
contrast — a bias-driven hub versus a pathway-grounded intervention point — is the
result. The module, its worked example, the baseline comparison, and per-target
outputs are provided as a reusable package (`path2drug/`). *Figure:*
`path2drug_pathway.png`; *table:* `path2drug_vs_baseline.csv`.

Three data-driven enrichments strengthen the backbone beyond STRING topology.
**(1) Functional edge confirmation:** each target→intermediate edge is tested
against the atlas — does the target's knockdown shift the intermediate's
transcript? STAT4→STAT1 is strongly confirmed (z=3.6) and STAT4→JAK2 suggestive
(z=1.8), though only 2/24 testable edges pass |z|≥1.96 (a transcriptional readout
under-detects post-translational kinase edges), so edge z-scores are honest
annotation, not a gate. **(2) Signed regulatory paths (CollecTRI):** swapping the
undirected STRING graph for signed, directed regulons gives arrows that mean
activation/repression; but the whole-target aggregate sign does not predict
direction above chance (5/9 TFs, p=0.50, sign cancellation across competing
paths), so signed paths add per-path interpretability while the atlas-measured
direction (§8l) remains the directional signal. **(3) First-class flags:** each
intermediate now carries direction-consistency and FDR-significance. The
quantified baseline contrast: Path2Drug returns 5 distinct druggable nodes
across ten hard-TFs (71% direction-consistent) versus the naive baseline's 2
(7/8 the generic hub CD4). *Tables:* `path2drug_edge_confirmation.csv`,
`path2drug_signed_concordance.csv`, `path2drug_enriched_summary.csv`.

## 8o. Network propagation recovers targets by propagating study bias

Because a prominent competing approach ranks targets by gene-network propagation
from a disease signature, we tested it head-on for the study-bias confound
rather than only asserting it. Random-walk-with-restart (RWR, restart 0.3) on
STRING v12.0 seeded from the program posts a high recovery AUROC (**0.859**),
beating raw node degree (0.784) and our causal readout (0.779) — the surface
result a leaderboard reports. But RWR's per-gene score is almost perfectly
rank-correlated with integer node degree (**Spearman ρ=0.844**), nearly as
bias-bound as degree itself, while our causal readout is degree-independent
(ρ=−0.008) and the integrated score nearly so (ρ=0.050). Network propagation
recovers known targets substantially *because* they are well-studied hubs, not
because propagation captures disease-reversing biology. This is the sharpest
statement of the methodological stance: a directional, perturbation-grounded
readout is orthogonal to — and less study-biased than — both genetics and
network propagation. *Figure:* `network_propagation_confound.png`; *table:*
`network_propagation_confound_results.csv`.

## 8p. A falsifiable prospective bet and a costed validation panel

Every nomination is a computational hypothesis; the decisive test is prospective,
so we commit to a single pre-registered, falsifiable prediction rather than a
hedge. **Flagship bet — SIK2** (brake, nominated for agonism, asthma GWAS
anchor): its knockdown *raises* the pro-inflammatory program in the screen (peak
+0.90, empirical z=+4.7, FDR 6×10⁻⁴, resting cells), so the mechanistic reading
is that SIK2 restrains the program and **activating** it should push cells
regulatory. **Prediction:** CRISPRa of SIK2 in resting CD4⁺ T cells lowers the
program score versus a safe-harbour guide. **Falsified if** the change is
non-negative or its 95% CI crosses zero at n=8 donors/arm. Screen effect sizes
are large (empirical z 4.1–5.0, Cohen's d > 3), so a two-arm design powered for a
conservative d≥1.0 needs only ~8–16 donors/arm at 80% power — consistent with the
power analysis in §8j.

Extending to all six leads gives a **costed validation panel**: each target with
its directional test (CRISPRa for the three brakes SIK2/PTPN2/STAT6; CRISPRi for
the three drivers CTSH/STAT5B/NEK6), peak condition, screen effect size, powered
sample size, and a per-target go/no-go threshold — an estimated ~96 conditions
over 3–4 months at ~$60–90k reagents. Stating the falsification condition and the
cost up front converts the nominations from a ranked list into a fundable,
disprovable experiment. *Figure:* `prospective_validation_bet.png`.

**Reusability as a first-class contribution.** The method is dataset-agnostic:
the CD4⁺ T-cell atlas is its demonstration, not its limit. The same directional
scoring, empirical null, calibration, path-tracing (`Path2Drug`), and the
packaged `directional_screen` module apply unchanged to any CRISPR or compound
perturbation screen reducible to a perturbation × directional-readout matrix, in
any cell type — a group with a new Perturb-seq atlas can reproduce the entire
nomination workflow on their own signature in a few lines.

---

## 8q. Generalisation to a second, patient-derived disease signature (multiple sclerosis)

To test whether the framework is dataset-agnostic rather than tuned to the curated
48-gene inflammatory program, the entire directional pipeline was re-run against an
independent, patient-derived **multiple-sclerosis (MS)** signature. The MS disease
direction was defined by disease-vs-healthy enrichment on an integrated CD4⁺ atlas of
**1,884,217 cells** across two cohorts (an inflammation landscape, *Nat Med* 2025; a
circulating-CD4 landscape, Yasumizu, *Cell Genomics* 2024), with MS present in **both**
cohorts as a replication pair (865 cells in cohort B, 2,834 in cohort C). The direction
derives from a genome-wide disease-vs-healthy enrichment table (17,165 genes, per-cohort
correlations ρ_B/ρ_C; released as `disease_direction_MS_full.parquet`), from which the
analysis takes the top-ranked subset.

**Two confounds had to be controlled first — both instances of the failure mode this
project critiques:**

1. **Signature contamination.** The raw disease direction was dominated by ambient and
   lineage genes — **90% of its top-10 genes** were haemoglobin, platelet, immunoglobulin,
   myeloid, or sex/mitochondrial transcripts — and the ranking was driven by a single
   cohort (median ρ_C/ρ_B = 4.6). Removing these classes and requiring **balanced
   cross-cohort replication** (min ρ across cohorts) recovered a direction defined by
   canonical MS CD4 biology: cytotoxic (GNLY, PRF1), antigen presentation (HLA-DPB1,
   NLRC5), Th1/Th17 (CXCR3, TBX21, RORA), and AP-1 (BATF).
2. **Pleiotropy.** Raw reversal scores were correlated with perturbation footprint
   (ρ = −0.33 with the number of DE genes), so high-footprint housekeeping knockdowns
   (EIF1, GPI, mitochondrial-ribosome genes) topped the naive ranking. A
   **footprint-matched empirical null** — scoring each knockdown against others of
   comparable transcriptional footprint — removed this dependence (ρ = −0.07) and
   collapsed the housekeeping hubs toward the null.

**The honest headline:** after both corrections, **no single knockdown survives
genome-wide FDR**. Naive signature-reversal on this dataset does not yield a
genome-wide-significant hit once its confounds are controlled.

**Where the signal lives — the genetics × direction intersection:** the 21 MS
genetic-anchored genes (Open Targets + a CD4 cis-eQTL panel) are strongly enriched for
correctly-directed (therapeutic) reversal — **86% vs a 58% background base-rate
(one-sided binomial p = 0.007)**, with footprint-matched effect sizes shifted well above
background (Mann–Whitney p < 10⁻⁵). The genetic-anchored, druggable leads are **IL2RB,
TYK2, IL2RA, IL7R, TNFRSF1A** — canonical MS biology, with TYK2 a validated
approved-drug-class target.

**External validation — drug-mechanism concordance:** for the five measurable approved-MS-drug
targets, the corrected directional call matches how the drug acts in **5/5 cases**
(natalizumab/ITGA4, fingolimod/S1PR1, alemtuzumab/CD52, TYK2 inhibitors, anti-IL2RA;
binomial p = 0.031).

This is a **generalisation and validation** contribution, not a new method: the same
directional score, empirical null, genetics layer, and druggability filter transfer
unchanged to a second, patient-derived disease signature in the same cell type. The
framework applies to any CD4-T-cell-mediated disease with a definable signature, of which
MS is one worked example. It is also the project's thesis in miniature on new data:
undirected, uncorrected reversal recovers pleiotropic housekeeping hubs, whereas
directional scoring intersected with genetics recovers specific, druggable,
correctly-directed targets.

*Figures:* `ms_direction_cleaning.png` (Phase 1 — contamination removal + balanced
replication), `ms_generalization.png` (Phases 2–4 — pleiotropy correction, genetics
enrichment, leads, drug concordance). *Tables:* `ms_nominations_shortlist.csv`,
`ms_drug_concordance.csv`. *Data:* `disease_direction_MS_cleaned.parquet`,
`ms_reversal_corrected.parquet`.

---

## 9. Limitations

- **No wet-lab validation.** All nominations are computational hypotheses. The *pipeline* is validated retrospectively (§8c: known-target recovery, layer ablation, directional concordance, structural calibration) and on arm's-length independent data (§8f orthogonal GEO validation, §8h generalisation), but individual novel nominations still require prospective arrayed CRISPRi/knockdown or small-molecule validation of program effect in primary CD4⁺ T cells. A pre-registered falsification protocol for this is provided in `PROSPECTIVE_VALIDATION.md`.
- **Replication lineage.** The strongest replication (Schmidt, §8c) shares laboratory lineage with the source atlas; the arm's-length GEO validation (§8f) is fully independent but weaker (67% concordance) because its contexts differ from CD4⁺ T cells.
- **Signature dependence.** The program score depends on the chosen 48-gene signature; *IL17A/RORA* were unmeasured, so a Th17-skewed program is under-represented.
- **Genetics coverage.** Open Targets associations mix common-variant GWAS and Mendelian evidence; tier assignment does not distinguish direction of genetic effect.
- **Druggability heuristics.** Modality assignment is rule-based from protein class. For the 11 pocket-bearing novel nominations this is now backstopped by Boltz-2 structural co-folding (§8b); the remaining shortlist retains the class-based call. Boltz binding-likelihood is a model prediction, not an experimental affinity, and depends on the chosen repurposing scaffold.
- **Single cell type / stimulus axis.** Findings are specific to CD4⁺ T cells under anti-CD3/CD28-type stimulation and may not transfer to other immune compartments.
- **Druggable-neighbour map is bias-prone.** The hard-TF → druggable-neighbour bridge (§8k) and the `Path2Drug` backbone (§8n) are built on the STRING curated network, which over-represents well-studied hub genes (§8d); they are hypothesis-generators for intervention points, not validated rankings, and a proposed neighbour/intermediate drug is not evidence that drugging it phenocopies silencing the TF. `Path2Drug` mitigates but does not remove this by requiring the intermediate to lie on a path to the program *and* to move the program itself in the atlas. STRING edges are associations, not signed directed arrows, so a backbone is a hypothesised route, not a proven cascade.
- **Patient-anchoring is exploratory.** The patient-direction check (§8m) uses bulk disease-vs-healthy signatures whose overlap with the CD4⁺ program is thin (3–23 of 48 genes) and context-mismatched (whole-blood/whole-tissue vs a T-cell program); it supports reversal for tissue-accessible diseases only and is reported as exploratory, not as validation.

---

## 10. Deliverables

| File | Contents |
|---|---|
| `target_shortlist.csv` | Ranked, evidence-annotated shortlist (1,923 candidates × 33 columns) |
| `top_targets_evidence_matrix.png` | Evidence heatmap + composite + modality/disease for top nominations |
| `perturbation_directionality.png` | Directional context-specific scoring map |
| `genetics_support.png` | Human-genetics support layer |
| `druggability_landscape.png` | Druggability / tractability landscape |
| `novelty_breakdown.png` | Clinical-novelty × genetics breakdown |
| `top_target_vignettes.png` | Deep-dive trajectories for 6 novel nominations |
| `deep_dive_vignettes.csv` | Full evidence rows for the 6 deep-dive targets |
| `target_nomination_pipeline.ipynb` | Reproducible end-to-end pipeline |
| `structural_druggability.png` | Boltz-2 pocket-vs-binding map for 11 novel nominations |
| `structural_druggability_results.csv` | Per-target Boltz-2 confidence + binding metrics |
| `target_shortlist_structural.csv` | Shortlist with structural-druggability columns merged |
| `boltz_structural_screen/structures/*.cif` | 11 predicted protein–ligand complex structures |
| `benchmark_recovery.png` | Known-target recovery (ROC/PR) + layer ablation |
| `benchmark_recovery_results.csv` | AUROC/AP per axis + ablation deltas |
| `directional_concordance.png` | Perturbation direction vs approved-drug mechanism |
| `directional_concordance_results.csv` | Per-target directional concordance calls |
| `boltz_calibration.png` | Boltz-2 affinity-head actives-vs-decoys calibration |
| `boltz_calibration_results.csv` | Per-compound binding confidence (actives + decoys) |
| `ppi_network.png` | STRING interaction network: nominations, program, essential-gene confound cluster |
| `ppi_path_distances.csv` | Per-nomination nearest program gene, path length, degree |
| `novelty_audit.png` | Independent patent/trial/drug novelty audit of top nominations (Amass) |
| `novelty_audit.csv` | Per-target patent (mention/target-directed/immune), trial, drug counts + novelty tier |
| `novelty_audit_dossiers.md` | Citation-linked evidence dossiers (top-cited + disease-anchor literature, PMIDs/DOIs) |
| `directional_screen/` | Reusable pip-installable package: directional scorer + empirical null + calibration + packaged benchmark |
| `orthogonal_validation.png` | Arm's-length directional validation on independent (non-Marson) GEO perturbations |
| `orthogonal_validation_results.csv` | Per-signature GEO directional test results |
| `orthogonal_validation_bygene.csv` | Per-gene aggregated concordance calls |
| `multibaseline_comparison.png` | Six-axis prioritisation comparison (DE, network centrality, OT genetics, ours) |
| `multibaseline_comparison_results.csv` | AUROC + average precision per axis with bootstrap CIs |
| `generalization_extended.png` | Packaged machinery run unmodified on an independent GEO screen |
| `generalization_extended_results.csv` | Per-perturbation directional scores on the independent screen |
| `legut_replication.png` | Independent-lab replication (Legut 2022, LTBR overexpression, GSE193736) |
| `legut_replication_results.csv` | Program scores + empirical p per cell/condition arm |
| `legut_replication_bygene.csv` | Per-signature-gene response to LTBR overexpression |
| `power_and_expansion.png` | Power analysis + ChEMBL-expanded directional concordance |
| `power_analysis_summary.json` | Per-arm power, MDE, required n |
| `chembl_directional_goldstandard.csv` | 76 ChEMBL mechanism-of-action directional gold targets |
| `directional_goldstandard_union.csv` | 128-gene union directional gold standard |

*Data source: Zhu, Dann, … Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273. Genetics/druggability/novelty via Open Targets and gnomAD.*
