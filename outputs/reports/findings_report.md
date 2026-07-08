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

## 9. Limitations

- **No wet-lab validation.** All nominations are computational hypotheses; the intended next step is arrayed CRISPRi/knockdown or small-molecule validation of program effect in primary CD4⁺ T cells.
- **Signature dependence.** The program score depends on the chosen 48-gene signature; *IL17A/RORA* were unmeasured, so a Th17-skewed program is under-represented.
- **Genetics coverage.** Open Targets associations mix common-variant GWAS and Mendelian evidence; tier assignment does not distinguish direction of genetic effect.
- **Druggability heuristics.** Modality assignment is rule-based from protein class; it approximates but does not replace structural tractability assessment.
- **Single cell type / stimulus axis.** Findings are specific to CD4⁺ T cells under anti-CD3/CD28-type stimulation and may not transfer to other immune compartments.

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

*Data source: Zhu, Dann, … Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273. Genetics/druggability/novelty via Open Targets and gnomAD.*
