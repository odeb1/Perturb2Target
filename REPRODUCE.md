# Reproduction ledger — claims → evidence

*Every headline number in the paper and report, mapped to the exact script that produces it, the
artifact it is read from, and its verified value. A reviewer should be able to check any single claim
in under a minute. All values below were re-read from the saved artifacts, not transcribed from prose.*

**Data source:** Zhu, Dann, … Pritchard & Marson (2025), *bioRxiv* 10.64898/2025.12.23.696273
(`GWCD4i.DE_stats.h5ad`; conditions Rest / Stim 8 h / Stim 48 h). 33,983 perturbation×condition rows,
10,282 readout genes.

---

## How to regenerate

```bash
pip install -r requirements.txt
# Step 1 streams signature-gene log2FC from the 16.8 GB S3 object (does not download the whole file).
# Steps 2–8 run from the checkpointed parquets in data/checkpoints/ if the raw stream is unavailable.
python src/01_acquire_validate_inputs.py     # → de_obs/de_var parquets, z-score matrix
python src/02_program_signature.py           # → program_signature.csv + heatmap
python src/03_directional_scoring.py         # → perturbation_scores + directionality figure
python src/04a_open_targets_fetch.py         # → Open Targets genetics/tractability/drug pull
python src/04_genetics_support.py            # → genetics layer
python src/05_druggability.py                # → druggability layer
python src/06_clinical_novelty.py            # → novelty layer
python src/07_integrated_ranking.py          # → target_shortlist.csv (1,923) + integrated score
python src/08_deep_dive_vignettes.py         # → deep-dive vignettes
```

The validation, structural, replication, PPI, generalization, sensitivity, and novelty-audit
analyses are downstream of `target_shortlist.csv` / `candidates_step7.parquet` and are documented in
`target_nomination_pipeline.ipynb`.

---

## Ledger

Legend for **Source**: `src/NN` = pipeline script; `nb` = `target_nomination_pipeline.ipynb`;
artifacts are under `outputs/`.

| # | Claim (as stated in paper/report) | Verified value | Source script | Evidence artifact |
|---|---|---|---|---|
| **Positive controls** |
| 1 | Established immune targets recovered at top of ranking | IL4R 0.869 (#1), IL2RA 0.844 (#2), STAT3 0.832 (#3), IL6R 0.790 | `src/07` | `tables/target_shortlist.csv` |
| 2 | Quality-gated directional shortlist size | 1,923 candidates | `src/07` | `tables/target_shortlist.csv` |
| **Benchmark (recovery of 73 known immune targets)** |
| 3 | Causal effect axis (clean) | AUROC 0.585, AP 0.025 | `nb` | `tables/benchmark_recovery_results.csv` |
| 4 | Disease genetics (OT, partial circularity) | AUROC 0.909, AP 0.228 | `nb` | `tables/benchmark_recovery_results.csv` |
| 5 | Causal + genetics (clean) | AUROC 0.838, AP 0.276 | `nb` | `tables/benchmark_recovery_results.csv` |
| 6 | Full integrated (4 layers, mixed) | AUROC 0.820, AP 0.370 | `nb` | `tables/benchmark_recovery_results.csv` |
| 7 | Layer ablation — genetics most load-bearing | drop-genetics ΔAP −0.223; drop-drug −0.275; drop-novelty +0.022 | `nb` | `tables/benchmark_recovery_results.csv` |
| **External replication (Schmidt 2022 CRISPRi)** |
| 8 | Directional calls replicate on IFN-γ axis | ρ=0.52, p=1.15e-08, n=106, 73% agreement (binom p=1.7e-06) | `nb` | `replication_summary.json` |
| 9 | Directional calls replicate on IL-2 axis | ρ=0.49, p=2.5e-06, n=82, 74% agreement | `nb` | `replication_summary.json` |
| 10 | Non-canonical genes still replicate | ρ=0.39, p=1.0e-04, n=94, 69% agreement | `nb` | `replication_summary.json` |
| 11 | Specificity control (Shifrut proliferation readout, correct fail) | ρ=0.06, p=0.83, 44% agreement | `nb` | `figures/replication_shifrut.*` (project artifact) |
| **Honest negatives** |
| 12 | Directional benchmark not above chance on full set | balanced acc 0.566, permutation p≈0.30 | `nb` | `directional_benchmark_honest.csv` (project artifact) |
| 13 | Learned integration does not beat genetics alone | learned 0.813 < raw OT-genetic 0.821 | `nb` | `learned_integration_results.csv` (project artifact) |
| 14 | Context-specificity dominated by viability confound | 0/35 direction-flips are immune; log-breadth OR=1.40 p<1e-27 | `nb` | `context_confound_summary.json` (project artifact) |
| **Method generalization (independent assay)** |
| 15 | Recovery transfers to independent CRISPRi assay | Schmidt AUROC 0.571 vs Perturb-seq 0.574 | `nb` | `tables/method_generalization_summary.json` |
| 16 | Top-50 nominations overlap across assays | 26/50 shared, hypergeometric p=3.3e-31 | `nb` | `tables/method_generalization_summary.json` |
| 17 | Directional agreement across assays (modest) | 59.2% | `nb` | `tables/method_generalization_summary.json` |
| **Integration-weight sensitivity** |
| 18 | Ranking robust to ±40% weight perturbation | Spearman ρ=0.995; top-50 Jaccard 0.82 | `nb` | `tables/weight_sensitivity_summary.json` |
| **Structural druggability (Boltz-2)** |
| 19 | High-confidence predicted binders | ADAM19 0.969, PRKX 0.890, IRAK3 0.889, GPR55 0.837 | `nb` | `tables/structural_druggability_results.csv` |
| 20 | Affinity-head separates actives from decoys | AUROC 0.833, Mann-Whitney p=0.018 | `nb` | `boltz_calibration_results.csv` |
| **Novelty audit (Amass, independent evidence)** |
| 21 | No approved/clinical drug against any nomination | 0/11 | `nb` | `tables/novelty_audit.csv` |
| 22 | Fully unencumbered nominations | 7/11 (CTSH, DAGLB, ADAM19, BRD1, PRKX, PIK3C3, NCOA2) | `nb` | `tables/novelty_audit.csv` |
| 23 | Emerging target-directed IP (honest flag) | 3/11: SIK2 (10 td/4 immune), GPR55 (3/2), IRAK3 (3/1 + 1 trial) | `nb` | `tables/novelty_audit.csv` |
| 24 | Arm's-length directional validation (non-Marson GEO) | 8/12 concordant (67%, one-sided binom p=0.19); canonical drivers replicate | `nb` | `tables/orthogonal_validation_summary.json` |
| 25 | Network centrality wins recovery via study bias | STRING degree AUROC 0.950; positives 5× better-connected (58.8 vs 12.4) | `nb` | `tables/multibaseline_comparison_results.csv` |
| 26 | Integrated score leads on average precision | AP 0.370 vs OT-genetics 0.228, centrality 0.170, DE-magnitude 0.043 | `nb` | `tables/multibaseline_comparison_results.csv` |
| 27 | Packaged machinery generalises to independent screen | 52 FDR-significant on 1,460×13,985 GEO matrix; IL2/AIRE/CTLA4 KO raise, STAT1 lowers; ρ=0.99 across seeds | `nb` | `tables/generalization_extended_summary.json` |

---

## Caveats attached to specific claims

- **Circularity (#4, #6):** genetics and drug/novelty layers are partly derived from drug data;
  recovery of *drug targets* by them is a labelled circular ceiling, not a clean test. Only the
  causal and causal+genetics axes (#3, #5) are circularity-clean.
- **Shared lineage (#8–#10):** the Schmidt replication set shares laboratory lineage with the source
  Perturb-seq data (both Marson lab), so it is a strong but not fully arm's-length replication. The
  GEO orthogonal validation (#24) and generalisation (#27) ARE fully arm's-length (multi-lab GEO), but
  weaker (67% concordance) because their contexts differ from CD4⁺ T cells.
- **Recovery vs precision (#25, #26):** high recovery-AUROC baselines (network centrality, genetics)
  track annotation density / study bias; average precision is the honest metric for a ranked
  shortlist, and there the integrated score leads. Centrality is used for interpretation only, never
  as a prioritisation axis.
- **Patent counts (#21–#23):** "mention" counts are capped at the query limit of 100 (a ceiling, not
  exact); target-directed and immune counts are exact.
- **All nominations are computational hypotheses (#19–#23):** no wet-lab validation. See
  `PROSPECTIVE_VALIDATION.md` for the falsification protocol.
