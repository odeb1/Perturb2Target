# Perturb2Target — Project Outputs

Context-specific, directional, genetics-anchored, druggability-filtered drug-target
nominations from the genome-scale CD4+ T-cell CRISPRi Perturb-seq dataset
(Zhu, Dann, ... Pritchard & Marson 2025).

## Folder layout

- `reports/`   — findings report (`findings_report.md`, `findings_report.pdf`)
- `figures/`   — all publication figures (PNG, 200 dpi)
- `tables/`    — CSV deliverables (shortlist, per-layer annotations, novelty, manifest)
- `data/checkpoints/` — parquet checkpoints for each pipeline step (reload points)
- `codebase/`  — reproducible pipeline: `src/01..08` scripts, README, requirements,
                 notebook, and the full `target_nomination_codebase.tar.gz` bundle

## Key deliverables

| File | What it is |
|---|---|
| tables/target_shortlist.csv | Ranked, evidence-annotated shortlist (1,923 candidates x 33 cols) |
| tables/deep_dive_vignettes.csv | 6 curated novel deep-dive nominations |
| figures/top_targets_evidence_matrix.png | Evidence matrix for top nominations |
| figures/perturbation_directionality.png | Directionality map (agonize brakes vs antagonize drivers) |
| figures/genetics_support.png | Human-genetics support layer |
| figures/druggability_landscape.png | Druggability / tractability layer |
| figures/novelty_breakdown.png | Clinical-novelty filter |
| figures/top_target_vignettes.png | Deep-dive trajectory vignettes |
| figures/ml_demo_perturbation_regression.png | Supplementary ML demo (strength -> effect magnitude) |
| reports/findings_report.pdf | Methods + findings writeup |
| codebase/target_nomination_pipeline.ipynb | Reproducible pipeline notebook |

## Pipeline order (codebase/src/)

01 acquire+validate inputs -> 02 program signature -> 03 directional scoring ->
04/04a genetics (Open Targets) -> 05 druggability -> 06 clinical novelty ->
07 integrated ranking (+07b evidence figure) -> 08 deep-dive vignettes

Nominations are COMPUTATIONAL hypotheses; no wet-lab validation. See findings report
limitations.
