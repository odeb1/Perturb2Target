# Context-specific, directional, genetics-anchored drug-target nomination — codebase

End-to-end pipeline that turns a genome-scale CD4+ T-cell CRISPRi Perturb-seq screen into a
ranked, evidence-annotated drug-target shortlist with **directionality** (antagonize drivers
vs. agonize brakes), **human-genetics anchoring**, **druggability/tractability**, and a
**clinical-novelty** filter.

**Data source:** Zhu, Dann, … Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273
(`GWCD4i.DE_stats.h5ad`; conditions Rest / Stim 8 h / Stim 48 h).

> These are **computational hypotheses** from multi-evidence triangulation — no wet-lab
> validation. Credibility rests on evidence convergence and on recovery of established immune
> drug targets (IL4R, IL2RA, STAT3, IL6R) at the top of the ranking as positive controls.

## Source files (run in order)

| File | Step | Produces |
|---|---|---|
| `src/01_acquire_validate_inputs.py` | 1 · acquire/validate | `inputs/de_obs.parquet`, `de_var.parquet`, `zscore_f32.npy` |
| `src/02_program_signature.py` | 2 · program signature | `program_signature.csv`, `program_signature_heatmap.png` |
| `src/03_directional_scoring.py` | 3 · directional scoring + empirical null | `perturbation_scores*.parquet/csv`, `perturbation_directionality.png` |
| `src/04a_open_targets_fetch.py` | 4a · OT GraphQL fetch | `handoff/ot_deep.json` |
| `src/04_genetics_support.py` | 4 · genetics support | `ot_annotation.parquet`, `genetics_support.png` |
| `src/05_druggability.py` | 5 · druggability/tractability | `candidates_step5.parquet`, `druggability_landscape.png` |
| `src/06_clinical_novelty.py` | 6 · clinical novelty | `candidates_step6.parquet`, `clinical_novelty.csv`, `novelty_breakdown.png` |
| `src/07_integrated_ranking.py` | 7 · integrated ranking | `target_shortlist.csv`, `candidates_step7.parquet` |
| `src/07b_evidence_matrix_figure.py` | 7 · figure | `top_targets_evidence_matrix.png` |
| `src/08_deep_dive_vignettes.py` | 8 · deep dive | `top_target_vignettes.png`, `deep_dive_vignettes.csv` |

`target_nomination_pipeline.ipynb` (top level) is the same pipeline as an annotated,
checkpoint-driven notebook.

## Layout & conventions

- **`inputs/`** — data checkpoints (parquet/npy) written and consumed across steps.
- **`handoff/`** — JSON passed between the MCP fetch loop and the analysis kernel
  (`deep_query_genes.json` → `ot_deep.json`).
- Paths in the scripts are **relative** to the project root; create `inputs/` and `handoff/`
  before running, or run the notebook which reads the saved artifact checkpoints directly.

## Provenance note

Each `src/NN_*.py` is the **reproduction code captured from the artifact lineage** of that
step's output — i.e. the code that actually ran to produce the saved artifact, with absolute
artifact-store paths rewritten to `inputs/` and `handoff/`. `04a_open_targets_fetch.py` is the
one reconstructed module: its Open Targets GraphQL calls ran through the Claude Science
`clinical-genomics` MCP connector (not captured in file lineage), so it ships both the
`host.mcp(...)` form and a portable public-API HTTP fallback.

## Environment

- Steps 2–8: Python 3.11, `numpy pandas scipy matplotlib` (+ the `figure-style` skill's
  `apply_figure_style()`/`panel_letter()` for figures).
- Step 1: streamed the 16.8 GB `.h5ad` from S3 (env `perturbseq`, `anndata`/`h5py`/`s3fs`/`boto3`).
- Step 4a: Open Targets GraphQL — via `clinical-genomics` MCP in-platform, or `requests` against
  `https://api.platform.opentargets.org/api/v4/graphql` externally.

See `requirements.txt`.
