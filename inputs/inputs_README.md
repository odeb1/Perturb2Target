# inputs/ — Curated reference tables fed into the pipeline

This folder holds the two hand-curated reference tables the analysis pipeline
**reads in but does not generate**. They are author-provided knowledge inputs,
distinct from the intermediate `data/checkpoints/` (produced and reloaded between
steps) and the final `tables/`, `figures/`, `reports/` outputs.

## Files

### program_signature.csv  (50 genes)
Defines the pathogenic pro-inflammatory CD4+ T-cell program the pipeline scores
perturbations against.

| Column | Meaning |
|---|---|
| gene | Gene symbol |
| weight | +1 = driver (pro-inflammatory); -1 = regulatory brake |
| arm | `pro_inflammatory` or `regulatory` |
| measured | Whether the gene is present in the readout (48/50 measured; IL17A, RORA absent) |

**Consumed by:** `src/02_program_signature.py`, then propagated into the
directional scoring in `src/03_directional_scoring.py`. This is the biological
definition the entire directional score is built on — the signed program a
perturbation is judged to push up (driver-like) or down (brake-like).

### druggable_class_map.csv  (5,037 genes)
Gene -> protein-class lookup used to assign drug modality and tractability.

| Column | Meaning |
|---|---|
| gene | Gene symbol |
| protein_classes | Semicolon-separated class tags (e.g. `catalytic_receptor;kinase`) |
| primary_class | Single dominant class (kinase, catalytic_receptor, GPCR, ...) |

**Consumed by:** `src/05_druggability.py` (modality assignment: kinase -> small
molecule, surface receptor -> antibody, etc.) and `src/06_clinical_novelty.py`.

## Not in this folder: the raw Perturb-seq data

The primary input to step 01 is the 16.8 GB `GWCD4i.DE_stats.h5ad` file
(Zhu, Dann, ... Pritchard & Marson 2025), streamed from
`s3://genome-scale-tcell-perturb-seq/marson2025_data/`. It is **too large to
store here** — `src/01_acquire_validate_inputs.py` streams it and writes the
small checkpoints `de_obs.parquet` / `de_var.parquet` into `data/checkpoints/`,
which every downstream step reads instead of the raw h5ad.

## Consumption map

| Input file | Read by |
|---|---|
| program_signature.csv | src/02_program_signature.py -> src/03_directional_scoring.py |
| druggable_class_map.csv | src/05_druggability.py, src/06_clinical_novelty.py |
| GWCD4i.DE_stats.h5ad (S3, not stored) | src/01_acquire_validate_inputs.py |
