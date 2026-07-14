# Perturb2Target

**🔗 Live demo: [perturb2target.streamlit.app](https://perturb2target.streamlit.app/)**

**From a genome-scale CD4⁺ T-cell CRISPRi Perturb-seq screen to directional, genetics-anchored, druggable drug-target nominations.**

Perturb2Target turns a genome-scale perturbation screen into a ranked, evidence-annotated
drug-target shortlist with four things most target-ID pipelines don't jointly provide:
**directionality** (antagonize a driver vs. agonize a brake), **human-genetics anchoring**,
**druggability / tractability** (including co-folded structures), and a **clinical-novelty**
filter. The method is dataset-agnostic and ships as a reusable package, an interactive
explorer, and a full write-up.

**Data source:** Zhu, Dann, … Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273
(`GWCD4i.DE_stats.h5ad`; conditions Rest / Stim 8 h / Stim 48 h; 33,983 perturbation×condition
profiles × 10,282 readout genes).

> **These are computational hypotheses** from multi-evidence triangulation — no wet-lab
> validation. Credibility rests on evidence convergence and on recovery of established immune
> drug targets (IL4R, IL2RA, STAT3, IL6R) at the top of the ranking as positive controls.

---

## What's in this repo

| Path | What it is |
|---|---|
| **`src/`** | The core pipeline as ordered scripts (steps 1–8). See table below. |
| **`target_nomination_pipeline.ipynb`** | The same pipeline as one annotated, checkpoint-driven notebook. |
| **`inputs/`** | Curated inputs: `program_signature.csv` (48-gene directional signature), `druggable_class_map.csv`, `inputs_README.md`. |
| **`data/checkpoints/`** | Parquet checkpoints written/consumed across steps (`de_obs`, `de_var`, `ot_annotation`, `candidates_step5/6/7`, `perturbation_scores_wide`). |
| **`outputs/tables/`** | Every results table (target shortlist 1,923×34, benchmark/validation CSVs, gold standards, summaries). |
| **`outputs/figures/`** | All standalone figures (33) from the analyses and validations. |
| **`outputs/paper/`** | The academic write-up: `main.tex` (ISMB-style), `paper.html`, `paper.pdf`, and `figs/` (19 numbered paper figures). |
| **`outputs/reports/`** | `findings_report.md`/`.pdf` (narrative), `novelty_audit_dossiers.md` (per-target patent/trial dossiers). |
| **`outputs/directional_screen/`** | **Reusable pip-installable package** — the scoring function, empirical null, calibration, and the versioned directional benchmark + loader, with a worked example. |
| **`explorer/`** | **Interactive Streamlit app** — directional map, context dynamics, method funnel, landscape, and the MS-generalization tab; co-folded 3-D structures. |
| **`boltz_structural_screen/`** | Boltz-2 co-folding of 11 pocket-bearing novel targets (YAMLs, `.cif` structures, `struct_druggability` results). |
| **`path2drug/`** | Network module that routes an undruggable TF hit to a druggable node on its regulatory pathway (STRING + CollecTRI). |
| **`MS_data/`** | Multiple-sclerosis generalization inputs + results (disease direction, corrected reversal, nominations, drug concordance, projection). |
| **`REPRODUCE.md`** | Step-by-step reproduction guide. |

## Core pipeline — `src/` (run in order)

| File | Step | Produces |
|---|---|---|
| `src/01_acquire_validate_inputs.py` | 1 · acquire/validate | `data/checkpoints/de_obs.parquet`, `de_var.parquet`, z-score matrix (streamed from S3) |
| `src/02_program_signature.py` | 2 · program signature | `inputs/program_signature.csv`, `program_signature_heatmap.png` |
| `src/03_directional_scoring.py` | 3 · directional scoring + empirical null | `perturbation_scores*.parquet/csv`, `perturbation_directionality.png` |
| `src/04a_open_targets_fetch.py` | 4a · Open Targets GraphQL fetch | `handoff/ot_deep.json` |
| `src/04_genetics_support.py` | 4 · genetics support | `data/checkpoints/ot_annotation.parquet`, `genetics_support.png` |
| `src/05_druggability.py` | 5 · druggability/tractability | `data/checkpoints/candidates_step5.parquet`, `druggability_landscape.png` |
| `src/06_clinical_novelty.py` | 6 · clinical novelty | `data/checkpoints/candidates_step6.parquet`, `clinical_novelty.csv`, `novelty_breakdown.png` |
| `src/07_integrated_ranking.py` | 7 · integrated ranking | `outputs/tables/target_shortlist.csv`, `candidates_step7.parquet` |
| `src/07b_evidence_matrix_figure.py` | 7 · figure | `top_targets_evidence_matrix.png` |
| `src/08_deep_dive_vignettes.py` | 8 · deep dive | `top_target_vignettes.png`, `deep_dive_vignettes.csv` |

**Integrated score** = 0.34·causal + 0.30·genetics + 0.22·druggability + 0.14·novelty, with
**direction kept as an explicit annotation** (never folded into the score). Quality gate →
1,923 candidates (1,290 driver→antagonize, 633 brake→agonize).

## Beyond the core pipeline

- **Validation** — retrospective concordance vs approved-drug mechanism (ChEMBL, 73%),
  known-target recovery, calibration, three-lab replication (Schmidt/Shifrut/Legut),
  orthogonal directional validation, multi-baseline comparison (vs Open Targets, DE ranking,
  network centrality), and power/effect-size analysis. Tables in `outputs/tables/`, figures in
  `outputs/figures/`.
- **MS generalization** — the entire pipeline re-run, unchanged, on an independent
  patient-derived multiple-sclerosis signature (1.9M CD4⁺ cells). After controlling
  contamination + pleiotropy confounds, MS genetic anchors are 86% correctly-directed
  (vs 58% background, p=0.007) and the direction matches approved-MS-drug mechanism 5/5
  (p=0.031). See `MS_data/` and the explorer's MS tab.
- **Structural druggability** — 11 novel pocket-bearing targets co-folded with Boltz-2
  (`boltz_structural_screen/`).
- **Path2Drug** — turns an undruggable TF hit into a druggable, direction-consistent target
  on its regulatory network (`path2drug/`).
- **Reusable resource** — `outputs/directional_screen/` is a standalone package
  (`pip install -e .`) with the scoring + empirical-null + calibration and a versioned
  directional benchmark, so the method runs on any perturbation × readout screen.

## Interactive explorer

```bash
cd explorer
pip install -r requirements.txt
streamlit run app.py
```
Opens at `http://localhost:8501`. Tabs: directional map · context dynamics · method funnel ·
landscape · MS generalization. All filtering runs on the real data.

## Layout & conventions

- Paths in `src/` scripts are **relative** to the project root. Create `handoff/` before
  running the OT fetch step, or run the notebook, which reads the saved checkpoints directly.
- **`handoff/`** — JSON passed between the Open Targets fetch loop and the analysis kernel
  (not committed; regenerated by step 4a).

## Provenance note

Each `src/NN_*.py` is the **reproduction code captured from artifact lineage** — the code that
actually produced each saved output, with absolute artifact-store paths rewritten to relative
`inputs/` / `data/` / `handoff/` paths. `04a_open_targets_fetch.py` is the one reconstructed
module: its Open Targets GraphQL calls ran through a platform MCP connector (not captured in
file lineage), so it ships both the connector form and a portable public-API HTTP fallback
against `https://api.platform.opentargets.org/api/v4/graphql`.

## Environment

- Core pipeline (steps 2–8): Python 3.11, `numpy pandas scipy matplotlib`.
- Step 1: streams the 16.8 GB `.h5ad` from S3 (`anndata`/`h5py`/`s3fs`/`boto3`).
- Explorer: `streamlit`, `plotly`, `py3Dmol` (see `explorer/requirements.txt`).
- Structural screen: hosted Boltz-2 API (see `boltz_structural_screen/README.md`).

See `requirements.txt`.
