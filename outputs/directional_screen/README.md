# directional_screen

**A small, assay-agnostic toolkit for _directional_ target nomination from massively-multiplexed perturbation screens.**

Standard screen analysis ranks perturbations by *effect magnitude*. `directional_screen` adds the axis that magnitude cannot supply: **direction** — is a gene a **driver** (knock it down and the phenotype goes down ⇒ a therapy should *antagonize* it) or a **brake** (knock it down and the phenotype goes up ⇒ a therapy should *agonize* it)? Direction is what decides drug modality, and human genetics is silent on it.

The toolkit was developed on CD4⁺ T-cell CRISPRi Perturb-seq (Zhu, Dann et al. 2025) but is deliberately assay-agnostic: it applies to any screen that yields a **perturbation × readout-gene** effect matrix — CRISPRi/a Perturb-seq, pooled cytokine-sort screens, or functional-metaviromics screens where the "program" is a signed host-response signature.

![How to reuse this on your own screen]({{artifact:art_3a6c8adf-627f-4af8-b4ef-d5361b014019}})

## Install

```bash
pip install -e .          # from this directory
# requires: numpy, pandas, scipy, scikit-learn
```

## Quickstart

```python
import directional_screen as ds

# zmatrix: DataFrame (perturbations × readout genes) of z-scores
# signature: Series gene -> +1 (pro-program) or -1 (anti-program)
scored = ds.score_screen(zmatrix, signature)      # signed score + empirical-null FDR + direction
print(scored[["program_score", "emp_z", "emp_fdr", "direction"]].head())

# score your directional calls against the packaged gold-standard benchmark
result = ds.evaluate_directions(scored["direction"])
print(result["balanced_accuracy"], result["permutation_p"])
```

Run the full worked example:

```bash
python examples/worked_example.py
```

## What's in the box

| Module | Functions | Purpose |
|---|---|---|
| `scoring` | `program_score`, `empirical_null`, `score_screen` | Signed program score + composition-matched empirical null → per-perturbation `emp_z`, BH-FDR, and a driver/brake call |
| `calibration` | `concordance_labels`, `fit_direction_confidence`, `reliability`, `brier_skill` | Calibrate P(direction correct) against an independent screen; report reliability + Brier skill so you learn the **directional resolution on your own data** |
| `benchmark` | `load_benchmark`, `load_dictionary`, `evaluate_directions` | A versioned, external directional gold standard (63 immune drug targets, labelled driver/brake from approved-drug mechanism) + an honest scorer |

## The packaged benchmark — and its honest negative

`data/directional_benchmark_v1.csv` labels 63 immune drug targets by the direction their approved drug implies (an inhibitor ⇒ the target is a driver to antagonize; an agonist/replacement ⇒ a brake to agonize). This is external ground truth **independent of any perturbation screen**.

It ships **with the honest negative** we found: on this set the reference framework's directional accuracy is 0.62 raw — *below* the 0.90 majority-class baseline — with balanced accuracy ~0.57 and permutation p ~0.30. Per-gene direction from a single screen has modest resolution. We release the benchmark **so others can measure directional resolution on their own screen rather than assume it**, and `evaluate_directions()` returns balanced accuracy + a permutation p-value as the honest metrics (raw accuracy is inflated by class imbalance).

## Design principles

- **Direction is a first-class output, not a magnitude by-product.**
- **Empirical, composition-matched null** — significance is judged against random signed sets matched in up/down composition, not a parametric assumption.
- **Calibrated honesty** — the calibration module makes the resolution ceiling measurable on the user's own data instead of over-promising.
- **Assay-agnostic** — anything reducible to a perturbation × readout matrix and a signed program.

## Citation

Perturb2Target (Claude Science hackathon). Source screen: Zhu, Dann, … Pritchard & Marson (2025), *bioRxiv* 10.64898/2025.12.23.696273. Directional benchmark labels curated from approved-drug mechanisms.
