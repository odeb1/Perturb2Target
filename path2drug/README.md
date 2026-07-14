# Path2Drug

*A small, reusable tool that explains **how** a perturbation's effect reaches a
target gene program, and exposes **druggable intermediate nodes** — for when the
nominated target itself is undruggable.*

Given a target gene + a signed gene signature + a STRING build, `Path2Drug`:

1. builds a confidence-weighted signalling graph from STRING high-confidence edges;
2. extracts the shortest / highest-confidence paths from the target to the
   signature (program) genes — the **deterministic mechanistic backbone**;
3. annotates each intermediate node with its **real perturbation-observed program
   effect** (from the Marson CRISPRi atlas) and its **ChEMBL druggability**;
4. optionally emits a plain-language narrative — strictly grounded in the
   extracted subnetwork.

## Why it is a tool, not a hallucination engine

Two guardrails, both enforced in code:

- **The LLM may only describe the extracted subnetwork.** The narration prompt is
  handed the exact node/edge set and is forbidden from introducing any gene,
  edge, or pathway not present. The deterministic backbone — not the narrative —
  is the ground-truth deliverable; if narration is unavailable (model refusal,
  offline), the backbone stands alone.
- **A deterministic baseline is always computed** (`path2drug_vs_baseline.csv`).
  A naive "highest-degree druggable neighbour" heuristic collapses to the generic
  hub **CD4** for almost every transcription factor; Path2Drug instead returns
  distinct, pathway-specific, direction-consistent druggable nodes (JAK1, JAK2,
  SRC, HDAC1, IFNAR1). That contrast *is* the result.

## Worked example (flagship)

`STAT4` (rank-9 driver, no small-molecule pocket) →
`STAT4 → IL12RB2 → JAK2 → IFNG`: the path exposes **JAK2** (approved JAK
inhibitors, phase 4) as a druggable intermediate whose *own* CRISPRi knockdown
lowers the program (−0.42, direction-consistent). This recovers the real-world
logic of treating STAT-driven inflammation with JAK inhibitors — derived, not
assumed.

## Enrichments (v2)

Four upgrades make the backbone harder to dismiss as STRING-topology alone:

1. **Functional edge confirmation from the atlas.** Each `target → intermediate`
   edge is tested against the perturbation data itself: does knocking down the
   target actually shift the intermediate's transcript in the readout? Every
   intermediate carries `edge_zscore_under_target_KD` and
   `edge_functionally_confirmed` (|z| ≥ 1.96). On the ten hard-TF backbones,
   2/24 testable edges pass strict significance — the flagship `STAT4 → STAT1`
   is strongly confirmed (z = 3.6) and `STAT4 → JAK2` is suggestive (z = 1.8).
   The low overall rate is itself informative: a transcriptional readout
   under-detects post-translational (kinase) edges, so this is honest
   *annotation*, not a pass/fail gate.
2. **Signed, directed regulatory paths (CollecTRI).** With a signed regulon
   network supplied, `explain()` also returns `signed_paths` whose arrows mean
   activation (+1) / repression (−1), not mere association, plus a
   `signed_direction_prediction`. Scope note, stated honestly: the whole-target
   aggregate sign does **not** beat chance for direction on this atlas (5/9,
   p = 0.50) because competing regulatory paths cancel — so signed paths are for
   per-path interpretability, and the atlas-measured direction (self-consistency,
   73%) remains the reliable directional signal.
3. **Direction-consistency and FDR gate as first-class fields.** Each
   intermediate reports `direction_consistent` (does its own knockdown move the
   program the way the target's call requires) and `fdr_significant` (its own KD
   effect at FDR < 0.05). Intermediates are ranked druggable → direction-consistent
   → effect magnitude.
4. **Quantified baseline contrast.** Path2Drug returns 5 distinct druggable
   nodes across the ten hard-TFs and is direction-consistent for 71% of them;
   the naive highest-degree baseline returns 2 distinct nodes, 7/8 of them the
   generic hub CD4 (0% pathway-specific).

## Usage

```python
from path2drug import PathToDrug
p2d = PathToDrug(
    edges_parquet   = "p2d_string_edges.parquet",          # STRING subnetwork
    perturb_parquet = "perturbation_scores_wide.parquet",  # atlas program effects
    druggable_csv   = "druggable_intermediate_nodes.csv",  # ChEMBL druggability
    signature_csv   = "program_signature.csv")             # the signed program

# deterministic (no LLM)
r = p2d.explain("STAT4", direction="antagonize", max_paths=4)
print(r["mechanism_backbone"])

# add a grounded narrative (optional)
r = p2d.explain("STAT4", direction="antagonize", narrate=True, host=host)
print(r.get("mechanism_narrative"))   # None if the model declined; backbone still valid
```

## Inputs / outputs

- **Inputs:** `p2d_string_edges.parquet` (columns `preferredName_A/B`, `score`),
  the atlas `perturbation_scores_wide.parquet` (indexed by gene, with
  `peak_signed_score`, `direction`, `peak_emp_fdr`), the ChEMBL neighbour table,
  and the signed `program_signature.csv`.
- **Outputs:** `p2d_<TF>.json` per target (paths, intermediate nodes with
  observed effect + druggability, backbone text, optional narrative);
  `path2drug_vs_baseline.csv` (the deterministic comparison);
  `path2drug_pathway.png` (backbone + comparison figure).

## Caveats

- STRING edges are **associations, not signed directed arrows** — the backbone is
  a hypothesis for the route, not a proven signed cascade. The narration is
  pinned to say so.
- The curated network carries **study bias** (well-studied hubs are
  over-represented); see the project's network-propagation confound analysis.
  Path2Drug mitigates this by requiring the intermediate to lie on a path to
  the program *and* to move the program itself in the atlas — not merely to be a
  hub — but it does not eliminate the bias.
- A proposed intermediate drug is a **repurposing hypothesis**, not evidence that
  drugging it phenocopies silencing the TF.
