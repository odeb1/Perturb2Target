# Structural druggability screen (Boltz-2) — Perturb2Target extension

Physical-evidence layer for the top novel small-molecule drug-target nominations.
For each target we co-fold the human protein with a class-appropriate drug-like
ligand and run Boltz-2's affinity head, upgrading "druggable by protein class"
to "druggable by predicted structure + predicted affinity."

## Candidate selection
From the 86 novel-priority targets, we kept the 11 with a genuine small-molecule
pocket (classes: kinase, enzyme, GPCR). Transcription factors and transporters
were excluded from the affinity analysis (no classic catalytic pocket / poor
membrane folding), though they remain in the shortlist as hard-modality targets.

## Files
- `<GENE>.yaml`         — Boltz-2 input: human protein chain A + ligand chain L + affinity head
- `manifest.csv/json`   — target metadata + ligand SMILES + rationale
- `run_boltz.sh`        — GPU run script (reads BOLTZ_API_KEY from env; never hardcoded)
- `collect_boltz_results.py` — gathers confidence + affinity JSONs into one CSV

## How to run (on your GPU machine)
```bash
pip install boltz
export BOLTZ_API_KEY=<your key>          # set at runtime; not stored in these files
bash run_boltz.sh boltz_inputs boltz_out
python collect_boltz_results.py boltz_out
```
Produces `structural_druggability_results.csv`. Send that back and I'll integrate
it into the shortlist + build the structural-druggability figure.

## Reading the output (per Boltz-2 conventions)
- `iptm` > 0.5  — interface (pocket-ligand) confidence pass line
- `complex_plddt` > 0.7 — fold confidence
- `affinity_probability_binary` (0-1) — binder vs non-binder; **rank hits by this**
- `affinity_pred_value` — log10(IC50 in uM); lower = tighter (0 -> 1 uM, -3 -> 1 nM)

## Ligand rationale
Each ligand is a known inhibitor of the target's nearest drugged homolog or a
validated chemical probe — a concrete scaffold-repurposing starting point, not a
random binder. See `rationale` in manifest.json.

## Compute notes
- ~24 GB VRAM recommended; if OOM, lower --diffusion_samples or --max_parallel_samples.
- --use_msa_server queries api.colabfold.com (30-90 s/chain, CPU-side).
- If cuequivariance kernels are missing, add --no_kernels (2x slower, identical results).
- All 11 ligands are 24-31 heavy atoms, well under Boltz's 128-atom affinity cap.
