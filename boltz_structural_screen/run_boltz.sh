#!/usr/bin/env bash
# ==========================================================================
# run_boltz.sh -- Structural druggability screen for 11 novel small-molecule
# drug-target nominations (Perturb2Target). Runs Boltz-2 co-folding + affinity
# head for each target-ligand pair, then collects confidence + affinity scores.
#
# USAGE (on a GPU machine, >= 1 GPU with >= 24 GB VRAM recommended):
#   1. pip install boltz                      # MIT, from PyPI
#   2. export BOLTZ_API_KEY=<your key>        # DO NOT hardcode; set at runtime
#      # (only needed if your Boltz build uses the hosted MSA/inference service;
#      #  the open-source `boltz` CLI + --use_msa_server does NOT require a key)
#   3. bash run_boltz.sh /path/to/boltz_inputs /path/to/out
# ==========================================================================
set -euo pipefail

IN_DIR="${1:-boltz_inputs}"
OUT_DIR="${2:-boltz_out}"
mkdir -p "$OUT_DIR"

# Recycling/samples: 3/5 is the skill's recommended default for binder validation.
RECYCLE=3
SAMPLES=5

for yaml in "$IN_DIR"/*.yaml; do
    gene=$(basename "$yaml" .yaml)
    echo "=== Boltz-2: $gene ==="
    boltz predict "$yaml" \
        --use_msa_server \
        --out_dir "$OUT_DIR/$gene" \
        --recycling_steps $RECYCLE \
        --diffusion_samples $SAMPLES \
        --output_format pdb \
        || echo "WARN: $gene failed, continuing"
done

echo "All runs done. Collect scores with: python collect_boltz_results.py $OUT_DIR"
