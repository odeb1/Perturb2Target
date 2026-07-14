#!/usr/bin/env python
"""
submit_boltz_hosted.py — Structural druggability screen via the hosted Boltz-2 API.

No local GPU required: jobs run on Boltz's GPU backend and results are polled/downloaded.

Usage:
    export BOLTZ_API_KEY=sk_bc_ws_live_...   # from Customize -> Credentials
    pip install boltz-api rdkit pandas
    python submit_boltz_hosted.py

Inputs : manifest.csv  (gene, uniprot, primary_class, direction, peak_condition,
                        disease_anchor, ligand, smiles, sequence)
Outputs: results/structural_druggability_results.csv
         structures/<GENE>_complex.cif   (best-sample protein-ligand complex)

Metrics captured per target:
    iptm, ligand_iptm, complex_plddt, structure_confidence, ptm   (structure)
    binding_confidence, optimization_score                        (affinity head)
    struct_druggability = sqrt(ligand_iptm * binding_confidence)  (composite, [0,1])
"""
import os, time, json, urllib.request
import pandas as pd
from boltz_api import Boltz

MODEL = "boltz-2.1"
client = Boltz(api_key=os.environ["BOLTZ_API_KEY"])

def build_input(seq, smiles):
    return {
        "entities": [
            {"type": "protein",       "chain_ids": ["A"], "value": seq},
            {"type": "ligand_smiles", "chain_ids": ["L"], "value": smiles},
        ],
        "binding": {"type": "ligand_protein_binding", "binder_chain_id": "L"},
    }

def main():
    man = pd.read_csv("manifest.csv")
    os.makedirs("results", exist_ok=True)
    os.makedirs("structures", exist_ok=True)

    # --- submit all jobs (non-blocking) ---
    jobs = {}
    for _, r in man.iterrows():
        g = r["gene"]
        inp = build_input(r["sequence"], r["smiles"])
        # cost check on first (optional): client.predictions.structure_and_binding.estimate_cost(input=inp, model=MODEL)
        job = client.predictions.structure_and_binding.start(
            input=inp, model=MODEL, idempotency_key=f"p2t_{g}_v1")
        jobs[g] = job.id
        print(f"submitted {g:8s} -> {job.id}")

    # --- poll to completion ---
    pending = set(jobs)
    while pending:
        time.sleep(30)
        for g in list(pending):
            st = client.predictions.structure_and_binding.retrieve(jobs[g]).status
            if st != "running":
                print(f"  {g}: {st}")
                pending.discard(g)

    # --- harvest metrics + structures ---
    rows = []
    for _, r in man.iterrows():
        g = r["gene"]
        out = client.predictions.structure_and_binding.retrieve(jobs[g]).output.to_dict()
        m  = out["best_sample"]["metrics"]
        bm = out.get("binding_metrics", {}) or {}
        li, bc = m["ligand_iptm"], bm.get("binding_confidence", float("nan"))
        rows.append({
            "gene": g, "uniprot": r["uniprot"], "class": r["primary_class"],
            "direction": r["direction"], "peak_condition": r["peak_condition"],
            "disease_anchor": r.get("disease_anchor", ""), "ligand": r["ligand"], "smiles": r["smiles"],
            "iptm": m["iptm"], "ligand_iptm": li, "complex_plddt": m["complex_plddt"],
            "structure_confidence": m["structure_confidence"], "ptm": m["ptm"],
            "binding_confidence": bc, "optimization_score": bm.get("optimization_score", float("nan")),
            "struct_druggability": round((li * bc) ** 0.5, 4) if bc == bc else float("nan"),
            "job_id": jobs[g],
        })
        # download best-sample complex
        url = out["best_sample"]["structure"]["url"]
        data = urllib.request.urlopen(url, timeout=90).read()
        open(f"structures/{g}_complex.cif", "wb").write(data)

    pd.DataFrame(rows).sort_values("struct_druggability", ascending=False)\
        .to_csv("results/structural_druggability_results.csv", index=False)
    print("done -> results/structural_druggability_results.csv")

if __name__ == "__main__":
    main()
