#!/usr/bin/env python
"""collect_boltz_results.py -- gather Boltz-2 confidence + affinity JSONs into one CSV.
Usage: python collect_boltz_results.py <boltz_out_dir>
Reads each target's confidence_*_model_0.json (iptm/ptm/plddt) and
affinity_*.json (affinity_pred_value=log10 IC50 uM, affinity_probability_binary),
merges with manifest.csv, writes structural_druggability_results.csv."""
import sys, os, json, glob, csv

out_dir = sys.argv[1] if len(sys.argv)>1 else "boltz_out"
man_path = os.path.join(os.path.dirname(out_dir.rstrip("/")) or ".", "manifest.csv")
if not os.path.exists(man_path):
    man_path = "boltz_inputs/manifest.csv"
manifest = {r["gene"]: r for r in csv.DictReader(open(man_path))} if os.path.exists(man_path) else {}

rows=[]
for gene_dir in sorted(glob.glob(os.path.join(out_dir,"*"))):
    gene=os.path.basename(gene_dir)
    # find prediction subdir
    conf = glob.glob(os.path.join(gene_dir,"**","confidence_*_model_0.json"), recursive=True)
    aff  = glob.glob(os.path.join(gene_dir,"**","affinity_*.json"), recursive=True)
    row={"gene":gene}
    row.update({k:manifest.get(gene,{}).get(k,"") for k in
                ("direction","primary_class","peak_condition","disease","ligand","integrated_score")})
    if conf:
        c=json.load(open(conf[0]))
        row["iptm"]=c.get("iptm"); row["ptm"]=c.get("ptm")
        row["complex_plddt"]=c.get("complex_plddt"); row["confidence_score"]=c.get("confidence_score")
    if aff:
        a=json.load(open(aff[0]))
        row["affinity_pred_value_log10IC50uM"]=a.get("affinity_pred_value")
        row["affinity_probability_binary"]=a.get("affinity_probability_binary")
    rows.append(row)

if rows:
    cols=["gene","direction","primary_class","peak_condition","disease","ligand","integrated_score",
          "iptm","ptm","complex_plddt","confidence_score",
          "affinity_pred_value_log10IC50uM","affinity_probability_binary"]
    with open("structural_druggability_results.csv","w",newline="") as f:
        w=csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({k:r.get(k,"") for k in cols})
    print(f"wrote structural_druggability_results.csv ({len(rows)} targets)")
else:
    print("no results found under", out_dir)
