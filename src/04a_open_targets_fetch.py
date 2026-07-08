#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 4 (companion) — Open Targets GraphQL fetcher.

This is the query loop that produced handoff/ot_deep.json, consumed by 04_genetics_support.py.
In the Claude Science environment these MCP calls run in the `repl` tool via host.mcp(...);
outside it, point OT_ENDPOINT at the public Open Targets GraphQL API
(https://api.platform.opentargets.org/api/v4/graphql) and use a plain HTTP POST.

Input : handoff/deep_query_genes.json  -> {gene_symbol: ensembl_id} for the top-N candidates
Output: handoff/ot_deep.json           -> {gene_symbol: <graphql `target` payload>}
"""
import json, time

OT_QUERY = """query($ensg:String!){ target(ensemblId:$ensg){
  approvedSymbol
  geneticConstraint{ constraintType score upperBin }
  tractability{ modality label value }
  associatedDiseases(page:{size:8,index:0}){ count rows{ score
     disease{ id name therapeuticAreas{ id name } } } }
  drugAndClinicalCandidates{ count rows{ maxClinicalStage
     drug{ name drugType maximumClinicalStage } } }
}}"""

def fetch_via_mcp(host, genes):
    """genes: {symbol: ensembl_id}. Returns {symbol: target_payload}.
    Run this inside the `repl` tool where `host.mcp` is available."""
    out = {}
    for sym, ensg in genes.items():
        try:
            res = host.mcp("clinical-genomics", "open_targets_graphql",
                           query=OT_QUERY, variables={"ensg": ensg})
            out[sym] = (res.get("data") or {}).get("target") or res
        except Exception as e:
            out[sym] = {"error": str(e)}
    return out

def fetch_via_http(genes, endpoint="https://api.platform.opentargets.org/api/v4/graphql"):
    """Portable fallback outside the Claude Science MCP environment."""
    import requests
    out = {}
    for sym, ensg in genes.items():
        try:
            r = requests.post(endpoint, json={"query": OT_QUERY, "variables": {"ensg": ensg}}, timeout=30)
            out[sym] = (r.json().get("data") or {}).get("target")
        except Exception as e:
            out[sym] = {"error": str(e)}
        time.sleep(0.1)
    return out

if __name__ == "__main__":
    genes = json.load(open("handoff/deep_query_genes.json"))
    result = fetch_via_http(genes)          # or fetch_via_mcp(host, genes) inside repl
    json.dump(result, open("handoff/ot_deep.json", "w"))
    print("fetched", len(result), "targets")
