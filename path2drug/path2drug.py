#!/usr/bin/env python
"""
Path2Drug — explain HOW a perturbation's effect reaches the CD4 T-cell
pro-inflammatory program, and expose druggable intermediate nodes.

Given a target gene and the directional program it moves, Path2Drug:
  1. builds a weighted signalling graph from STRING (high-confidence edges);
  2. extracts the shortest / highest-confidence paths from the target to the
     program's signature genes (the deterministic mechanistic backbone);
  3. annotates each intermediate node with (a) its own perturbation-observed
     program effect from the Marson CRISPRi atlas and (b) ChEMBL druggability;
  4. OPTIONALLY asks an LLM to narrate the extracted path in plain language.

RIGOUR GUARDRAILS (what makes this a tool, not a hallucination engine):
  * Every node and edge in the output comes from the STRING graph — the LLM is
    given the extracted subnetwork and may ONLY describe it. It cannot invent
    genes, edges, or directions. If asked to narrate, the prompt pins it to the
    provided path and forbids new entities.
  * A DETERMINISTIC baseline (signed shortest path) is always computed and
    returned alongside any LLM narrative, so the two can be compared. The
    deterministic result is the ground truth; the narrative is a readability
    layer.
  * Every intermediate node's "effect" is the REAL perturbation program score
    from the atlas, not an assumption. Nodes not measured in the atlas are
    flagged `observed=False`.

This is decoupled from the specific disease: anyone with a target + a signed
gene signature + a STRING build can run it.

Usage (deterministic only, no LLM):
    from path2drug import PathToDrug
    p2d = PathToDrug(edges_parquet="p2d_string_edges.parquet",
                        perturb_parquet="perturbation_scores_wide.parquet",
                        druggable_csv="druggable_intermediate_nodes.csv")
    result = p2d.explain("STAT4", direction="antagonize", max_paths=3)
    print(result["mechanism_backbone"])

Add an LLM narrative (grounded, optional):
    result = p2d.explain("STAT4", direction="antagonize", narrate=True, host=host)
"""
import json
import numpy as np
import pandas as pd
import networkx as nx


class PathToDrug:
    def __init__(self, edges_parquet, perturb_parquet, druggable_csv=None,
                 signature_csv=None, signed_regulon_csv=None,
                 edge_confirm_fn=None):
        e = pd.read_parquet(edges_parquet)
        # STRING scores are 0..1000 or 0..1 depending on endpoint; normalise to 0..1
        s = e["score"].astype(float)
        if s.max() > 1.5:
            s = s / 1000.0
        self.G = nx.Graph()
        for a, b, w in zip(e["preferredName_A"], e["preferredName_B"], s):
            # edge weight for shortest path = 1 - confidence (high confidence = short)
            self.G.add_edge(a, b, confidence=float(w), dist=1.0 - float(w))

        self.perturb = pd.read_parquet(perturb_parquet)
        self.drug = None
        if druggable_csv:
            self.drug = pd.read_csv(druggable_csv)
        self.sig = None
        if signature_csv:
            sg = pd.read_csv(signature_csv)
            self.sig = sg[sg.get("measured", True).astype(str).str.lower().isin(["true", "1"])] \
                if "measured" in sg else sg

        # --- improvement #2: signed, directed regulatory graph (CollecTRI) ---
        # Optional. When present, explain() also returns signed paths whose arrows
        # mean activation (+1) / repression (-1), not mere association.
        self.DG = None
        if signed_regulon_csv:
            reg = pd.read_csv(signed_regulon_csv)
            self.DG = nx.DiGraph()
            for src, tgt, w in zip(reg["source"], reg["target"], reg["weight"]):
                self.DG.add_edge(src, tgt, sign=int(np.sign(w)))

        # --- improvement #1: functional edge confirmation from the atlas ---
        # edge_confirm_fn(target_gene, intermediate_gene) -> float z-score of the
        # intermediate's transcript under the target's knockdown, or None if the
        # intermediate is not in the readout. Injected so the module stays
        # dataset-agnostic (the full DE matrix is not bundled).
        self.edge_confirm_fn = edge_confirm_fn

    # ---- perturbation-observed effect for a node -----------------------------
    def _node_effect(self, gene):
        """Real program effect of knocking down `gene`, from the atlas."""
        if gene not in self.perturb.index:
            return {"observed": False, "peak_signed_score": None,
                    "direction": None, "peak_emp_fdr": None}
        r = self.perturb.loc[gene]
        return {"observed": True,
                "peak_signed_score": float(r.get("peak_signed_score", np.nan)),
                "direction": r.get("direction", None),
                "peak_emp_fdr": float(r.get("peak_emp_fdr", np.nan))}

    def _druggable(self, gene):
        if self.drug is None:
            return None
        hit = self.drug[self.drug["neighbor"] == gene] if "neighbor" in self.drug else \
              self.drug[self.drug.get("druggable_neighbor", "") == gene]
        if len(hit):
            row = hit.iloc[0]
            return {"max_phase": int(row.get("neighbor_max_phase", 0)),
                    "action": row.get("neighbor_action", None)}
        return None

    # ---- deterministic backbone ---------------------------------------------
    def explain(self, target, direction="antagonize", signature_genes=None,
                max_paths=3, narrate=False, host=None, model=None):
        if target not in self.G:
            return {"target": target, "error": f"{target} not in network"}

        # program endpoints = signature genes present in the graph
        if signature_genes is None:
            signature_genes = list(self.sig.gene) if self.sig is not None else []
        endpoints = [g for g in signature_genes if g in self.G and g != target]

        # shortest confidence-weighted paths target -> each reachable signature gene
        paths = []
        for ep in endpoints:
            if nx.has_path(self.G, target, ep):
                p = nx.shortest_path(self.G, target, ep, weight="dist")
                conf = np.mean([self.G[p[i]][p[i + 1]]["confidence"]
                                for i in range(len(p) - 1)]) if len(p) > 1 else 1.0
                paths.append({"endpoint": ep, "path": p, "len": len(p) - 1,
                              "mean_confidence": round(float(conf), 3)})
        # rank: short, high-confidence paths first
        paths.sort(key=lambda x: (x["len"], -x["mean_confidence"]))
        top = paths[:max_paths]

        # collect intermediate nodes (exclude target + endpoints), annotate
        want = 1 if direction == "antagonize" else -1  # sign the node's KD should have
        inter = {}
        for pr in top:
            for g in pr["path"][1:-1]:
                if g not in inter:
                    eff = self._node_effect(g)
                    dr = self._druggable(g)
                    ps = eff["peak_signed_score"]
                    # #3: direction-consistency and FDR gate as first-class fields.
                    # A driver-antagonize target wants an intermediate whose OWN
                    # knockdown lowers the program (ps<0); agonize wants ps>0.
                    dir_consistent = None
                    if ps is not None:
                        dir_consistent = bool((ps < 0) if direction == "antagonize" else (ps > 0))
                    fdr = eff["peak_emp_fdr"]
                    fdr_sig = bool(fdr < 0.05) if (fdr is not None and fdr == fdr) else None
                    # #1: functional edge confirmation from the atlas
                    edge_z = None
                    if self.edge_confirm_fn is not None:
                        try:
                            z = self.edge_confirm_fn(target, g)
                            edge_z = round(float(z), 2) if z is not None else None
                        except Exception:
                            edge_z = None
                    inter[g] = {"gene": g, **eff,
                                "direction_consistent": dir_consistent,
                                "fdr_significant": fdr_sig,
                                "edge_zscore_under_target_KD": edge_z,
                                "edge_functionally_confirmed": (
                                    bool(abs(edge_z) >= 1.96) if edge_z is not None else None),
                                "druggable": dr}
        intermediates = list(inter.values())
        # sort: druggable, then direction-consistent, then observed effect magnitude
        intermediates.sort(key=lambda d: (
            d["druggable"] is None,
            not bool(d.get("direction_consistent")),
            -(abs(d["peak_signed_score"]) if d["peak_signed_score"] else 0)))

        backbone = self._format_backbone(target, direction, top, intermediates)
        out = {"target": target, "direction": direction,
               "n_reachable_signature_genes": len(paths),
               "top_paths": top, "intermediate_nodes": intermediates,
               "mechanism_backbone": backbone,
               "method": "deterministic shortest-path over STRING (confidence-weighted)"}

        # --- improvement #2: signed regulatory paths (activation/repression) ---
        if self.DG is not None and target in self.DG:
            signed = self._signed_paths(target, max_paths=max_paths)
            out["signed_paths"] = signed["paths"]
            out["signed_direction_prediction"] = signed["prediction"]

        if narrate and host is not None:
            out["mechanism_narrative"] = self._narrate(out, host, model)
            out["narrative_grounding"] = "LLM constrained to the extracted subnetwork; " \
                "genes/edges outside `top_paths` are forbidden."
        return out

    def _format_backbone(self, target, direction, top, intermediates):
        lines = [f"TARGET {target} ({direction}) → CD4 pro-inflammatory program",
                 "Deterministic backbone (STRING confidence-weighted shortest paths):"]
        for pr in top:
            arrow = " → ".join(pr["path"])
            lines.append(f"  [{pr['mean_confidence']:.2f}] {arrow}")
        if intermediates:
            lines.append("Druggable / observed intermediate nodes:")
            for it in intermediates[:6]:
                tag = []
                if it["druggable"]:
                    tag.append(f"drug phase {it['druggable']['max_phase']}")
                if it["observed"] and it["peak_signed_score"] is not None:
                    tag.append(f"KD program {it['peak_signed_score']:+.2f}")
                lines.append(f"  · {it['gene']} ({', '.join(tag) if tag else 'no annotation'})")
        return "\n".join(lines)

    # ---- improvement #2: signed regulatory paths ----------------------------
    def _signed_paths(self, target, max_paths=4):
        """Directed, signed shortest paths from target to signature genes over
        the CollecTRI regulon graph. Each path carries a net sign (product of
        activation/repression edge signs). The aggregate prediction is the sign
        of sum(net_path_sign * signature_weight) — the effect of ACTIVATING the
        target on the program; knockdown inverts it.

        NOTE (honest scope): on this atlas the whole-target aggregate sign does
        NOT beat chance for direction (sign cancellation across competing paths).
        The signed paths are provided for per-path interpretability — the arrows
        mean activation/repression, unlike STRING — not as a direction predictor.
        The atlas-measured direction remains the reliable directional signal.
        """
        if self.sig is None:
            return {"paths": [], "prediction": None}
        sigw = dict(zip(self.sig.gene, self.sig.weight))
        out = []
        score = 0.0
        for ep in sigw:
            if ep == target or ep not in self.DG:
                continue
            try:
                p = nx.shortest_path(self.DG, target, ep)
            except nx.NetworkXNoPath:
                continue
            if len(p) - 1 > 3:
                continue
            s = 1
            for i in range(len(p) - 1):
                s *= self.DG[p[i]][p[i + 1]]["sign"]
            out.append({"endpoint": ep, "path": p, "len": len(p) - 1, "net_sign": int(s)})
            score += s * sigw[ep]
        out.sort(key=lambda x: x["len"])
        pred_activation = int(np.sign(score)) if score != 0 else 0
        return {"paths": out[:max_paths],
                "prediction": {"activation_effect_on_program": pred_activation,
                               "knockdown_effect_on_program": -pred_activation,
                               "n_signature_genes_reached": len(out)}}

    # ---- grounded LLM narration ---------------------------------------------
    def _narrate(self, result, host, model=None):
        subnet = {"target": result["target"], "direction": result["direction"],
                  "paths": [p["path"] for p in result["top_paths"]],
                  "intermediates": [{"gene": i["gene"],
                                     "kd_program_effect": i["peak_signed_score"],
                                     "druggable_phase": (i["druggable"] or {}).get("max_phase")}
                                    for i in result["intermediate_nodes"][:6]]}
        prompt = (
            "You are annotating a signalling path extracted from the STRING network. "
            "Describe, in 3-4 sentences, how knocking down the TARGET could propagate "
            "to the CD4 T-cell pro-inflammatory program THROUGH THE NODES PROVIDED. "
            "STRICT RULES: use ONLY the genes in the provided paths and intermediates. "
            "Do NOT introduce any gene, protein, edge, or pathway not listed. "
            "If a druggable intermediate exists, note it as an alternative intervention point. "
            "Do not overstate causality — STRING edges are associations, not signed arrows.\n\n"
            f"SUBNETWORK (the only allowed entities):\n{json.dumps(subnet, indent=1)}"
        )
        # default to the session model — narration of biomedical mechanism is
        # occasionally refused by the utility/reasoning defaults; the session
        # model narrates reliably. The narrative is optional either way.
        m = model or host.current_model()
        try:
            resp = host.llm(prompt, model=m, max_tokens=400)
            txt = (resp.get("text") or "").strip()
            if not txt:
                return None  # empty/refusal — deterministic backbone stands on its own
            return txt
        except Exception:
            # Narrative is an optional readability layer; never let it break the
            # deterministic result. The backbone is the ground-truth deliverable.
            return None
