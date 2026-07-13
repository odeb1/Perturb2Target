"""
Perturb2Target — interactive directional target explorer
=========================================================

A self-contained, visually rich Streamlit app over the 1,923-target directional
shortlist from genome-scale CD4+ T-cell CRISPRi Perturb-seq.

Layers:
  1. DETERMINISTIC core (always on, no network): sidebar filters + a natural-language
     box that maps to structured filters. Filtering runs entirely on the real CSV.
  2. OPTIONAL Claude layer (only if ANTHROPIC_API_KEY is set): (a) parses a free-form
     query into the SAME structured filter spec via a tool schema, and (b) writes a
     plain-language explanation of a selected gene STRICTLY from that gene's computed
     evidence row. The model never invents a gene or a number.
  3. VISUALS: KPI cards, an interactive "directional map" (Plotly), a per-gene evidence
     bar, and a live 3D AlphaFold structure (py3Dmol, fetched on demand from EBI).

Run:
    pip install -r requirements.txt
    export ANTHROPIC_API_KEY=sk-...      # optional; app works fully without it
    streamlit run app.py
"""
import os
import json
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "target_shortlist.csv")

# brand palette
C_DRIVER = "#c1272d"   # driver_antagonize (block) — red
C_BRAKE = "#1a7f9e"    # brake_agonize (activate) — teal
C_BG = "#f5f7fa"       # app background (light)
C_PANEL = "#ffffff"    # plot / card background
C_INK = "#1b2733"      # primary text
C_STRUCT_BG = "0xeef2f7"  # 3D viewer background (light)

# ----------------------------------------------------------------------------- data
@st.cache_data
def load_shortlist():
    df = pd.read_csv(CSV)
    return df

@st.cache_data
def load_signature_response():
    """Per-signature-gene z-scores per nomination per condition (49 demo genes × 48 sig × 3)."""
    fp = os.path.join(HERE, "signature_response.parquet")
    return pd.read_parquet(fp) if os.path.exists(fp) else None

@st.cache_data
def load_map_by_condition():
    """Per-condition directional scores for the whole shortlist (for the condition slider)."""
    fp = os.path.join(HERE, "map_by_condition.parquet")
    return pd.read_parquet(fp) if os.path.exists(fp) else None

@st.cache_data
def load_ms_nominations():
    fp = os.path.join(HERE, "ms_nominations_shortlist.csv")
    return pd.read_csv(fp) if os.path.exists(fp) else None

@st.cache_data
def load_ms_concordance():
    fp = os.path.join(HERE, "ms_drug_concordance.csv")
    return pd.read_csv(fp) if os.path.exists(fp) else None

@st.cache_data
def load_ms_summary():
    fp = os.path.join(HERE, "ms_generalization_summary.json")
    return json.load(open(fp)) if os.path.exists(fp) else None

LIGAND_DIR = os.path.join(HERE, "structures_ligand")
def ligand_genes():
    if os.path.isdir(LIGAND_DIR):
        return sorted(f[:-4] for f in os.listdir(LIGAND_DIR) if f.endswith(".cif"))
    return []

# ------------------------------------------------------------- external structure
STRUCT_DIR = os.path.join(HERE, "structures")

def _http_get(url, timeout=30):
    """stdlib urllib fetch — avoids the requests/urllib3+macOS-SSL segfault seen in some
    miniforge setups. Returns text or None."""
    import urllib.request
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

@st.cache_data(show_spinner=False)
def gene_to_uniprot(sym):
    try:
        txt = _http_get("https://rest.uniprot.org/uniprotkb/search?"
                        + "query=gene_exact:%s+AND+organism_id:9606+AND+reviewed:true" % sym
                        + "&fields=accession&format=tsv&size=1")
        lines = txt.strip().split("\n")
        return lines[1] if len(lines) > 1 else None
    except Exception:
        return None

@st.cache_data(show_spinner=False)
def load_structure(gene):
    """Return (pdb_text, source_label). Prefers a bundled local structure (offline-safe,
    no native-SSL crash); falls back to a live AlphaFold fetch via stdlib urllib."""
    import json as _json
    local = os.path.join(STRUCT_DIR, f"{gene}.pdb")
    if os.path.exists(local) and os.path.getsize(local) > 1000:
        with open(local) as fh:
            return fh.read(), f"bundled · {gene}.pdb"
    # network fallback (stdlib only)
    try:
        acc = gene_to_uniprot(gene)
        if not acc:
            return None, None
        meta = _json.loads(_http_get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=20))
        if not meta:
            return None, None
        pdb = _http_get(meta[0]["pdbUrl"], timeout=40)
        return pdb, f"AlphaFold DB · {acc}"
    except Exception:
        return None, None

def structure_html(pdb, style="confidence", height=420, spin=True):
    import py3Dmol
    view = py3Dmol.view(width=440, height=height)
    view.addModel(pdb, "pdb")
    if style == "confidence":
        # AlphaFold stores per-residue pLDDT in the B-factor column
        view.setStyle({"cartoon": {"colorscheme": {"prop": "b",
                       "gradient": "roygb", "min": 50, "max": 90}}})
    else:
        view.setStyle({"cartoon": {"color": "spectrum"}})
    view.setBackgroundColor(C_STRUCT_BG)
    view.zoomTo()
    if spin:
        view.spin("y", 0.4)   # gentle auto-rotation on load (0.4 = slow)
    return view._make_html()

def ligand_structure_html(gene, height=440, spin=True, staged=True, reveal_ms=900):
    """Render a Boltz co-folded complex: protein cartoon + the docked ligand.

    staged=True → a genuine timed reveal: the protein cartoon renders first (spinning), then
    after `reveal_ms` a setTimeout injected into the 3Dmol.js callback makes the pre-computed
    ligand appear (green sticks + spheres) while the camera animates a 1.2 s zoom into the
    pocket. This shows/hides the Boltz co-folded pose on a timer — it does NOT invent any
    binding motion; the ligand does not move, it becomes visible in its predicted position.
    staged=False → protein and ligand both visible from the first frame (no timed reveal).
    """
    import py3Dmol, re, json
    cif = open(os.path.join(LIGAND_DIR, f"{gene}.cif")).read()
    view = py3Dmol.view(width=440, height=height)
    view.addModel(cif, "cif")
    view.setStyle({"cartoon": {"color": "#7fa8c9"}})   # protein cartoon
    if not staged:
        view.addStyle({"hetflag": True}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.22}})
        view.addStyle({"hetflag": True}, {"sphere": {"scale": 0.28, "colorscheme": "greenCarbon"}})
    view.setBackgroundColor(C_STRUCT_BG)
    view.zoomTo()                                       # frame the whole protein first
    if not staged:
        try:
            view.zoomTo({"hetflag": True}, 1200)
        except Exception:
            view.zoomTo({"hetflag": True})
    if spin:
        view.spin("y", 0.4)
    view.render()
    html = view._make_html()
    if not staged:
        return html
    # --- inject a real timed reveal of the ligand into the 3Dmol.js .then() callback ---
    m = re.search(r'(viewer_\d+)', html)
    if not m:
        return html
    vv = m.group(1)
    stick = json.dumps({"stick": {"colorscheme": "greenCarbon", "radius": 0.22}})
    sphere = json.dumps({"sphere": {"scale": 0.28, "colorscheme": "greenCarbon"}})
    inject = (f'\nsetTimeout(function(){{'
              f'{vv}.addStyle({{"hetflag":true}},{stick});'
              f'{vv}.addStyle({{"hetflag":true}},{sphere});'
              f'{vv}.zoomTo({{"hetflag":true}},1200);'
              f'{vv}.render();}},{reveal_ms});\n')
    idx = html.rfind("});")   # closing brace of the .then(function(viewer){...}) callback
    if idx != -1:
        html = html[:idx] + inject + html[idx:]
    return html

# ------------------------------------------------------------------ filter engine
# The single canonical filter spec that BOTH the sidebar and the NL parser emit.
def apply_filters(df, spec):
    d = df.copy()
    if spec.get("direction") in ("driver_antagonize", "brake_agonize"):
        d = d[d.direction == spec["direction"]]
    if spec.get("primary_class"):
        classes = spec["primary_class"] if isinstance(spec["primary_class"], list) else [spec["primary_class"]]
        d = d[d.primary_class.isin(classes)]
    if spec.get("genetics_tier"):
        tiers = spec["genetics_tier"] if isinstance(spec["genetics_tier"], list) else [spec["genetics_tier"]]
        d = d[d.genetics_tier.isin(tiers)]
    if spec.get("druggable_only"):
        d = d[d.is_druggable == True]  # noqa: E712
    if spec.get("novel_only"):
        d = d[d.novelty_class == "novel_undrugged"]
    if spec.get("novelty_class"):
        nc = spec["novelty_class"] if isinstance(spec["novelty_class"], list) else [spec["novelty_class"]]
        d = d[d.novelty_class.isin(nc)]
    if spec.get("disease_contains"):
        q = str(spec["disease_contains"]).lower()
        d = d[d.ot_top_disease.fillna("").str.lower().str.contains(q)
              | d.ot_autoimmune_disease.fillna("").str.lower().str.contains(q)]
    if spec.get("condition") in ("Rest", "Stim8hr", "Stim48hr"):
        d = d[d.peak_condition == spec["condition"]]
    if spec.get("min_integrated_score") is not None:
        d = d[d.integrated_score >= float(spec["min_integrated_score"])]
    if spec.get("max_rank") is not None:
        d = d[d["rank"] <= int(spec["max_rank"])]
    return d.sort_values("integrated_score", ascending=False)

# ------------------------------------------------------------------- Claude layer
FILTER_TOOL = {
    "name": "set_filters",
    "description": "Set the shortlist filter spec from the user's query. Omit any field you cannot infer.",
    "input_schema": {
        "type": "object",
        "properties": {
            "direction": {"type": "string", "enum": ["driver_antagonize", "brake_agonize"]},
            "primary_class": {"type": "array", "items": {"type": "string", "enum": [
                "transcription_factor", "enzyme", "transporter", "kinase", "catalytic_receptor",
                "GPCR", "cytokine", "ion_channel", "nuclear_receptor", "cytokine_receptor"]}},
            "genetics_tier": {"type": "array", "items": {"type": "string", "enum": ["strong", "moderate", "weak"]}},
            "condition": {"type": "string", "enum": ["Rest", "Stim8hr", "Stim48hr"]},
            "druggable_only": {"type": "boolean"},
            "novel_only": {"type": "boolean"},
            "disease_contains": {"type": "string", "description": "a SPECIFIC disease name only (asthma, rheumatoid arthritis, type 1 diabetes, lupus, IBD, psoriasis). NOT generic words like 'inflammation' or 'immune'."},
            "max_rank": {"type": "integer"},
            "min_integrated_score": {"type": "number"},
        },
    },
}
PARSE_SYS = (
    "Translate the user's natural-language query about a drug-target shortlist into a call to set_filters. "
    "'antagonize'/'driver'/'block or inhibit a driver' => direction=driver_antagonize. "
    "'agonize'/'brake'/'activate a brake' => direction=brake_agonize. "
    "'novel'/'undrugged'/'clean patent space' => novel_only=true. "
    "'druggable'/'small-molecule'/'antibody' => druggable_only=true. 'top N' => max_rank=N. "
    "Only set disease_contains for a specific named disease. Never invent gene names."
)

def get_client():
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        return None
    try:
        import anthropic
        return anthropic.Anthropic(api_key=key)
    except Exception:
        return None

MODEL = os.environ.get("ANTHROPIC_MODEL", "claude-3-5-haiku-latest")

def nl_to_spec(client, query):
    msg = client.messages.create(
        model=MODEL, max_tokens=400, system=PARSE_SYS,
        tools=[FILTER_TOOL], tool_choice={"type": "tool", "name": "set_filters"},
        messages=[{"role": "user", "content": query}],
    )
    for block in msg.content:
        if block.type == "tool_use":
            return dict(block.input)
    return {}

def explain_gene(client, row):
    """Explanation generated STRICTLY from the gene's own evidence row."""
    ev = {
        "gene": row.target_gene, "rank": int(row["rank"]),
        "direction": row.direction, "integrated_score": round(float(row.integrated_score), 3),
        "peak_condition": row.peak_condition,
        "peak_emp_z": round(float(row.peak_emp_z), 2), "peak_emp_fdr": float(row.peak_emp_fdr),
        "genetics_tier": row.genetics_tier, "top_disease": row.ot_top_disease,
        "primary_class": row.primary_class, "suggested_modality": row.suggested_modality,
        "novelty_class": row.novelty_class, "n_drugs": int(row.n_drugs) if pd.notna(row.n_drugs) else 0,
    }
    sys = (
        "You explain ONE drug-target nomination to an immunologist, using ONLY the evidence JSON provided. "
        "Do not add facts, numbers, diseases, or mechanisms not present in the JSON. 3-4 sentences. "
        "Explain what the directional call means (driver_antagonize = knockdown lowers the inflammatory "
        "program, so BLOCK it; brake_agonize = knockdown raises it, so ACTIVATE it), and tie together the "
        "genetics tier, disease anchor, and suggested modality. If a value is missing, say so; never guess."
    )
    msg = client.messages.create(
        model=MODEL, max_tokens=350, system=sys,
        messages=[{"role": "user", "content": "Evidence JSON:\n" + json.dumps(ev, indent=1)}],
    )
    return "".join(b.text for b in msg.content if b.type == "text")

# -------------------------------------------------------------------- plots
def directional_map(dfp, highlight=None):
    """Hero scatter: causal strength (x) vs genetics (y), size=integrated, color=direction."""
    import plotly.graph_objects as go
    fig = go.Figure()
    for dirval, col, name in [("driver_antagonize", C_DRIVER, "Driver → antagonize (block)"),
                              ("brake_agonize", C_BRAKE, "Brake → agonize (activate)")]:
        d = dfp[dfp.direction == dirval]
        fig.add_trace(go.Scatter(
            x=d.causal_component, y=d.genetics_component, mode="markers",
            name=name,
            marker=dict(size=6 + 22 * d.integrated_score.clip(0, 1),
                        color=col, opacity=0.72, line=dict(width=0.5, color="white")),
            text=d.target_gene,
            customdata=d[["rank", "integrated_score", "genetics_tier", "ot_top_disease"]].values,
            hovertemplate="<b>%{text}</b><br>rank %{customdata[0]}<br>"
                          "integrated %{customdata[1]:.3f}<br>causal %{x:.2f} · genetics %{y:.2f}"
                          "<br>%{customdata[2]} genetics · %{customdata[3]}<extra></extra>"))
    if highlight is not None and len(highlight):
        fig.add_trace(go.Scatter(
            x=highlight.causal_component, y=highlight.genetics_component, mode="markers+text",
            name="selection", text=highlight.target_gene, textposition="top center",
            textfont=dict(color=C_INK, size=12, family="Arial Black"),
            marker=dict(size=16, color="rgba(0,0,0,0)", line=dict(width=2.6, color="#e08a1e")),
            hoverinfo="skip", showlegend=False))
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK),
        height=500, margin=dict(l=60, r=20, t=70, b=90),
        title=dict(text="Directional map — causal effect vs human genetics<br>"
                        "<sup>bubble size = integrated score</sup>",
                   font=dict(size=15, color=C_INK), x=0.01, xanchor="left", y=0.97, yanchor="top"),
        xaxis_title="causal directional strength", yaxis_title="human-genetics support",
        xaxis=dict(gridcolor="#e3e8ef"), yaxis=dict(gridcolor="#e3e8ef"),
        legend=dict(orientation="h", yanchor="top", y=-0.16, x=0.5, xanchor="center",
                    bgcolor="rgba(0,0,0,0)"))
    return fig

def evidence_bar(row):
    import plotly.graph_objects as go
    comps = {"Causal": row.causal_component, "Genetics": row.genetics_component,
             "Druggability": row.drug_component, "Novelty": row.novelty_component}
    fig = go.Figure(go.Bar(
        x=list(comps.values()), y=list(comps.keys()), orientation="h",
        marker=dict(color=[C_DRIVER if row.direction == "driver_antagonize" else C_BRAKE] * 4,
                    line=dict(width=0)),
        text=[f"{v:.2f}" for v in comps.values()], textposition="outside"))
    fig.update_layout(template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
                      font=dict(color=C_INK),
                      height=230, margin=dict(l=10, r=30, t=36, b=10),
                      title=dict(text=f"{row.target_gene}: evidence components", font=dict(size=13, color=C_INK)),
                      xaxis=dict(range=[0, 1.05], gridcolor="#e3e8ef"))
    return fig

def signature_heatmap(sigdf, gene, condition):
    """#1 Mechanism reveal: 48 signature genes' z-scores under this KD, split by arm."""
    import plotly.graph_objects as go
    d = sigdf[(sigdf.target_gene == gene) & (sigdf.condition == condition)].copy()
    if not len(d):
        return None
    # order: pro-inflammatory block then regulatory block, each sorted by z
    d["arm_order"] = (d.arm == "regulatory").astype(int)
    d = d.sort_values(["arm_order", "zscore"], ascending=[True, False])
    colors = d.zscore.tolist()
    fig = go.Figure(go.Bar(
        x=d.zscore, y=d.signature_gene, orientation="h",
        marker=dict(color=colors, colorscale="RdBu_r", cmid=0, cmin=-6, cmax=6,
                    colorbar=dict(title="z", thickness=12, len=0.6)),
        hovertemplate="<b>%{y}</b><br>z = %{x:.2f}<extra></extra>"))
    # arm separator annotation
    n_pro = int((d.arm == "pro_inflammatory").sum())
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK), height=760, margin=dict(l=10, r=10, t=48, b=10),
        title=dict(text=f"{gene} KD @ {condition}: signature-gene response<br>"
                        f"<sup>top {n_pro} = pro-inflammatory arm · bottom = regulatory arm</sup>",
                   font=dict(size=13, color=C_INK)),
        xaxis=dict(title="per-gene z-score (KD vs control)", gridcolor="#e3e8ef", zeroline=True,
                   zerolinecolor="#9fb3c8"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=9)))
    return fig

def signature_heatmap_animated(sigdf, gene):
    """#3 Animated mechanism: signature-gene response replayed Rest → Stim8hr → Stim48hr.
    Gene order is fixed (from the peak condition) so bars stay in place and only lengths move."""
    import plotly.graph_objects as go
    conds = ["Rest", "Stim8hr", "Stim48hr"]
    sub = sigdf[sigdf.target_gene == gene]
    if not len(sub):
        return None
    # fixed row order from the Stim48hr frame (arm block, then z within arm)
    order_src = sub[sub.condition == "Stim48hr"].copy()
    if not len(order_src):
        order_src = sub[sub.condition == conds[0]].copy()
    order_src["arm_order"] = (order_src.arm == "regulatory").astype(int)
    order_src = order_src.sort_values(["arm_order", "zscore"], ascending=[True, False])
    gene_order = order_src.signature_gene.tolist()
    n_pro = int((order_src.arm == "pro_inflammatory").sum())

    def bars(cond):
        d = sub[sub.condition == cond].set_index("signature_gene").reindex(gene_order).reset_index()
        return go.Bar(x=d.zscore, y=d.signature_gene, orientation="h",
                      marker=dict(color=d.zscore, colorscale="RdBu_r", cmid=0, cmin=-6, cmax=6,
                                  colorbar=dict(title="z", thickness=12, len=0.6)),
                      hovertemplate="<b>%{y}</b><br>z = %{x:.2f}<extra></extra>")

    fig = go.Figure(data=[bars("Rest")],
                    frames=[go.Frame(data=[bars(c)], name=c) for c in conds])
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK), height=780, margin=dict(l=10, r=10, t=52, b=10),
        title=dict(text=f"{gene}: signature-gene response as the T cell activates<br>"
                        f"<sup>top {n_pro} = pro-inflammatory arm · bottom = regulatory · press ▶</sup>",
                   font=dict(size=13, color=C_INK)),
        xaxis=dict(title="per-gene z-score (KD vs control)", range=[-8, 8],
                   gridcolor="#e3e8ef", zeroline=True, zerolinecolor="#9fb3c8"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=9)),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.06, xanchor="right",
                          buttons=[dict(label="▶ Play", method="animate",
                                        args=[None, {"frame": {"duration": 900, "redraw": True},
                                                     "fromcurrent": True,
                                                     "transition": {"duration": 500}}]),
                                   dict(label="❚❚", method="animate",
                                        args=[[None], {"frame": {"duration": 0, "redraw": False},
                                                       "mode": "immediate"}])])],
        sliders=[dict(active=0, currentvalue={"prefix": "condition: "}, pad={"t": 30}, x=0.1, len=0.8,
                      steps=[dict(method="animate", label=c,
                                  args=[[c], {"frame": {"duration": 500, "redraw": True},
                                              "mode": "immediate"}]) for c in conds])])
    return fig

def funnel_animated(stages):
    """#2 The method as one motion: genome → directional → druggable → genetics → novel → leads.
    stages: list of (label, count). Bars grow in one at a time top-to-bottom."""
    import plotly.graph_objects as go
    labels = [s[0] for s in stages]
    counts = [s[1] for s in stages]
    palette = ["#9fb3c8", C_BRAKE, "#2b8fb0", C_DRIVER, "#a8322f", C_ACCENT if "C_ACCENT" in globals() else "#e08a1e"]
    palette = palette[:len(stages)]
    # frames: reveal one more bar each step
    def frame_bars(k):
        xs = counts[:k] + [None] * (len(counts) - k)
        txt = [f"{c:,}" for c in counts[:k]] + [""] * (len(counts) - k)
        return go.Bar(y=labels, x=[c if c is not None else 0 for c in xs], orientation="h",
                      marker=dict(color=palette), text=txt, textposition="outside",
                      cliponaxis=False, hovertemplate="%{y}: %{x:,}<extra></extra>")
    fig = go.Figure(data=[frame_bars(1)],
                    frames=[go.Frame(data=[frame_bars(k)], name=str(k)) for k in range(1, len(stages) + 1)])
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL, font=dict(color=C_INK),
        height=420, margin=dict(l=10, r=60, t=56, b=30),
        title=dict(text="From genome to six leads — one filter at a time<br>"
                        "<sup>press ▶ to watch 11,526 genes narrow to 6 nominations</sup>",
                   font=dict(size=15, color=C_INK)),
        xaxis=dict(title="targets remaining (log)", type="log", gridcolor="#e3e8ef"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=12)),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.10, xanchor="right",
                          buttons=[dict(label="▶ Play", method="animate",
                                        args=[None, {"frame": {"duration": 800, "redraw": True},
                                                     "fromcurrent": True}])])])
    return fig

def score_buildup(row):
    """#4 Evidence components animate in as a growing stacked contribution to the integrated score."""
    import plotly.graph_objects as go
    W = {"causal_component": 0.34, "genetics_component": 0.30,
         "drug_component": 0.22, "novelty_component": 0.14}
    labs = {"causal_component": "Causal", "genetics_component": "Genetics",
            "drug_component": "Druggability", "novelty_component": "Novelty"}
    cols = {"causal_component": C_DRIVER, "genetics_component": "#2b8fb0",
            "drug_component": C_BRAKE, "novelty_component": "#8b9aa8"}
    keys = list(W)
    contrib = [float(getattr(row, k)) * W[k] for k in keys]
    dircol = C_DRIVER if row.direction == "driver_antagonize" else C_BRAKE

    def frame_traces(k):
        traces = []
        base = 0.0
        for i, key in enumerate(keys):
            val = contrib[i] if i < k else 0.0
            traces.append(go.Bar(x=[val], y=["integrated"], orientation="h", base=base,
                                 name=labs[key], marker=dict(color=cols[key]),
                                 hovertemplate=f"{labs[key]}: {val:.3f}<extra></extra>"))
            base += val
        return traces
    fig = go.Figure(data=frame_traces(0),
                    frames=[go.Frame(data=frame_traces(k), name=str(k)) for k in range(1, len(keys) + 1)])
    total = sum(contrib)
    fig.update_layout(
        barmode="stack", template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK), height=200, margin=dict(l=10, r=30, t=46, b=10),
        title=dict(text=f"{row.target_gene}: how the integrated score ({total:.3f}) is built",
                   font=dict(size=13, color=C_INK)),
        xaxis=dict(range=[0, 1.02], gridcolor="#e3e8ef"), yaxis=dict(showticklabels=False),
        legend=dict(orientation="h", y=-0.4, x=0.5, xanchor="center"),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.5, xanchor="right",
                          buttons=[dict(label="▶ Build", method="animate",
                                        args=[None, {"frame": {"duration": 650, "redraw": True},
                                                     "fromcurrent": True}])])])
    return fig

def condition_slider_map(mapd, genes=None):
    """#2 Context dynamics: directional map animated across Rest → Stim8hr → Stim48hr."""
    import plotly.graph_objects as go
    conds = ["Rest", "Stim8hr", "Stim48hr"]
    d = mapd.dropna(subset=["prog_score__Rest"]).copy()
    if genes is not None:
        d = d[d.target_gene.isin(genes)]
    if len(d) > 400:
        d = d.nlargest(400, "integrated_score")

    def frame_traces(cond):
        traces = []
        for dirval, col in [("driver_antagonize", C_DRIVER), ("brake_agonize", C_BRAKE)]:
            dd = d[d.direction == dirval]
            traces.append(go.Scatter(
                x=dd[f"emp_z__{cond}"], y=dd["genetics_component"], mode="markers",
                name=("Driver → antagonize" if dirval == "driver_antagonize" else "Brake → agonize"),
                marker=dict(size=6 + 20 * dd.integrated_score.clip(0, 1), color=col,
                            opacity=0.72, line=dict(width=0.5, color="#5a6b7b")),
                text=dd.target_gene,
                hovertemplate="<b>%{text}</b><br>emp-z %{x:.1f} · genetics %{y:.2f}<extra></extra>"))
        return traces

    fig = go.Figure(data=frame_traces("Rest"),
                    frames=[go.Frame(data=frame_traces(c), name=c) for c in conds])
    xmax = max(abs(d[[f"emp_z__{c}" for c in conds]].to_numpy().min()),
               abs(d[[f"emp_z__{c}" for c in conds]].to_numpy().max()))
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL, font=dict(color=C_INK),
        height=520, margin=dict(l=60, r=20, t=70, b=90),
        title=dict(text="Context dynamics — directional signal as CD4+ T cells activate<br>"
                        "<sup>play to move Rest → 8 h → 48 h post-stimulation</sup>",
                   font=dict(size=15, color=C_INK), x=0.01, xanchor="left", y=0.97, yanchor="top"),
        xaxis=dict(title="causal directional strength (empirical z)", gridcolor="#e3e8ef",
                   range=[-xmax * 1.05, xmax * 1.05], zeroline=True, zerolinecolor="#9fb3c8"),
        yaxis=dict(title="human-genetics support", gridcolor="#e3e8ef"),
        legend=dict(orientation="h", yanchor="top", y=-0.16, x=0.5, xanchor="center"),
        updatemenus=[dict(type="buttons", showactive=False, x=0.01, y=1.02, xanchor="left",
                          buttons=[dict(label="▶ Play", method="animate",
                                        args=[None, {"frame": {"duration": 900, "redraw": True},
                                                     "fromcurrent": True,
                                                     "transition": {"duration": 500}}]),
                                   dict(label="❚❚ Pause", method="animate",
                                        args=[[None], {"frame": {"duration": 0, "redraw": False},
                                                       "mode": "immediate"}])])],
        sliders=[dict(active=0, currentvalue={"prefix": "condition: "},
                      pad={"t": 40}, x=0.12, len=0.8,
                      steps=[dict(method="animate", label=c,
                                  args=[[c], {"frame": {"duration": 500, "redraw": True},
                                              "mode": "immediate"}]) for c in conds])])
    return fig

def disease_sunburst(dfp):
    """#4 Drill: direction → protein class → disease anchor."""
    import plotly.express as px
    d = dfp.copy()
    d["disease"] = d.ot_top_disease.fillna("(none)").str.slice(0, 28)
    d["dir_lbl"] = d.direction.map({"driver_antagonize": "Antagonize drivers",
                                    "brake_agonize": "Agonize brakes"})
    d["cls"] = d.primary_class.fillna("other")
    # cap to keep the sunburst legible
    top_disease = d.disease.value_counts().nlargest(14).index
    d = d[d.disease.isin(top_disease)]
    fig = px.sunburst(d, path=["dir_lbl", "cls", "disease"],
                      color="dir_lbl",
                      color_discrete_map={"Antagonize drivers": C_DRIVER, "Agonize brakes": C_BRAKE})
    fig.update_layout(template="plotly_white", paper_bgcolor=C_PANEL, font=dict(color=C_INK),
                      height=480, margin=dict(l=10, r=10, t=50, b=10),
                      title=dict(text="Where the nominations sit — direction → class → disease",
                                 font=dict(size=14, color=C_INK)))
    return fig

def query_trace_map(dfp, match_genes):
    """#5 Trace a query on the map: all points fade, then matches light up (2-frame animation)."""
    import plotly.graph_objects as go
    match = set(match_genes)
    def traces(revealed):
        out = []
        for dirval, col, name in [("driver_antagonize", C_DRIVER, "Driver → antagonize"),
                                  ("brake_agonize", C_BRAKE, "Brake → agonize")]:
            d = dfp[dfp.direction == dirval]
            if revealed:
                op = [0.92 if g in match else 0.06 for g in d.target_gene]
                sz = [(9 + 22 * s) if g in match else 4
                      for g, s in zip(d.target_gene, d.integrated_score.clip(0, 1))]
            else:
                op = [0.5] * len(d); sz = (5 + 10 * d.integrated_score.clip(0, 1)).tolist()
            out.append(go.Scatter(
                x=d.causal_component, y=d.genetics_component, mode="markers", name=name,
                marker=dict(size=sz, color=col, opacity=op, line=dict(width=0.5, color="white")),
                text=d.target_gene,
                hovertemplate="<b>%{text}</b><br>causal %{x:.2f} · genetics %{y:.2f}<extra></extra>"))
        return out
    fig = go.Figure(data=traces(False),
                    frames=[go.Frame(data=traces(False), name="all"),
                            go.Frame(data=traces(True), name="matches")])
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL, font=dict(color=C_INK),
        height=500, margin=dict(l=60, r=20, t=64, b=80),
        title=dict(text=f"Query trace — {len(match)} matches light up in the directional landscape",
                   font=dict(size=15, color=C_INK), x=0.01, xanchor="left"),
        xaxis=dict(title="causal directional strength", gridcolor="#e3e8ef"),
        yaxis=dict(title="human-genetics support", gridcolor="#e3e8ef"),
        legend=dict(orientation="h", yanchor="top", y=-0.14, x=0.5, xanchor="center"),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.10, xanchor="right",
                          buttons=[dict(label="▶ Trace query", method="animate",
                                        args=[["matches"], {"frame": {"duration": 700, "redraw": True},
                                                            "transition": {"duration": 500},
                                                            "mode": "immediate"}]),
                                   dict(label="⟲ Reset", method="animate",
                                        args=[["all"], {"frame": {"duration": 300, "redraw": True},
                                                        "mode": "immediate"}])])])
    return fig

def gauge_row(row):
    """#6 Genetics-constraint + tractability gauges."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    tier_val = {"strong": 0.9, "moderate": 0.55, "weak": 0.2}.get(row.genetics_tier, 0.1)
    tract = {"clinical_precedent": 1.0, "discovery_precedent": 0.7, "predicted_tractable": 0.5,
             "difficult": 0.25}.get(str(row.tractability_tier), 0.4)
    fig = make_subplots(rows=1, cols=3, specs=[[{"type": "indicator"}] * 3])
    dircol = C_DRIVER if row.direction == "driver_antagonize" else C_BRAKE
    for i, (val, lab) in enumerate([(row.integrated_score, "integrated"),
                                    (tier_val, "genetics"), (tract, "tractability")], 1):
        fig.add_trace(go.Indicator(
            mode="gauge+number", value=round(float(val), 2),
            gauge=dict(axis=dict(range=[0, 1]), bar=dict(color=dircol),
                       bgcolor="#eef2f7"),
            title=dict(text=lab, font=dict(size=12))), row=1, col=i)
    fig.update_layout(paper_bgcolor=C_PANEL, font=dict(color=C_INK),
                      height=180, margin=dict(l=10, r=10, t=10, b=10))
    return fig

# --------------------------------------------------------------------------- UI
st.set_page_config(page_title="Perturb2Target Explorer", layout="wide", page_icon="🧬")
st.markdown(f"""
<style>
  .stApp {{ background: {C_BG}; }}
  .hero {{ background: linear-gradient(105deg,#1a7f9e 0%,#12506a 100%);
           border-radius:14px; padding:20px 26px; margin-bottom:14px;
           box-shadow:0 2px 10px rgba(18,80,106,0.18); }}
  .hero h1 {{ margin:0; font-size:30px; color:#ffffff; letter-spacing:-0.5px; }}
  .hero p  {{ margin:6px 0 0 0; color:#dbeef4; font-size:14px; }}
  .kpi {{ background:{C_PANEL}; border:1px solid #dfe6ee; border-radius:12px;
          padding:14px 16px; text-align:center; box-shadow:0 1px 4px rgba(27,39,51,0.06); }}
  .kpi .v {{ font-size:26px; font-weight:700; color:#12506a; }}
  .kpi .l {{ font-size:12px; color:#5a6b7b; margin-top:2px; }}
  .pill {{ display:inline-block; padding:2px 10px; border-radius:12px; font-size:12px;
           font-weight:600; color:white; }}
</style>
""", unsafe_allow_html=True)

df = load_shortlist()
client = get_client()

st.markdown("""
<div class="hero">
  <h1>🧬 Perturb2Target — directional drug-target explorer</h1>
  <p>1,923 <b>directional</b> nominations from genome-scale CD4+ T-cell CRISPRi Perturb-seq —
  block a driver or activate a brake. Filtering runs on the real data; Claude only translates
  your words into filters and narrates a gene's own evidence, never inventing a gene or a number.</p>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.header("Filters")
    spec = {}
    dsel = st.selectbox("Direction", ["(any)", "driver_antagonize (block a driver)", "brake_agonize (activate a brake)"])
    if dsel.startswith("driver"): spec["direction"] = "driver_antagonize"
    elif dsel.startswith("brake"): spec["direction"] = "brake_agonize"
    classes = st.multiselect("Protein class", sorted(df.primary_class.dropna().unique()))
    if classes: spec["primary_class"] = classes
    tiers = st.multiselect("Genetics tier", ["strong", "moderate", "weak"])
    if tiers: spec["genetics_tier"] = tiers
    cond = st.selectbox("Peak condition", ["(any)", "Rest", "Stim8hr", "Stim48hr"])
    if cond != "(any)": spec["condition"] = cond
    if st.checkbox("Druggable only"): spec["druggable_only"] = True
    if st.checkbox("Novel / undrugged only"): spec["novel_only"] = True
    disease = st.text_input("Disease anchor contains")
    if disease.strip(): spec["disease_contains"] = disease.strip()
    maxrank = st.number_input("Max rank (0 = no limit)", min_value=0, max_value=1923, value=0)
    if maxrank > 0: spec["max_rank"] = int(maxrank)
    st.markdown("---")
    st.caption("Structures load from the bundled set (offline); missing ones fall back to AlphaFold DB. "
               "Natural-language box needs ANTHROPIC_API_KEY.")

# natural-language query box
st.markdown("#### 🔎 Ask in plain language")
if client is None:
    st.info("Set ANTHROPIC_API_KEY before launching to enable the natural-language box. "
            "The sidebar filters work fully without it.")
nlq = st.text_input(
    "e.g. 'brake-agonize kinases with strong asthma genetics and a clean patent space'",
    disabled=(client is None), label_visibility="collapsed")
if nlq and client is not None:
    try:
        nlspec = nl_to_spec(client, nlq)
        spec = nlspec  # NL query overrides sidebar when present
        chips = " ".join(f"<span class='pill' style='background:#2a4a6a'>{k}={v}</span>"
                         for k, v in nlspec.items())
        st.markdown("**Parsed filters:** " + (chips or "<i>none</i>"), unsafe_allow_html=True)
        # #5 query trace: light up the matches against the full landscape
        match_res = apply_filters(df, nlspec)
        if 0 < len(match_res) <= df.shape[0]:
            st.plotly_chart(query_trace_map(df, set(match_res.target_gene)),
                            width="stretch", config={"displayModeBar": False})
            st.caption("Press ▶ Trace query — the matches light up while the rest of the landscape fades back.")
    except Exception as e:
        st.error(f"Parse failed ({e}); using sidebar filters.")

# results
res = apply_filters(df, spec)

# KPI cards
n = len(res)
n_novel = int((res.novelty_class == "novel_undrugged").sum())
n_strong = int((res.genetics_tier == "strong").sum())
n_drug = int((res.is_druggable == True).sum())  # noqa: E712
k1, k2, k3, k4 = st.columns(4)
for col, val, lab in [(k1, n, "nominations"), (k2, n_novel, "novel / undrugged"),
                      (k3, n_strong, "strong genetics"), (k4, n_drug, "druggable")]:
    col.markdown(f"<div class='kpi'><div class='v'>{val}</div><div class='l'>{lab}</div></div>",
                 unsafe_allow_html=True)

st.markdown("")

# ------- hero: tabbed views (static map / animated context / landscape sunburst) -------
mapd = load_map_by_condition()
tab_map, tab_ctx, tab_funnel, tab_land, tab_ms = st.tabs(
    ["🗺 Directional map", "⏱ Context dynamics", "🔻 Method funnel", "🌅 Landscape",
     "🧭 MS generalization"])
with tab_map:
    left, right = st.columns([1.15, 1])
    with left:
        st.plotly_chart(directional_map(res), width="stretch", config={"displayModeBar": False})
    with right:
        st.markdown(f"**{n} matching nominations**")
        cols = ["rank", "target_gene", "direction", "integrated_score", "genetics_tier",
                "ot_top_disease", "suggested_modality", "novelty_class"]
        st.dataframe(res[cols].reset_index(drop=True), width="stretch", height=420)
with tab_ctx:
    if mapd is not None:
        st.plotly_chart(condition_slider_map(mapd, genes=set(res.target_gene)),
                        width="stretch", config={"displayModeBar": False})
        st.caption("Press ▶ Play. Targets move left/right as their causal directional signal changes "
                   "from resting to 8 h and 48 h post-stimulation — the project's context-specificity claim, live.")
    else:
        st.info("Condition-dynamics data (map_by_condition.parquet) not found next to app.py.")
with tab_funnel:
    funnel_stages = [("Genes screened", 11526), ("Directional signal (FDR<0.05)", 1923),
                     ("+ Druggable", 374), ("+ Genetics anchor", 150),
                     ("+ Novel / undrugged", 86), ("Deep-dive leads", 6)]
    st.plotly_chart(funnel_animated(funnel_stages), width="stretch", config={"displayModeBar": False})
    st.caption("Press ▶ to watch the genome-scale screen narrow, one filter at a time, to the six leads. "
               "Each stage is a strict subset of the one above it.")
with tab_land:
    st.plotly_chart(disease_sunburst(res), width="stretch", config={"displayModeBar": False})
    st.caption("Click a wedge to zoom in: direction → protein class → disease anchor.")

# per-gene deep dive: evidence bar + live 3D structure + grounded explanation
if n:
    st.markdown("---")
    st.markdown("#### 🔬 Nomination deep-dive")
    gene = st.selectbox("Select a gene", res.target_gene.tolist())
    row = res[res.target_gene == gene].iloc[0]
    dircol = C_DRIVER if row.direction == "driver_antagonize" else C_BRAKE
    action = "BLOCK (antagonize a driver)" if row.direction == "driver_antagonize" else "ACTIVATE (agonize a brake)"
    st.markdown(f"### {gene} &nbsp; <span class='pill' style='background:{dircol}'>{action}</span>",
                unsafe_allow_html=True)

    # gauges row (#6)
    st.plotly_chart(gauge_row(row), width="stretch", config={"displayModeBar": False})

    # score buildup (#4) — components animate in as a growing stacked bar
    st.plotly_chart(score_buildup(row), width="stretch", config={"displayModeBar": False})
    st.caption("Press ▶ Build — the four weighted components stack up to the integrated score.")

    cA, cB, cC = st.columns([1, 1, 1])
    with cA:
        st.plotly_chart(evidence_bar(row), width="stretch", config={"displayModeBar": False})
        st.markdown(
            f"- **rank** {int(row['rank'])} of 1,923\n"
            f"- **integrated score** {row.integrated_score:.3f}\n"
            f"- **causal** emp-z {row.peak_emp_z:.1f} (FDR {row.peak_emp_fdr:.1e}) @ {row.peak_condition}\n"
            f"- **genetics** {row.genetics_tier} — {row.ot_top_disease}\n"
            f"- **modality** {row.suggested_modality}\n"
            f"- **novelty** {row.novelty_class}")
    with cB:
        lig = ligand_genes()
        has_lig = gene in lig
        struct_choice = st.radio(
            "Structure view", (["Ligand in pocket (Boltz)"] if has_lig else []) +
            ["AlphaFold (confidence)"],
            key="structmode", horizontal=True)
        # renders AUTOMATICALLY on gene / view change — no button click needed.
        # (uncheck to suppress the WebGL viewer for a lighter page)
        show_3d = st.checkbox("Show 3D structure", value=True, key="show3d")
        if show_3d:
            if struct_choice.startswith("Ligand"):
                with st.spinner(f"Rendering Boltz complex for {gene}…"):
                    components.html(ligand_structure_html(gene), height=460)
                st.caption(f"Boltz-2 co-folded complex · the protein cartoon renders first, then the "
                           f"pre-computed ligand appears (green sticks) as the camera zooms into the pocket. "
                           f"{gene} is one of 11 pocket-bearing nominations.")
            else:
                with st.spinner(f"Loading structure for {gene}…"):
                    pdb, src = load_structure(gene)
                if pdb:
                    components.html(structure_html(pdb, "confidence"), height=440)
                    st.caption(f"{src} · red→blue = low→high model confidence (pLDDT)")
                else:
                    st.warning(f"No structure available for {gene}.")
        else:
            st.info("3D viewer hidden — tick 'Show 3D structure' to display it.")
    with cC:
        st.caption("Plain-language rationale — generated only from this gene's evidence row")
        if client is None:
            st.info("Enable ANTHROPIC_API_KEY for a grounded explanation.")
        elif st.button("💬 Explain with Claude", key="explain"):
            with st.spinner("Generating from the gene's evidence…"):
                st.markdown(explain_gene(client, row))
        with st.expander("Raw evidence row"):
            st.json({"integrated_score": round(float(row.integrated_score), 3),
                     "peak_emp_z": round(float(row.peak_emp_z), 2),
                     "peak_emp_fdr": float(row.peak_emp_fdr),
                     "genetics_tier": row.genetics_tier, "top_disease": row.ot_top_disease,
                     "suggested_modality": row.suggested_modality, "novelty_class": row.novelty_class})

    # signature-response heatmap (#1) — the mechanism reveal
    sigdf = load_signature_response()
    if sigdf is not None and gene in set(sigdf.target_gene):
        st.markdown("##### Why this call? — signature-gene response")
        hc1, hc2 = st.columns([3, 1])
        with hc2:
            animate_hm = st.checkbox("▶ Animate across conditions", key="hm_anim",
                                     help="Replay the signature response Rest → 8 h → 48 h")
            hcond = st.selectbox("Condition", ["Rest", "Stim8hr", "Stim48hr"],
                                 index=["Rest", "Stim8hr", "Stim48hr"].index(
                                     row.peak_condition if row.peak_condition in
                                     ["Rest", "Stim8hr", "Stim48hr"] else "Stim48hr"),
                                 key="hcond", disabled=animate_hm)
            st.caption("Each bar is one of the 48 signature genes' z-score under this knockdown. "
                       "A driver KD pushes the pro-inflammatory arm negative; a brake KD pushes it positive. "
                       "This is the evidence behind the directional call — not a black-box score.")
        with hc1:
            hm = (signature_heatmap_animated(sigdf, gene) if animate_hm
                  else signature_heatmap(sigdf, gene, hcond))
            if hm is not None:
                st.plotly_chart(hm, width="stretch", config={"displayModeBar": False})
    elif sigdf is not None:
        st.caption(f"(Signature-response detail is bundled for the demo gene set; {gene} is not in it.)")

    # compare mode (#5)
    with st.expander("⚖️ Compare mode — put 2–3 nominations side by side"):
        picks = st.multiselect("Choose up to 3 genes", res.target_gene.tolist(),
                               default=[gene], max_selections=3, key="cmp")
        if len(picks) >= 2:
            ccols = st.columns(len(picks))
            for cc, gp in zip(ccols, picks):
                rp = res[res.target_gene == gp].iloc[0]
                with cc:
                    dcol = C_DRIVER if rp.direction == "driver_antagonize" else C_BRAKE
                    st.markdown(f"**{gp}** <span class='pill' style='background:{dcol};font-size:10px'>"
                                f"{'block' if rp.direction=='driver_antagonize' else 'activate'}</span>",
                                unsafe_allow_html=True)
                    st.plotly_chart(evidence_bar(rp), width="stretch",
                                    config={"displayModeBar": False}, key=f"cmp_{gp}")
                    st.markdown(f"rank **{int(rp['rank'])}** · {rp.genetics_tier} genetics<br>"
                                f"<small>{rp.ot_top_disease}</small>", unsafe_allow_html=True)

with tab_ms:
    ms_sum = load_ms_summary()
    ms_nom = load_ms_nominations()
    ms_cc = load_ms_concordance()
    st.markdown(
        "#### The whole pipeline, re-run unchanged on a second disease (multiple sclerosis)")
    st.caption(
        "A generalization + validation test: the same directional score, empirical null, genetics "
        "layer and druggability filter, applied to an independent patient-derived MS signature "
        "(1.9M CD4⁺ cells across two cohorts). Two confounds this project critiques had to be "
        "controlled first — signature contamination and perturbation-footprint pleiotropy.")

    if ms_sum:
        k = st.columns(4)
        k[0].markdown(f"<div class='kpi'><div class='v'>86%</div><div class='l'>MS genetic anchors "
                      f"correctly-directed<br>(vs 58% background · p={ms_sum['enrich_binom_p']})</div></div>",
                      unsafe_allow_html=True)
        k[1].markdown(f"<div class='kpi'><div class='v'>{ms_sum['drug_concordance']}</div>"
                      f"<div class='l'>approved-MS-drug mechanism<br>concordant (p={ms_sum['drug_concordance_p']})</div></div>",
                      unsafe_allow_html=True)
        k[2].markdown(f"<div class='kpi'><div class='v'>−0.33→−0.07</div><div class='l'>pleiotropy confound "
                      f"removed<br>(reversal × footprint ρ)</div></div>", unsafe_allow_html=True)
        k[3].markdown(f"<div class='kpi'><div class='v'>{ms_sum['n_survive_fdr']}</div><div class='l'>knockdowns survive "
                      f"genome-wide FDR<br>(honest headline)</div></div>", unsafe_allow_html=True)
        st.markdown("<br>", unsafe_allow_html=True)

    fig_ms = os.path.join(HERE, "ms_generalization.png")
    if os.path.exists(fig_ms):
        st.image(fig_ms, width="stretch",
                 caption="(a) footprint-matched null removes the pleiotropy confound; "
                         "(b) MS genetic anchors enriched for therapeutic direction (86% vs 58%, p=0.007); "
                         "(c) genetic-anchored leads by footprint-matched effect size (teal = druggable); "
                         "(d) directional call matches approved-MS-drug mechanism 5/5 (p=0.031).")

    cL, cR = st.columns([1, 1])
    with cL:
        st.markdown("**MS nomination shortlist** — genetic-anchored, therapeutic-direction")
        if ms_nom is not None:
            show = ms_nom[ms_nom.ms_nomination].copy() if "ms_nomination" in ms_nom else ms_nom
            cols = [c for c in ["sym", "context", "emp_z", "genetics_tier", "primary_class", "druggable"]
                    if c in show.columns]
            st.dataframe(show[cols].round(2).reset_index(drop=True), width="stretch", height=340)
            st.caption("Druggable genetic-anchored leads: **IL2RB, TYK2, IL2RA, IL7R, TNFRSF1A** — "
                       "TYK2 is a validated approved-drug-class target.")
    with cR:
        st.markdown("**Drug-mechanism concordance** — direction vs how the approved drug acts")
        if ms_cc is not None:
            ccols = [c for c in ["gene", "drug", "drug_action", "our_call", "concordant"]
                     if c in ms_cc.columns]
            st.dataframe(ms_cc[ms_cc.concordant.notna()][ccols].reset_index(drop=True)
                         if "concordant" in ms_cc else ms_cc[ccols],
                         width="stretch", height=340)
            st.caption("Natalizumab/ITGA4 · fingolimod/S1PR1 · alemtuzumab/CD52 · TYK2 inhibitors · "
                       "anti-IL2RA — the corrected direction calls all five as *antagonize*.")

    st.info("**The thesis in miniature on new data:** undirected, uncorrected reversal recovers "
            "pleiotropic housekeeping hubs; directional scoring intersected with genetics recovers "
            "specific, druggable, correctly-directed targets. The framework applies to any "
            "CD4-T-cell-mediated disease with a definable signature.", icon="🧭")

st.markdown("---")
st.caption("All nominations are computational hypotheses; see PROSPECTIVE_VALIDATION.md for the "
           "pre-registered falsification protocol. Directional map axes are the pipeline's causal and "
           "genetics score components; bubble size is the integrated score.")
