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
import re
import json
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "target_shortlist.csv")

# brand palette
# Directional call. These two encode the therapeutic direction and are NOT part of the
# Oxford/brass brand palette — do not restyle them for visual harmony.
# ANTAGONIZE (block a driver) = Cambridge red; AGONIZE (activate a brake) = green.
# Green is Cambridge #55a51c darkened to 75% so it clears WCAG AA (4.76:1 on the app
# background); the pure brand green measures 2.88:1 and fails even the large-text floor.
C_DRIVER = "#d6083b"   # driver_antagonize — antagonize / block — Cambridge red
C_BRAKE = "#407c15"    # brake_agonize — agonize / activate — Cambridge green (AA-safe)
# Agreement colours for the agentic panel. Deliberately NOT the direction colours: with
# red now meaning "brake", reusing it for "disagrees" would read as an error flag, and
# reusing green for "agrees" would collide with the driver call.
C_AGREE = "#01305f"    # navy — reasoning layer agrees with the screen
# "our method" vs "external baseline" in the benchmark is a provenance distinction,
# not a therapeutic direction — it gets its own colour so the swap cannot desync it.
C_OURS = "#01305f"     # navy — this method, in baseline comparisons
C_DISAGREE = "#8c6608"  # dark amber — disagrees (interesting, not wrong);
                        # darker than brand amber so white pill text clears AA
C_BG = "#f5f7fa"       # app background (light)
C_PANEL = "#ffffff"    # plot / card background
C_INK = "#1b2733"      # primary text
C_OXFORD = "#002147"   # Oxford Blue — banner / brand
C_OXFORD_LT = "#01305f"  # lighter Oxford — gradients, hover states
C_OXFORD_DK = "#00152e"  # deepest Oxford — shadows, sidebar footer
C_BRASS = "#b08d57"      # muted brass — for use ON Oxford Blue (5.2:1 there)
# Darker brass for brass-on-light-background text. The lighter #b08d57 only reaches
# 2.9:1 on the page background, below the 3:1 WCAG floor even for large text, so
# anything brass sitting on white or C_BG uses this instead.
C_BRASS_DK = "#8a6a38"
C_LINE = "#e3e9f0"       # hairline borders
C_MUTED = "#5a6b7b"      # secondary text
# UI font stack. Inter is loaded from Google Fonts with a full system fallback, so the
# app still renders correctly if the webfont is unavailable (offline, or blocked).
FONT_UI = ("Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, "
           "'Helvetica Neue', Arial, sans-serif")
FONT_DISPLAY = ("'Source Serif 4', 'Iowan Old Style', 'Palatino Linotype', "
                "Palatino, Georgia, serif")
C_STRUCT_BG = "0xeef2f7"  # 3D viewer background (light)

# ----------------------------------------------------------------------------- data
# ---------- display labels ----------
# One source of truth for how a direction is worded in tables and figures. The stored
# values stay `driver_antagonize` / `brake_agonize`; only the presentation changes.
DIRECTION_LABEL = {
    "driver_antagonize": "Antagonize (Block Driver)",
    "brake_agonize": "Agonize (Activate Brake)",
}

# Column titles for user-facing tables: Camel Case words, no underscores.
COLUMN_TITLE = {
    "rank": "Rank",
    "target_gene": "Target Gene",
    "direction": "Direction",
    "integrated_score": "Integrated Score",
    "genetics_tier": "Genetics Tier",
    "ot_top_disease": "Top Disease",
    "suggested_modality": "Suggested Modality",
    "novelty_class": "Novelty Class",
    "sym": "Gene",
    "context": "Context",
    "emp_z": "Empirical Z",
    "primary_class": "Primary Class",
    "druggable": "Druggable",
    "gene": "Gene",
    "drug": "Drug",
    "drug_action": "Drug Action",
    "our_call": "Our Call",
    "concordant": "Concordant",
    "reversal_score_MS": "MS Reversal Score",
    "n_de_genes": "DE Genes",
    "kd_gene": "Knockdown Gene",
    "axis": "Method",
    "group": "Provenance",
    "auroc": "AUROC",
    "auroc_lo": "AUROC Low",
    "auroc_hi": "AUROC High",
    "ap": "AP",
    "n_positives": "Positives",
    "n_universe": "Universe",
}

def pretty_col(c):
    """Title for a column: mapped if known, else Camel Case with underscores removed."""
    if c in COLUMN_TITLE:
        return COLUMN_TITLE[c]
    return " ".join(w.capitalize() if not w.isupper() else w for w in str(c).split("_"))

def for_display(frame, cols=None):
    """Table-ready copy: direction values reworded, column titles Camel Case."""
    d = frame[cols].copy() if cols is not None else frame.copy()
    if "direction" in d.columns:
        d["direction"] = d["direction"].map(lambda v: DIRECTION_LABEL.get(v, v))
    return d.rename(columns={c: pretty_col(c) for c in d.columns})

# ---------------------------------------------------------------------------
# Cached readers.
#
# `@st.cache_data` on a loader that returns None when its file is missing will
# REMEMBER that None. If the app is running when a data file is added, the tab
# keeps reporting "not found" until the server restarts -- the file is there, the
# cache is stale. So existence is checked UNCACHED on every call, and only the
# actual read is cached, keyed on (path, mtime) so replacing a file also busts it.
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def _read_csv_cached(fp, _mtime):
    return pd.read_csv(fp)


@st.cache_data(show_spinner=False)
def _read_parquet_cached(fp, _mtime):
    return pd.read_parquet(fp)


@st.cache_data(show_spinner=False)
def _read_json_cached(fp, _mtime):
    with open(fp) as fh:
        return json.load(fh)


def _read_if_present(fp, kind="csv"):
    """Read a data file, or return None if it isn't there. Never caches the None."""
    if not os.path.exists(fp):
        return None
    mt = os.path.getmtime(fp)
    if kind == "csv":
        return _read_csv_cached(fp, mt)
    if kind == "parquet":
        return _read_parquet_cached(fp, mt)
    return _read_json_cached(fp, mt)


@st.cache_data
def load_shortlist():
    df = pd.read_csv(CSV)
    return df

def load_signature_response():
    """Per-signature-gene z-scores per nomination per condition (49 demo genes × 48 sig × 3)."""
    return _read_if_present(os.path.join(HERE, "signature_response.parquet"), "parquet")

def load_map_by_condition():
    """Per-condition directional scores for the whole shortlist (for the condition slider)."""
    return _read_if_present(os.path.join(HERE, "map_by_condition.parquet"), "parquet")

def load_agentic_calls():
    """Per-gene calls from the agentic reasoning layer (92 of the 1,923 shortlist genes).

    Each row is one gene the reasoning layer triaged: its own therapeutic_action call, a
    confidence, how many of the four evidence layers were concordant, whether it agreed
    with the screen's directional call, and -- where it disagreed -- an explicit
    primary_inconsistency. Coverage is deliberately partial: only 92 genes were triaged.
    """
    return _read_if_present(os.path.join(HERE, "agentic_triage_calls.csv"))

def load_ms_projection():
    """KD x context effects placed on the patient CD4 manifold (7,874 rows).

    Each row is one knockdown in one activation context, positioned where its
    transcriptional effect lands in the integrated patient UMAP, with the MS reversal
    score attached. Produced by the MS generalization analysis.
    """
    return _read_if_present(os.path.join(HERE, "ms_projection_app.parquet"), "parquet")

def load_ms_background():
    """Downsampled patient-cell UMAP coordinates -- the grey cloud behind the overlay.

    25,000 of 188,422 cells from the integrated CD4 atlas. Downsampled because the
    cloud only conveys the manifold's shape, and a browser scatter of 188k points is
    slow to no benefit.
    """
    return _read_if_present(os.path.join(HERE, "ms_manifold_background.parquet"), "parquet")

def load_ms_nominations():
    return _read_if_present(os.path.join(HERE, "ms_nominations_shortlist.csv"))

def load_ms_concordance():
    return _read_if_present(os.path.join(HERE, "ms_drug_concordance.csv"))

def load_ms_summary():
    return _read_if_present(os.path.join(HERE, "ms_generalization_summary.json"), "json")

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
def load_ms_anchored():
    """Genes carrying MS evidence from the patient-signature analysis (not Open Targets).
    These are a separate line of evidence: the two sources do not overlap.

    Uncached guard for the same reason as the other loaders: caching an empty set when
    the file is absent would survive the file being added.
    """
    m = _read_if_present(os.path.join(HERE, "ms_nominations_shortlist.csv"))
    if m is None or "ms_nomination" not in m.columns:
        return set()
    return set(m.loc[m.ms_nomination == True, "sym"])  # noqa: E712


MS_ANCHORED_GENES = load_ms_anchored()

# ------------------------------------------------------------------ disease filter
# Curated disease vocabulary. Values are the exact lower-cased strings that appear in
# the Open Targets annotation columns (ot_top_disease / ot_autoimmune_disease), so the
# match is an equality test on a controlled vocabulary rather than a substring guess.
# Built by enumerating the 145 distinct disease strings actually present in the data.
DISEASE_GROUPS = {
    "Multiple sclerosis": [
        "autoimmune disorder of central nervous system", "multiple sclerosis"],
    "Rheumatoid arthritis": [
        "rheumatoid arthritis"],
    "Type 1 diabetes": [
        "type 1 diabetes mellitus"],
    "Inflammatory bowel disease": [
        "colitis", "crohn disease", "inflammatory bowel disease", "inflammatory bowel disease 25",
        "neonatal inflammatory skin and bowel disease", "ulcerative colitis"],
    "Asthma / allergy": [
        "allergic rhinitis", "asthma", "asthma, nasal polyps, and aspirin intolerance",
        "atopic eczema", "childhood onset asthma"],
    "Systemic lupus erythematosus": [
        "familial chilblain lupus", "systemic lupus erythematosus"],
    "Psoriasis / psoriatic arthritis": [
        "psoriasis", "psoriasis vulgaris", "psoriatic arthritis"],
    "Autoimmune thyroid disease": [
        "autoimmune thyroid disease", "graves disease", "hashimoto thyroiditis",
        "hypoparathyroidism-deafness-renal disease syndrome", "hypothyroidism"],
    "Ankylosing spondylitis": [
        "ankylosing spondylitis"],
    "Vitiligo / alopecia areata": [
        "alopecia areata", "vitiligo"],
    "Other autoimmune / inflammatory": [
        "autoinflammatory syndrome, familial, behcet-like 1", "immunodeficiency 87 and autoimmunity",
        "mixed connective tissue disease"],
    "Primary immunodeficiency": [
        "autosomal dominant hyper-ige syndrome", "chronic mucocutaneous candidosis",
        "combined immunodeficiency due to stim1 deficiency",
        "cryptosporidiosis - chronic cholangitis - liver disease",
        "growth hormone insensitivity with immune dysregulation 1, autosomal recessive",
        "hepatic veno-occlusive disease-immunodeficiency syndrome", "herpetic encephalitis",
        "hyper-ige syndrome 6, autosomal dominant, with recurrent infections", "immunodeficiency 32b",
        "immunodeficiency 39", "immunodeficiency 69", "immunodeficiency due to cd25 deficiency",
        "neutropenia", "neutropenia, severe congenital, 8, autosomal dominant"],
    "Neurodegenerative disease": [
        "early-onset parkinsonism-intellectual disability syndrome",
        "global developmental delay, progressive ataxia, and elevated glutamine",
        "neurodegenerative disease"],
    "Neuropsychiatric / behavioural": [
        "alcohol drinking", "mathematical ability", "schizophrenia", "smoking initiation"],
    "Haematological malignancy": [
        "acute myeloid leukemia", "acute myeloid leukemia with minimal differentiation",
        "acute promyelocytic leukemia", "b-cell chronic lymphocytic leukemia", "burkitt lymphoma",
        "chronic myelogenous leukemia, bcr-abl1 positive", "diffuse large b-cell lymphoma",
        "juvenile myelomonocytic leukemia", "lymphoproliferative syndrome 1", "plasma cell myeloma"],
    "Cancer (solid tumour)": [
        "cancer", "gastric carcinoma", "head and neck squamous cell carcinoma", "melanoma",
        "ovarian neoplasm", "urinary bladder cancer"],
    "Metabolic trait / disease": [
        "congenital disorder of glycosylation with defective fucosylation",
        "congenital disorder of glycosylation, type 2v", "congenital generalized lipodystrophy type 4",
        "metabolic disease", "metabolic syndrome", "rft1-congenital disorder of glycosylation",
        "type 2 diabetes mellitus"],
    "Cardiovascular / haemostatic": [
        "atrial fibrillation", "coronary artery disorder", "hypertensive disorder", "preeclampsia",
        "pulmonary arterial hypertension", "stroke disorder", "varicose veins",
        "venous thromboembolism"],
    "Other rare / Mendelian": [
        "abnormality of the skeletal system", "aging", "azoospermia", "bardet-biedl syndrome 12",
        "bardet-biedl syndrome 7", "brooke-spiegler syndrome", "cataract",
        "charcot-marie-tooth disease type 1c", "combined oxidative phosphorylation defect type 13",
        "combined oxidative phosphorylation defect type 4", "cone rod dystrophy", "costello syndrome",
        "cystathioninuria", "cytosolic phospholipase-a2 alpha deficiency associated bleeding disorder",
        "darier disease", "decreased total leukocyte count",
        "disabling pansclerotic morphea of childhood", "distal hereditary motor neuropathy type 2",
        "emery-dreifuss muscular dystrophy", "familial cerebral saccular aneurysm",
        "fontaine progeroid syndrome", "galloway-mowat syndrome 10", "gastrointestinal disease",
        "glycogen storage disease vii", "glycosuria", "gout",
        "hereditary combined deficiency of vitamin k-dependent clotting factors",
        "hereditary spastic paraplegia 15", "hiv infectious disease", "holoprosencephaly",
        "houge-janssens syndrome 4", "intellectual disability, autosomal dominant 14",
        "isolated coq-cytochrome c reductase deficiency", "isolated sulfite oxidase deficiency",
        "joubert syndrome", "juvenile polyposis syndrome", "lymphatic system disorder",
        "lynch syndrome", "mitochondrial complex i deficiency", "multiple endocrine neoplasia type 1",
        "myopathy, tubular aggregate, 1", "neurofibromatosis type 1", "noonan syndrome",
        "osteoarthritis, hip", "osteoarthritis, knee", "osteogenesis imperfecta, type 21",
        "osteopetrosis, autosomal dominant 3", "poisoning", "polycystic liver disease 1",
        "polycystic ovary syndrome", "precordial pain", "primary biliary cholangitis",
        "primary ciliary dyskinesia",
        "pulmonary fibrosis and/or bone marrow failure syndrome, telomere-related, 7",
        "ribose-5-p isomerase deficiency", "rubinstein-taybi syndrome due to crebbp mutations",
        "schinzel-giedion syndrome", "spina bifida",
        "sting-associated vasculopathy with onset in infancy",
        "syndromic x-linked intellectual disability snyder type", "synovium disorder",
        "woodhouse-sakati syndrome"],
}

# Free-text synonyms -> canonical group, so "MS", "T1D", "IBD", "SLE" all work.
DISEASE_ALIASES = {
    "ms": "Multiple sclerosis", "multiple sclerosis": "Multiple sclerosis",
    "ra": "Rheumatoid arthritis", "rheumatoid": "Rheumatoid arthritis",
    "rheumatoid arthritis": "Rheumatoid arthritis",
    "t1d": "Type 1 diabetes", "type 1 diabetes": "Type 1 diabetes",
    "type i diabetes": "Type 1 diabetes", "iddm": "Type 1 diabetes",
    "ibd": "Inflammatory bowel disease", "crohn": "Inflammatory bowel disease",
    "crohns": "Inflammatory bowel disease", "crohn's": "Inflammatory bowel disease",
    "ulcerative colitis": "Inflammatory bowel disease", "uc": "Inflammatory bowel disease",
    "colitis": "Inflammatory bowel disease",
    "asthma": "Asthma / allergy", "allergy": "Asthma / allergy",
    "allergic": "Asthma / allergy", "atopy": "Asthma / allergy", "eczema": "Asthma / allergy",
    "sle": "Systemic lupus erythematosus", "lupus": "Systemic lupus erythematosus",
    "psoriasis": "Psoriasis / psoriatic arthritis", "psa": "Psoriasis / psoriatic arthritis",
    "thyroid": "Autoimmune thyroid disease", "hashimoto": "Autoimmune thyroid disease",
    "graves": "Autoimmune thyroid disease", "hypothyroidism": "Autoimmune thyroid disease",
    "as": "Ankylosing spondylitis", "ankylosing spondylitis": "Ankylosing spondylitis",
    "vitiligo": "Vitiligo / alopecia areata", "alopecia": "Vitiligo / alopecia areata",
    # broader disease categories
    "neurodegenerative": "Neurodegenerative disease", "neurodegeneration": "Neurodegenerative disease",
    "alzheimer": "Neurodegenerative disease", "parkinson": "Neurodegenerative disease",
    "als": "Neurodegenerative disease", "dementia": "Neurodegenerative disease",
    "leukemia": "Haematological malignancy", "leukaemia": "Haematological malignancy",
    "lymphoma": "Haematological malignancy", "myeloma": "Haematological malignancy",
    "cancer": "Cancer (solid tumour)", "tumour": "Cancer (solid tumour)",
    "tumor": "Cancer (solid tumour)", "immunodeficiency": "Primary immunodeficiency",
    "pid": "Primary immunodeficiency", "t2d": "Metabolic trait / disease",
    "type 2 diabetes": "Metabolic trait / disease", "obesity": "Metabolic trait / disease",
    "metabolic": "Metabolic trait / disease", "cardiovascular": "Cardiovascular / haemostatic",
    "hypertension": "Cardiovascular / haemostatic", "schizophrenia": "Neuropsychiatric / behavioural",
    "celiac": "Other autoimmune / inflammatory", "coeliac": "Other autoimmune / inflammatory",
    "sjogren": "Other autoimmune / inflammatory", "scleroderma": "Other autoimmune / inflammatory"
}


def resolve_disease_text(q):
    """Map a free-text disease query to a canonical group, or None if unrecognised."""
    if not q:
        return None
    k = str(q).strip().lower().rstrip("?.")
    if k in DISEASE_ALIASES:
        return DISEASE_ALIASES[k]
    for alias, grp in DISEASE_ALIASES.items():
        if len(alias) > 3 and alias in k:
            return grp
    for grp, terms in DISEASE_GROUPS.items():
        if k in grp.lower() or any(k in t for t in terms):
            return grp
    return None


def disease_mask(d, groups):
    """Rows whose Open Targets disease annotation falls in any selected group."""
    terms = {t for g in groups for t in DISEASE_GROUPS.get(g, [])}
    top = d.ot_top_disease.fillna("").str.lower()
    auto = d.ot_autoimmune_disease.fillna("").str.lower()
    return top.isin(terms) | auto.isin(terms)


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
    # Disease evidence is a UNION across sources, not an intersection. The Open Targets
    # MS annotation and the patient-signature MS set are disjoint (0 genes in common), so
    # AND-ing them would always return an empty table -- which is what a plain-language
    # query like "which MS targets should I block?" produces, since it legitimately sets
    # both flags. Both mean "carries MS evidence", so they are OR-ed.
    
    _dgroups = spec.get("disease_groups")
    if _dgroups and not isinstance(_dgroups, list):
        _dgroups = [_dgroups]
    if _dgroups or spec.get("ms_anchored_only"):
        _mask = None
        if _dgroups:
            _mask = disease_mask(d, _dgroups)
        if spec.get("ms_anchored_only"):
            _ms = d.target_gene.isin(MS_ANCHORED_GENES)
            _mask = _ms if _mask is None else (_mask | _ms)
        d = d[_mask]
    if spec.get("disease_contains"):
        # Free text: resolve to a curated group when possible (so "MS", "T1D", "IBD" work),
        # otherwise fall back to a word-boundary substring match. A plain `in` test would
        # make "RA" match "random", so the fallback is anchored on word boundaries.
        q = str(spec["disease_contains"]).strip()
        grp = resolve_disease_text(q)
        if grp:
            d = d[disease_mask(d, [grp])]
        else:
            pat = r"\b" + re.escape(q.lower())
            d = d[d.ot_top_disease.fillna("").str.lower().str.contains(pat, regex=True)
                  | d.ot_autoimmune_disease.fillna("").str.lower().str.contains(pat, regex=True)]
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
            "disease_groups": {"type": "array", "items": {"type": "string", "enum": [
                "Multiple sclerosis", "Rheumatoid arthritis",
                "Type 1 diabetes", "Inflammatory bowel disease",
                "Asthma / allergy", "Systemic lupus erythematosus",
                "Psoriasis / psoriatic arthritis", "Autoimmune thyroid disease",
                "Ankylosing spondylitis", "Vitiligo / alopecia areata",
                "Other autoimmune / inflammatory", "Primary immunodeficiency",
                "Neurodegenerative disease", "Neuropsychiatric / behavioural",
                "Haematological malignancy", "Cancer (solid tumour)",
                "Metabolic trait / disease", "Cardiovascular / haemostatic",
                "Other rare / Mendelian"]},
                "description": "Curated disease anchor groups covering every disease term in the data. Prefer this over disease_contains."},
            "ms_anchored_only": {"type": "boolean",
                "description": "True only when the user asks for multiple-sclerosis genes from the patient-signature analysis."},
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
    "For diseases prefer disease_groups with one of the enumerated values: "
    "MS/multiple sclerosis => 'Multiple sclerosis'; RA => 'Rheumatoid arthritis'; "
    "T1D => 'Type 1 diabetes'; IBD/Crohn/colitis => 'Inflammatory bowel disease'; "
    "asthma/allergy/atopy => 'Asthma / allergy'; SLE/lupus => 'Systemic lupus erythematosus'. "
    "Use disease_contains only for a disease outside that list. Never invent gene names."
)

def _secret(name, default=None):
    """Read a setting from Streamlit secrets first, then the environment.

    Both paths are needed. Locally the natural way to pass a key is an environment
    variable; on Streamlit Community Cloud the key is entered in the app's Secrets
    panel, which populates st.secrets. Streamlit does mirror root-level secrets into
    os.environ, but that has been unreliable across Cloud versions, so we read
    st.secrets explicitly rather than depending on the mirroring.
    """
    try:
        if name in st.secrets:
            return st.secrets[name]
    except Exception:
        # No secrets.toml present at all -- normal for a plain local run.
        pass
    return os.environ.get(name, default)


def get_client():
    key = _secret("ANTHROPIC_API_KEY")
    if not key:
        return None
    try:
        import anthropic
        return anthropic.Anthropic(api_key=str(key).strip())
    except Exception:
        return None


# Preference order for automatic model selection: cheap-and-fast first, since the job is
# short structured extraction, then larger models as fallbacks. Matching is by PREFIX
# against the live model list, so a dated release (e.g. claude-haiku-4-5-20251001) is
# picked up without this list needing to name the exact snapshot.
MODEL_PREFERENCE = ("claude-haiku-4-5", "claude-sonnet-4-5", "claude-sonnet-4-6",
                    "claude-haiku", "claude-sonnet", "claude-opus")


@st.cache_data(show_spinner=False, ttl=3600)
def _discover_model(_key_fingerprint):
    """Ask the API which models this key can actually use, and pick the cheapest suitable one.

    Hardcoding a model id is a maintenance trap: ids are retired, and the app then fails
    with a 404 that looks like a bug in the query parser rather than a stale constant.
    The argument is only a cache key -- the real key is read inside via _secret().
    """
    try:
        import anthropic
        key = _secret("ANTHROPIC_API_KEY")
        if not key:
            return None, []
        ids = [m.id for m in anthropic.Anthropic(api_key=str(key).strip()).models.list(limit=100).data]
        for pref in MODEL_PREFERENCE:
            for mid in ids:
                if mid.startswith(pref):
                    return mid, ids
        return (ids[0] if ids else None), ids
    except Exception:
        return None, []


def resolve_model():
    """Explicit ANTHROPIC_MODEL wins; otherwise discover a valid id from the API."""
    pinned = _secret("ANTHROPIC_MODEL")
    if pinned:
        return str(pinned).strip(), "pinned via ANTHROPIC_MODEL"
    key = _secret("ANTHROPIC_API_KEY")
    fp = (str(key)[-6:] if key else "none")
    mid, _ = _discover_model(fp)
    if mid:
        return mid, "auto-selected from the models available to your key"
    return None, "could not reach the models endpoint"


MODEL, MODEL_SOURCE = resolve_model()

def _create_with_fallback(client, **kw):
    """Call the API, and if the model id is rejected, retry once with a discovered one.

    This matters when ANTHROPIC_MODEL is pinned to an id that has since been retired:
    without the retry the app surfaces a bare 404 that reads like a parser failure.
    """
    try:
        return client.messages.create(model=MODEL, **kw)
    except Exception as e:
        if "not_found_error" not in str(e) and "model:" not in str(e):
            raise
        alt, _ = _discover_model("retry")
        if not alt or alt == MODEL:
            raise
        return client.messages.create(model=alt, **kw)


def nl_to_spec(client, query):
    msg = _create_with_fallback(
        client, max_tokens=400, system=PARSE_SYS,
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
    msg = _create_with_fallback(
        client, max_tokens=350, system=sys,
        messages=[{"role": "user", "content": "Evidence JSON:\n" + json.dumps(ev, indent=1)}],
    )
    return "".join(b.text for b in msg.content if b.type == "text")

# -------------------------------------------------------------------- plots
def directional_map(dfp, highlight=None):
    """Hero scatter: causal strength (x) vs genetics (y), size=integrated, color=direction."""
    import plotly.graph_objects as go
    fig = go.Figure()
    for dirval, col, name in [("driver_antagonize", C_DRIVER, "Antagonize (Block Driver)"),
                              ("brake_agonize", C_BRAKE, "Agonize (Activate Brake)")]:
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
        font=dict(color=C_INK, family=FONT_UI),
        height=500, margin=dict(l=60, r=20, t=70, b=90),
        title=dict(text="Directional map — causal effect vs human genetics<br>"
                        "<sup>bubble size = integrated score</sup>",
                   font=dict(size=15, color=C_OXFORD), x=0.01, xanchor="left", y=0.97, yanchor="top"),
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
                      font=dict(color=C_INK, family=FONT_UI),
                      height=230, margin=dict(l=10, r=30, t=36, b=10),
                      title=dict(text=f"{row.target_gene}: evidence components", font=dict(size=13, color=C_OXFORD)),
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
        font=dict(color=C_INK, family=FONT_UI), height=760, margin=dict(l=10, r=10, t=48, b=10),
        title=dict(text=f"{gene} KD @ {condition}: signature-gene response<br>"
                        f"<sup>top {n_pro} = pro-inflammatory arm · bottom = regulatory arm</sup>",
                   font=dict(size=13, color=C_OXFORD)),
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
        font=dict(color=C_INK, family=FONT_UI), height=780, margin=dict(l=10, r=10, t=52, b=10),
        title=dict(text=f"{gene}: signature-gene response as the T cell activates<br>"
                        f"<sup>top {n_pro} = pro-inflammatory arm · bottom = regulatory · press ▶</sup>",
                   font=dict(size=13, color=C_OXFORD)),
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
    palette = ["#9fb3c8", "#7f9cba", "#5f86ab", "#40709c", "#21598d", "#01305f"]
    palette = palette[:len(stages)]
    # frames: reveal one more bar each step
    def frame_bars(k):
        xs = counts[:k] + [None] * (len(counts) - k)
        txt = [f"{c:,}" for c in counts[:k]] + [""] * (len(counts) - k)
        return go.Bar(y=labels, x=[c if c is not None else 0 for c in xs], orientation="h",
                      marker=dict(color=palette), text=txt, textposition="outside",
                      cliponaxis=False, hovertemplate="%{y}: %{x:,}<extra></extra>")
    # Initialise with every stage shown. Starting at frame_bars(1) meant only the first bar
    # was visible until someone pressed Play, so the funnel -- the whole point of the panel --
    # was hidden behind a click. Frames are kept so the reveal can still be replayed.
    fig = go.Figure(data=[frame_bars(len(stages))],
                    frames=[go.Frame(data=[frame_bars(k)], name=str(k)) for k in range(1, len(stages) + 1)])
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL, font=dict(color=C_INK, family=FONT_UI),
        height=420, margin=dict(l=10, r=60, t=56, b=30),
        title=dict(text="From genome to six leads — one filter at a time<br>"
                        "<sup>11,526 genes narrow to 6 nominations · ↻ replays the reveal</sup>",
                   font=dict(size=15, color=C_OXFORD)),
        xaxis=dict(title="targets remaining (log)", type="log", gridcolor="#e3e8ef"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=12)),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.10, xanchor="right",
                          buttons=[dict(label="↻ Replay", method="animate",
                                        args=[None, {"frame": {"duration": 800, "redraw": True},
                                                     "fromcurrent": False,
                                                     "mode": "immediate"}])])])
    return fig

def score_buildup(row):
    """#4 Evidence components animate in as a growing stacked contribution to the integrated score."""
    import plotly.graph_objects as go
    W = {"causal_component": 0.34, "genetics_component": 0.30,
         "drug_component": 0.22, "novelty_component": 0.14}
    labs = {"causal_component": "Causal", "genetics_component": "Genetics",
            "drug_component": "Druggability", "novelty_component": "Novelty"}
    cols = {"causal_component": "#01305f", "genetics_component": "#2b8fb0",
            "drug_component": "#b08d57", "novelty_component": "#8b9aa8"}
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
    # Initialise with the FULL stack so the bar is populated on load -- the score is the
    # point of the panel, and an empty bar until someone finds the button reads as broken.
    # The frames are kept so the build-up can still be replayed as an explanation.
    fig = go.Figure(data=frame_traces(len(keys)),
                    frames=[go.Frame(data=frame_traces(k), name=str(k)) for k in range(1, len(keys) + 1)])
    total = sum(contrib)
    fig.update_layout(
        barmode="stack", template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK, family=FONT_UI), height=200, margin=dict(l=10, r=30, t=46, b=10),
        title=dict(text=f"{row.target_gene}: how the integrated score ({total:.3f}) is built",
                   font=dict(size=13, color=C_OXFORD)),
        xaxis=dict(range=[0, 1.02], gridcolor="#e3e8ef"), yaxis=dict(showticklabels=False),
        legend=dict(orientation="h", y=-0.4, x=0.5, xanchor="center"),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.5, xanchor="right",
                          buttons=[dict(label="↻ Replay build", method="animate",
                                        args=[None, {"frame": {"duration": 650, "redraw": True},
                                                     "fromcurrent": False,
                                                     "mode": "immediate"}])])])
    return fig

def p2d_backbone_fig(d):
    """Deterministic path backbone for one TF, laid out by hop distance.

    Left-to-right: the TF, then intermediate nodes at their hop position, then the
    signature-gene endpoints. Node colour marks role (target / druggable intermediate /
    other intermediate / endpoint), NOT direction these are pathway positions, not
    therapeutic calls, so the direction palette would misread here.
    """
    import plotly.graph_objects as go
    paths = d.get("top_paths", [])
    if not paths:
        return None
    inter = {n["gene"]: n for n in d.get("intermediate_nodes", [])}
    tf = d["target"]

    # hop index for every gene across all paths (earliest position wins)
    hop, endpoints = {}, set()
    for p in paths:
        seq = p["path"]
        endpoints.add(seq[-1])
        for k, g in enumerate(seq):
            hop[g] = min(hop.get(g, 99), k)
    hop[tf] = 0
    maxhop = max(hop.values())

    # vertical slot per column, ordered so edges cross as little as possible
    bycol = {}
    for g, h in hop.items():
        bycol.setdefault(h, []).append(g)
    pos = {}
    for h, genes in bycol.items():
        genes = sorted(genes)
        n = len(genes)
        for k, g in enumerate(genes):
            pos[g] = (h, (n - 1) / 2 - k)

    fig = go.Figure()
    seen = set()
    for p in paths:
        seq = p["path"]
        for a, b in zip(seq[:-1], seq[1:]):
            if (a, b) in seen:
                continue
            seen.add((a, b))
            x0, y0 = pos[a]; x1, y1 = pos[b]
            fig.add_trace(go.Scatter(
                x=[x0, x1], y=[y0, y1], mode="lines", hoverinfo="skip",
                line=dict(color="#c9d4e0", width=1.6), showlegend=False))

    groups = {"target": ([], [], [], C_OXFORD, 30),
              "druggable": ([], [], [], "#b08d57", 24),
              "intermediate": ([], [], [], "#7f9cba", 19),
              "signature gene": ([], [], [], "#2b8fb0", 16)}
    for g, (x, y) in pos.items():
        if g == tf:
            key, txt = "target", f"<b>{g}</b><br>target TF ({d.get('direction','')})"
        elif g in inter:
            n = inter[g]
            dr = n.get("druggable") or {}
            key = "druggable" if dr else "intermediate"
            bits = [f"<b>{g}</b>", "intermediate node"]
            if dr:
                bits.append(f"ChEMBL max phase {dr.get('max_phase','?')} &middot; {dr.get('action','?')}")
            if n.get("observed"):
                bits.append(f"own program effect {n.get('peak_signed_score', float('nan')):+.2f}")
                bits.append("direction-consistent" if n.get("direction_consistent")
                            else "direction-INconsistent")
            ez = n.get("edge_zscore_under_target_KD")
            if ez is not None:
                bits.append(f"edge z under {tf} KD = {ez:+.2f}"
                            + ("  (confirmed)" if n.get("edge_functionally_confirmed") else ""))
            txt = "<br>".join(bits)
        elif g in endpoints:
            key, txt = "signature gene", f"<b>{g}</b><br>program signature gene"
        else:
            key, txt = "intermediate", f"<b>{g}</b>"
        groups[key][0].append(x); groups[key][1].append(y); groups[key][2].append(txt)

    for key, (xs, ys, txts, col, size) in groups.items():
        if not xs:
            continue
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="markers", name=key,
            marker=dict(size=size, color=col, line=dict(color="white", width=1.6)),
            text=txts, hovertemplate="%{text}<extra></extra>"))
    for g, (x, y) in pos.items():
        fig.add_annotation(x=x, y=y, text=f"<b>{g}</b>", showarrow=False,
                           yshift=-20, font=dict(size=10.5, color=C_INK))

    fig.update_layout(
        height=380, margin=dict(l=10, r=10, t=54, b=34),
        paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        title=dict(text=(f"Mechanistic backbone: {tf} to the program signature"
                         f"<br><sup>{d.get('n_reachable_signature_genes','?')} signature genes "
                         f"reachable &middot; gold = druggable intermediate &middot; "
                         f"hover any node</sup>"),
                   font=dict(size=14, color=C_OXFORD), x=0.01),
        legend=dict(orientation="h", y=1.10, x=0.34, font=dict(size=10.5)),
        xaxis=dict(visible=False, range=[-0.45, maxhop + 0.45]),
        yaxis=dict(visible=False))
    return fig


def p2d_baseline_fig(bl):
    """Path2Drug's node choice vs the naive highest-degree heuristic, per TF."""
    import plotly.graph_objects as go
    d = bl.dropna(subset=["cell2network_node"]).copy()
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[0] * len(d), y=d.tf, mode="markers+text",
        marker=dict(size=15, color="#9fb3c8", line=dict(color="white", width=1.4)),
        text=d.baseline_hub_node.fillna("(none)"), textposition="middle left",
        textfont=dict(size=10.5, color=C_MUTED),
        name="naive highest-degree neighbour",
        hovertemplate="%{y}: baseline picks %{text}<extra></extra>"))
    fig.add_trace(go.Scatter(
        x=[1] * len(d), y=d.tf, mode="markers+text",
        marker=dict(size=15, color=C_BRASS, line=dict(color="white", width=1.4)),
        text=d.cell2network_node, textposition="middle right",
        textfont=dict(size=10.5, color=C_OXFORD),
        name="Path2Drug", hovertemplate="%{y}: Path2Drug picks %{text}<extra></extra>"))
    for r in d.itertuples():
        fig.add_shape(type="line", x0=0, x1=1, y0=r.tf, y1=r.tf,
                      line=dict(color="#e3e9f0", width=1.2), layer="below")
    nb = d.baseline_hub_node.nunique(dropna=True)
    np_ = d.cell2network_node.nunique(dropna=True)
    fig.update_layout(
        height=340, margin=dict(l=10, r=10, t=62, b=30),
        paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        title=dict(text=("Pathway-aware vs naive node choice"
                         f"<br><sup>the degree heuristic collapses to {nb} distinct node(s); "
                         f"Path2Drug returns {np_}</sup>"),
                   font=dict(size=14, color=C_OXFORD), x=0.01),
        legend=dict(orientation="h", y=1.13, x=0.30, font=dict(size=10.5)),
        xaxis=dict(visible=False, range=[-0.62, 1.62]),
        yaxis=dict(autorange="reversed", tickfont=dict(size=11, color=C_INK),
                   showgrid=False, zeroline=False))
    return fig


def p2d_edge_fig(ec):
    """Functional edge test: does knocking the TF down move the intermediate's transcript?"""
    import plotly.graph_objects as go
    d = ec.dropna(subset=["edge_zscore"]).copy()
    d = d.reindex(d.edge_zscore.abs().sort_values().index)
    lab = d.tf + " \u2192 " + d.intermediate
    # `edge_confirmed` is an object column holding Python bools plus NaN for untestable
    # edges. `.astype(bool)` would make every NaN truthy and quadruple the count, so
    # compare to True explicitly.
    conf = d.edge_confirmed == True  # noqa: E712
    fig = go.Figure(go.Bar(
        x=d.edge_zscore, y=lab, orientation="h",
        marker=dict(color=["#2e7d32" if c else "#9fb3c8" for c in conf]),
        hovertemplate="%{y}<br>edge z = %{x:.2f}<extra></extra>"))
    for v in (-1.96, 1.96):
        fig.add_vline(x=v, line=dict(color="#8b9aa8", width=1, dash="dot"))
    fig.update_layout(
        height=max(320, 17 * len(d)), margin=dict(l=10, r=20, t=62, b=40),
        paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        title=dict(text=(f"Functional edge confirmation &mdash; {int(conf.sum())} of {len(d)} testable edges pass"
                         "<br><sup>z-score of the intermediate's transcript under its target's "
                         "knockdown; dotted lines |z| = 1.96</sup>"),
                   font=dict(size=14, color=C_OXFORD), x=0.01),
        xaxis=dict(title="edge z-score", gridcolor=C_LINE, zerolinecolor="#8b9aa8"),
        yaxis=dict(tickfont=dict(size=9.5, color=C_INK)), showlegend=False)
    return fig


def multibaseline_fig(mb):
    """Our scores against three external baselines, with bootstrap CIs.

    Bars are ordered by AUROC so the reader sees immediately that two baselines beat the
    full integrated score. Error bars are the bootstrap CIs already in the results table --
    they overlap heavily, which is the actual finding.
    """
    import plotly.graph_objects as go
    d = mb.sort_values("auroc")
    cols = [C_OURS if g.startswith("ours") else "#9fb3c8" for g in d.group]
    fig = go.Figure(go.Bar(
        x=d.auroc, y=d.axis, orientation="h", marker=dict(color=cols),
        error_x=dict(type="data", symmetric=False,
                     array=(d.auroc_hi - d.auroc).tolist(),
                     arrayminus=(d.auroc - d.auroc_lo).tolist(),
                     color="#5a6b7b", thickness=1.2, width=4),
        customdata=d[["auroc_lo", "auroc_hi", "ap", "n_positives"]].values,
        hovertemplate=("<b>%{y}</b><br>AUROC %{x:.3f} "
                       "(95%% CI %{customdata[0]:.3f}&ndash;%{customdata[1]:.3f})<br>"
                       "AP %{customdata[2]:.3f} &middot; %{customdata[3]} positives"
                       "<extra></extra>")))
    # Value labels sit past the END OF THE WHISKER, not at the bar end. textposition="outside"
    # places them at auroc, which is 0.02-0.16 short of auroc_hi -- so the number printed on
    # top of its own error bar. Anchoring to auroc_hi guarantees clearance for every row.
    for _, r in d.iterrows():
        fig.add_annotation(x=r.auroc_hi + 0.008, y=r.axis, text=f"{r.auroc:.3f}",
                           showarrow=False, xanchor="left", yanchor="middle",
                           font=dict(size=11, color=C_INK))
    fig.add_vline(x=0.5, line=dict(color="#8b9aa8", width=1, dash="dot"))
    # Pin the "chance" label to the top of the plot area rather than below the lowest
    # category, where it competed with the axis title.
    fig.add_annotation(x=0.5, yref="paper", y=1.02, text="chance", showarrow=False,
                       xanchor="center", yanchor="bottom",
                       font=dict(size=9, color=C_MUTED))
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK, family=FONT_UI), height=340,
        margin=dict(l=10, r=80, t=64, b=46),
        title=dict(text=("Recovering known drug targets: ours vs external baselines"
                         "<br><sup>navy = this method &middot; network centrality and Open Targets "
                         "genetics both score higher</sup>"),
                   font=dict(size=14, color=C_OXFORD), x=0.01, xanchor="left"),
        xaxis=dict(title="AUROC (19 known targets among 1,923)", range=[0.35, 1.09],
                   gridcolor="#e3e8ef"),
        yaxis=dict(tickfont=dict(size=11)))
    return fig

def load_p2d_summary():
    """Per-TF backbone summary: intermediates found, how many druggable / consistent."""
    return _read_if_present(os.path.join(HERE, "path2drug_enriched_summary.csv"))


def load_p2d_baseline():
    """Path2Drug's chosen node vs the naive highest-degree-druggable-neighbour heuristic."""
    return _read_if_present(os.path.join(HERE, "path2drug_vs_baseline.csv"))


def load_p2d_edges():
    """Functional test of each target->intermediate edge against the atlas itself."""
    return _read_if_present(os.path.join(HERE, "path2drug_edge_confirmation.csv"))


def load_p2d_target(tf):
    """Full deterministic backbone for one TF: paths, intermediates, druggability."""
    return _read_if_present(os.path.join(HERE, "path2drug", f"p2d_{tf}.json"), "json")


def p2d_available():
    """TFs with a released backbone, in the summary's own order (best-supported first).

    Deliberately uncached: this is a directory listing, and caching an empty result
    would keep the tab reporting "not found" after the data is copied in.
    """
    d = os.path.join(HERE, "path2drug")
    if not os.path.isdir(d):
        return []
    got = sorted(f[4:-5] for f in os.listdir(d) if f.startswith("p2d_") and f.endswith(".json"))
    s = load_p2d_summary()
    if s is not None:
        order = [t for t in s.tf.tolist() if t in got]
        return order + [t for t in got if t not in order]
    return got


def load_multibaseline():
    """AUROC/AP for our scores vs three external baselines (6 rows)."""
    return _read_if_present(os.path.join(HERE, "multibaseline_comparison_results.csv"))

def load_weight_sensitivity():
    """Summary of the +/-40% weight-perturbation robustness analysis."""
    return _read_if_present(os.path.join(HERE, "weight_sensitivity_summary.json"), "json")

def weight_reranking_fig(dfp, w_causal, w_gen, w_drug, w_nov, top_n=15):
    """Live re-ranking under user-chosen integration weights.

    Recomputes the integrated score from the four stored components and shows how the
    top-N shifts against the published weights (0.34/0.30/0.22/0.14). Weights are
    normalised to sum to 1 so the score stays on its original scale and the comparison
    is like-for-like.
    """
    import plotly.graph_objects as go
    tot = w_causal + w_gen + w_drug + w_nov
    if tot <= 0:
        tot = 1.0
    wc, wg, wd, wn = (w_causal / tot, w_gen / tot, w_drug / tot, w_nov / tot)
    d = dfp.copy()
    d["custom_score"] = (wc * d.causal_component + wg * d.genetics_component
                         + wd * d.drug_component + wn * d.novelty_component)
    d["custom_rank"] = d.custom_score.rank(ascending=False, method="min").astype(int)
    pub = d.nsmallest(top_n, "rank")
    top = d.nsmallest(top_n, "custom_rank").sort_values("custom_rank")
    # Colour by whether a gene is in the published top-N: entrants are the interesting cases.
    pubset = set(pub.target_gene)
    cols = [C_OXFORD if g in pubset else C_BRASS_DK for g in top.target_gene]
    fig = go.Figure(go.Bar(
        x=top.custom_score, y=top.target_gene, orientation="h",
        marker=dict(color=cols),
        text=[f"{s:.3f}  (was #{int(r)})" for s, r in zip(top.custom_score, top["rank"])],
        textposition="outside", cliponaxis=False,
        customdata=top[["rank", "custom_rank", "direction"]].values,
        hovertemplate=("<b>%{y}</b><br>custom score %{x:.3f}<br>"
                       "published rank %{customdata[0]} &rarr; custom rank %{customdata[1]}"
                       "<extra></extra>")))
    n_new = int(sum(1 for g in top.target_gene if g not in pubset))
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK, family=FONT_UI), height=440,
        margin=dict(l=10, r=130, t=62, b=40),
        title=dict(text=(f"Top {top_n} under your weights"
                         f"<br><sup>{n_new} gene(s) not in the published top {top_n} "
                         f"&middot; brass = new entrant, navy = also published</sup>"),
                   font=dict(size=14, color=C_OXFORD), x=0.01, xanchor="left"),
        xaxis=dict(title="re-weighted integrated score", gridcolor="#e3e8ef"),
        yaxis=dict(autorange="reversed", tickfont=dict(size=11)))
    return fig, d

def agentic_layers_fig(call):
    """Four-layer concordance for one gene, as a compact horizontal indicator.

    layers_concordant is a count (1-4), not per-layer detail, so this shows filled vs
    empty slots rather than naming which layer dissented -- claiming more resolution than
    the data has would be a fabrication.
    """
    import plotly.graph_objects as go
    n = int(call["layers_concordant"])
    agree = str(call["agrees_with_screen"])
    col = {"yes": C_AGREE, "no": C_DISAGREE, "uncertain": "#8b9aa8"}.get(agree, "#8b9aa8")
    fig = go.Figure()
    for i in range(4):
        filled = i < n
        fig.add_trace(go.Bar(
            x=[1], y=["layers"], orientation="h", showlegend=False,
            marker=dict(color=col if filled else "#eaeef4",
                        line=dict(width=1, color="#ffffff")),
            hovertemplate=(f"layer {i+1}: concordant<extra></extra>" if filled
                           else f"layer {i+1}: not concordant<extra></extra>")))
    fig.update_layout(
        barmode="stack", template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        height=72, margin=dict(l=6, r=6, t=30, b=6),
        title=dict(text=f"{n} of 4 evidence layers concordant",
                   font=dict(size=11.5, color=C_OXFORD), x=0.0, xanchor="left"),
        xaxis=dict(visible=False), yaxis=dict(visible=False))
    return fig

def agentic_overview_fig(ag):
    """Distribution of the reasoning layer's calls across all 92 triaged genes."""
    import plotly.graph_objects as go
    order = ["high", "medium", "low"]
    agree_order = ["yes", "no", "uncertain"]
    cols = {"yes": C_AGREE, "no": C_DISAGREE, "uncertain": "#8b9aa8"}
    labs = {"yes": "agrees with screen", "no": "disagrees", "uncertain": "uncertain"}
    fig = go.Figure()
    for a in agree_order:
        counts = [int(((ag.confidence == c) & (ag.agrees_with_screen == a)).sum()) for c in order]
        fig.add_trace(go.Bar(x=order, y=counts, name=labs[a], marker=dict(color=cols[a]),
                             text=[c if c else "" for c in counts], textposition="inside",
                             hovertemplate="%{x} confidence · " + labs[a] + ": %{y}<extra></extra>"))
    fig.update_layout(
        barmode="stack", template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK, family=FONT_UI), height=300,
        margin=dict(l=50, r=20, t=58, b=40),
        title=dict(text=("Reasoning layer: 92 genes triaged"
                         "<br><sup>disagreements are the useful cases &mdash; every one carries an "
                         "explicit inconsistency flag</sup>"),
                   font=dict(size=14, color=C_OXFORD), x=0.01, xanchor="left"),
        xaxis=dict(title="confidence"), yaxis=dict(title="genes"),
        legend=dict(orientation="h", y=-0.22, x=0.5, xanchor="center", font=dict(size=10)))
    return fig

def ms_manifold_overlay(proj, bg, context="All", top_n=20):
    """Interactive version of the perturbation-overlay figure.

    Panel a of the static figure showed every KDxcontext coloured by MS reversal score
    over the patient cloud; panel b labelled the top reversing knockdowns. Here both are
    one figure with the top-N called out, and the context is selectable -- the static
    figure had to pool contexts, which hid the fact that a knockdown's position moves
    as the T cell activates.
    """
    import plotly.graph_objects as go
    d = proj if context == "All" else proj[proj.context == context]
    top = d.nlargest(top_n, "reversal_score_MS")
    fig = go.Figure()
    if bg is not None:
        fig.add_trace(go.Scattergl(
            x=bg.umap_x, y=bg.umap_y, mode="markers", name="patient CD4 cells",
            marker=dict(size=2.4, color="#d3dae2", opacity=0.55), hoverinfo="skip"))
    # Diverging scale centred on zero: sign is the whole point (reversing vs worsening),
    # so the midpoint must sit at 0 rather than at the data median.
    lim = float(max(abs(d.reversal_score_MS.min()), abs(d.reversal_score_MS.max())))
    fig.add_trace(go.Scattergl(
        x=d.umap_x, y=d.umap_y, mode="markers", name="knockdown x context",
        marker=dict(size=5.5, color=d.reversal_score_MS, colorscale="RdBu_r",
                    cmin=-lim, cmax=lim, showscale=True,
                    colorbar=dict(title=dict(text="MS reversal<br>score", side="right"),
                                  thickness=12, len=0.62, x=1.015,
                                  tickfont=dict(size=9)),
                    line=dict(width=0.3, color="rgba(255,255,255,0.5)")),
        text=d.kd_gene, customdata=d[["context", "n_de_genes"]].values,
        hovertemplate=("<b>%{text}</b><br>%{customdata[0]}<br>"
                       "reversal %{marker.color:.4f}<br>"
                       "%{customdata[1]} DE genes<extra></extra>")))
    fig.add_trace(go.Scattergl(
        x=top.umap_x, y=top.umap_y, mode="markers+text", name=f"top {top_n} reversing",
        marker=dict(size=11, color="rgba(0,0,0,0)",
                    line=dict(width=1.8, color=C_OXFORD)),
        text=top.kd_gene, textposition="top center",
        textfont=dict(size=9, color=C_OXFORD, family=FONT_UI),
        hovertemplate="<b>%{text}</b><br>top reversing<extra></extra>"))
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL,
        font=dict(color=C_INK, family=FONT_UI), height=560,
        margin=dict(l=52, r=90, t=62, b=48),
        title=dict(text=(f"Candidate MS-reversing knockdowns on the patient CD4 manifold"
                         f"<br><sup>{len(d):,} knockdown x context effects &middot; "
                         f"{context.lower() if context != 'All' else 'all contexts'} &middot; "
                         f"grey = 25,000 patient cells</sup>"),
                   font=dict(size=15, color=C_OXFORD), x=0.01, xanchor="left"),
        xaxis=dict(title="UMAP1", showgrid=False, zeroline=False),
        yaxis=dict(title="UMAP2", showgrid=False, zeroline=False),
        legend=dict(orientation="h", y=-0.13, x=0.5, xanchor="center",
                    font=dict(size=10)))
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
                name=DIRECTION_LABEL[dirval],
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
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL, font=dict(color=C_INK, family=FONT_UI),
        height=520, margin=dict(l=60, r=20, t=70, b=90),
        title=dict(text="Context dynamics — directional signal as CD4+ T cells activate<br>"
                        "<sup>play to move Rest → 8 h → 48 h post-stimulation</sup>",
                   font=dict(size=15, color=C_OXFORD), x=0.01, xanchor="left", y=0.97, yanchor="top"),
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
    d["dir_lbl"] = d.direction.map(DIRECTION_LABEL)
    d["cls"] = d.primary_class.fillna("other")
    # cap to keep the sunburst legible
    top_disease = d.disease.value_counts().nlargest(14).index
    d = d[d.disease.isin(top_disease)]
    fig = px.sunburst(d, path=["dir_lbl", "cls", "disease"],
                      color="dir_lbl",
                      color_discrete_map={DIRECTION_LABEL["driver_antagonize"]: C_DRIVER,
                                            DIRECTION_LABEL["brake_agonize"]: C_BRAKE})
    fig.update_layout(template="plotly_white", paper_bgcolor=C_PANEL, font=dict(color=C_INK, family=FONT_UI),
                      height=480, margin=dict(l=10, r=10, t=50, b=10),
                      title=dict(text="Where the nominations sit — direction → class → disease",
                                 font=dict(size=14, color=C_OXFORD)))
    return fig

def query_trace_map(dfp, match_genes):
    """#5 Trace a query on the map: matches highlighted, non-matches faded back.

    Renders in the revealed state so the answer is visible immediately; the two frames
    remain so the reveal can be replayed or the full landscape shown undimmed.
    """
    import plotly.graph_objects as go
    match = set(match_genes)
    def traces(revealed):
        out = []
        for dirval, col, name in [("driver_antagonize", C_DRIVER, DIRECTION_LABEL["driver_antagonize"]),
                                  ("brake_agonize", C_BRAKE, DIRECTION_LABEL["brake_agonize"])]:
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
    # Initialise in the REVEALED state: matches highlighted, non-matches faded back. The
    # panel's whole content is which targets your query picked out, and traces(False) is a
    # uniform cloud that shows none of it -- so the default has to be the revealed frame,
    # not the pre-reveal one. Both frames stay available as replay / show-all.
    fig = go.Figure(data=traces(True),
                    frames=[go.Frame(data=traces(False), name="all"),
                            go.Frame(data=traces(True), name="matches")])
    fig.update_layout(
        template="plotly_white", paper_bgcolor=C_PANEL, plot_bgcolor=C_PANEL, font=dict(color=C_INK, family=FONT_UI),
        height=500, margin=dict(l=60, r=20, t=64, b=80),
        title=dict(text=f"Query trace — {len(match)} matches light up in the directional landscape",
                   font=dict(size=15, color=C_OXFORD), x=0.01, xanchor="left"),
        xaxis=dict(title="causal directional strength", gridcolor="#e3e8ef"),
        yaxis=dict(title="human-genetics support", gridcolor="#e3e8ef"),
        legend=dict(orientation="h", yanchor="top", y=-0.14, x=0.5, xanchor="center"),
        updatemenus=[dict(type="buttons", showactive=False, x=0.98, y=1.10, xanchor="right",
                          buttons=[dict(label="↻ Replay reveal", method="animate",
                                        args=[["all", "matches"],
                                              {"frame": {"duration": 700, "redraw": True},
                                               "transition": {"duration": 500},
                                               "mode": "immediate"}]),
                                   dict(label="◌ Show all", method="animate",
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
    fig.update_layout(paper_bgcolor=C_PANEL, font=dict(color=C_INK, family=FONT_UI),
                      height=180, margin=dict(l=10, r=10, t=10, b=10))
    return fig

# --------------------------------------------------------------------------- UI
st.set_page_config(page_title="Perturb2Target Explorer", layout="wide", page_icon="🧬")
# NOTE ON THE STYLE BLOCK BELOW: it must contain NO blank lines. Streamlit renders
# markdown, and in markdown a blank line terminates a raw-HTML block -- everything after
# the first blank line inside <style> gets emitted as visible text instead of being applied.
# Use the /* ---------- */ comments as section separators, not an empty line.  
st.markdown(f"""
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&display=swap" rel="stylesheet">
<style>
  /* ---------- foundation ---------- */
  html, body, .stApp, [class*="css"] {{ font-family:{FONT_UI}; }}
  .stApp {{ background:{C_BG}; color:{C_INK}; }}
  .block-container {{ padding-top:1.6rem; padding-bottom:3rem; max-width:1500px; }}
  h1,h2,h3,h4,h5 {{ font-family:{FONT_UI}; color:{C_OXFORD}; letter-spacing:-0.15px; }}
  a {{ color:{C_OXFORD_LT}; text-decoration:none; border-bottom:1px solid rgba(1,48,95,0.25); }}
  a:hover {{ border-bottom-color:{C_BRASS}; }}
  /* ---------- masthead ---------- */
  .hero {{ background:linear-gradient(118deg,{C_OXFORD} 0%,{C_OXFORD_LT} 62%,#013b73 100%);
           border-radius:4px; padding:30px 34px 28px 34px; margin-bottom:20px;
           border-top:3px solid {C_BRASS};
           box-shadow:0 10px 30px -12px rgba(0,21,46,0.55); position:relative; overflow:hidden; }}
  /* Faint engraved rule under the wordmark — the only ornament, kept subtle. */
  .hero:after {{ content:""; position:absolute; left:34px; right:34px; top:0; height:1px;
                 background:linear-gradient(90deg,rgba(176,141,87,0.55),transparent 60%); }}
  .hero .eyebrow {{ font-size:10.5px; font-weight:600; letter-spacing:2.4px;
                    text-transform:uppercase; color:{C_BRASS}; margin:0 0 9px 0; }}
  .hero h1 {{ font-family:{FONT_DISPLAY}; font-weight:700; margin:0; font-size:35px;
              color:#ffffff; letter-spacing:-0.4px; line-height:1.12; }}
  /* Author byline: set in the display serif at a size that reads as attribution under the
     wordmark rather than as a second heading. 13.7:1 on Oxford Blue. */
  .hero .byline {{ font-family:{FONT_DISPLAY}; font-size:14.5px; font-weight:600;
                   color:#e7eef6; margin:9px 0 0 0; letter-spacing:0.2px; }}
  .hero .rule {{ width:52px; height:2px; background:{C_BRASS}; margin:14px 0 13px 0;
                 border-radius:2px; }}
  .hero p {{ margin:0; color:#c6d6e6; font-size:13.5px; line-height:1.62; max-width:66em;
             font-weight:400; }}
  .hero p b {{ color:#ffffff; font-weight:600; }}
  /* ---------- metric cards ---------- */
  .kpi {{ background:{C_PANEL}; border:1px solid {C_LINE};
          border-top:2px solid {C_OXFORD}; border-radius:3px;
          padding:15px 16px 13px 16px; text-align:center;
          box-shadow:0 1px 3px rgba(0,33,71,0.05);
          transition:box-shadow 160ms ease, transform 160ms ease; }}
  .kpi:hover {{ box-shadow:0 6px 18px -8px rgba(0,33,71,0.28); transform:translateY(-1px); }}
  .kpi .v {{ font-family:{FONT_DISPLAY}; font-size:31px; font-weight:700;
             color:{C_OXFORD}; line-height:1.05; font-variant-numeric:tabular-nums; }}
  .kpi .l {{ font-size:10.5px; color:{C_MUTED}; margin-top:5px; font-weight:600;
             letter-spacing:0.9px; text-transform:uppercase; }}
  /* ---------- section headings ---------- */
  .sec {{ display:flex; align-items:baseline; gap:11px; margin:30px 0 13px 0;
          padding-bottom:8px; border-bottom:1px solid {C_LINE}; }}
  /* Section numeral uses the darker brass: on the light page background the lighter
     brass measures 2.9:1, under the 3:1 floor. */
  .sec .n {{ font-family:{FONT_DISPLAY}; font-size:12.5px; font-weight:700;
             color:{C_BRASS_DK}; letter-spacing:0.5px; }}
  .sec h3 {{ margin:0; font-size:17.5px; font-weight:600; color:{C_OXFORD}; }}
  .sec .sub {{ font-size:12px; color:{C_MUTED}; margin-left:auto; text-align:right; }}
  /* ---------- direction pills ---------- */
  .pill {{ display:inline-block; padding:3px 11px; border-radius:2px; font-size:11px;
           font-weight:600; color:white; letter-spacing:0.5px; text-transform:uppercase; }}
  /* ---------- sidebar ---------- */
  section[data-testid="stSidebar"] {{ background:{C_OXFORD}; border-right:1px solid {C_OXFORD_DK}; }}
  section[data-testid="stSidebar"] * {{ color:#dbe6f1; }}
  section[data-testid="stSidebar"] h1,
  section[data-testid="stSidebar"] h2,
  section[data-testid="stSidebar"] h3 {{ color:#ffffff; font-size:16px; font-weight:600;
        letter-spacing:0.3px; }}
  /* Group label: small caps in brass, so the panel reads as sections not a widget dump.
     Sized at 13px -- 10.5px uppercase was legible up close but too quiet as a heading,
     and it has to out-rank the 12.5px widget labels below it to read as a group title. */
  section[data-testid="stSidebar"] h3 {{ font-size:13px; letter-spacing:1.4px;
        font-weight:700; text-transform:uppercase; color:{C_BRASS};
        margin-top:24px; margin-bottom:4px; }}
  section[data-testid="stSidebar"] label {{ font-size:12.5px; color:#b8cade; font-weight:500; }}
  section[data-testid="stSidebar"] hr {{ border-color:rgba(255,255,255,0.13); margin:18px 0; }}
  /* Filter fields are LIGHT with Oxford-blue text. They were a translucent dark fill with
     white text, which made typed text and selected values invisible: the multiselect search
     input and the BaseWeb dropdown popover both render on a white surface, so white-on-white.
     Recolouring the text alone is not enough -- Oxford text on the old dark fill measures
     1.14:1. The fill has to be light for the text to read, so field and popover now match. */
  section[data-testid="stSidebar"] [data-baseweb="select"] > div,
  section[data-testid="stSidebar"] [data-baseweb="input"],
  section[data-testid="stSidebar"] input,
  section[data-testid="stSidebar"] textarea {{
        background:#ffffff !important;
        border-color:rgba(255,255,255,0.55) !important; color:{C_OXFORD} !important;
        border-radius:2px !important; }}
  /* The `* {{ color }}` sidebar rule above would otherwise repaint the value/search text
     inside the control, so name those elements explicitly. */
  section[data-testid="stSidebar"] [data-baseweb="select"] input,
  section[data-testid="stSidebar"] [data-baseweb="select"] div[class*="Value"],
  section[data-testid="stSidebar"] [data-baseweb="select"] div[class*="value"],
  section[data-testid="stSidebar"] input,
  section[data-testid="stSidebar"] [data-testid="stNumberInputField"] {{
        color:{C_OXFORD} !important; -webkit-text-fill-color:{C_OXFORD} !important; }}
  section[data-testid="stSidebar"] [data-baseweb="select"] svg {{ fill:{C_OXFORD_LT} !important; }}
  section[data-testid="stSidebar"] ::placeholder {{
        color:#5a6b7b !important; -webkit-text-fill-color:#5a6b7b !important; opacity:1 !important; }}
  section[data-testid="stSidebar"] [data-testid="stNumberInput"] button {{
        background:#eef3f9 !important; color:{C_OXFORD} !important; }}
  section[data-testid="stSidebar"] [data-testid="stNumberInput"] button svg {{ fill:{C_OXFORD} !important; }}
  /* Dropdown menus render in a portal at document root, not inside the sidebar, so these
     are global rules -- a sidebar-scoped selector would never reach them. */
  div[data-baseweb="popover"] li, div[data-baseweb="popover"] [role="option"],
  ul[role="listbox"] li, div[role="listbox"] [role="option"] {{
        color:{C_OXFORD} !important; -webkit-text-fill-color:{C_OXFORD} !important;
        background:#ffffff !important; font-size:12.5px !important; }}
  div[data-baseweb="popover"] li:hover, div[data-baseweb="popover"] [role="option"]:hover,
  ul[role="listbox"] li[aria-selected="true"], div[role="listbox"] [role="option"][aria-selected="true"] {{
        background:#eef3f9 !important; color:{C_OXFORD} !important; }}
  div[data-baseweb="popover"] ul, div[data-baseweb="popover"] div[role="listbox"] {{
        background:#ffffff !important; border:1px solid #cdd8e4 !important; }}
  /* Selected-value chips in the multiselects. Brass fill with near-black text was left over
     from when the fields were dark; on the white field it read as a highlighter marking.
     Now a pale tint with Oxford-blue text, plus a border so the chip boundary still reads
     against white (tint alone separates only 1.12:1). Explicit padding and overflow:visible
     because the label was clipping its first characters ("ultiple sclerosis", "rong"). */
  section[data-testid="stSidebar"] [data-baseweb="tag"] {{
        background:#eef3f9 !important; border:1px solid #cdd8e4 !important;
        border-radius:2px !important; color:{C_OXFORD} !important;
        padding:1px 5px 1px 7px !important; margin:2px 4px 2px 0 !important;
        max-width:none !important; overflow:visible !important; }}
  section[data-testid="stSidebar"] [data-baseweb="tag"] span,
  section[data-testid="stSidebar"] [data-baseweb="tag"] div {{
        color:{C_OXFORD} !important; -webkit-text-fill-color:{C_OXFORD} !important;
        overflow:visible !important; text-overflow:clip !important;
        max-width:none !important; font-size:12px !important; }}
  section[data-testid="stSidebar"] [data-baseweb="tag"] svg {{
        fill:{C_OXFORD_LT} !important; }}
  section[data-testid="stSidebar"] [data-baseweb="tag"]:hover svg {{
        fill:{C_DRIVER} !important; }}
  /* ---------- tabs ---------- */
  .stTabs [data-baseweb="tab-list"] {{ gap:2px; border-bottom:1px solid {C_LINE}; }}
  .stTabs [data-baseweb="tab"] {{ height:41px; padding:0 17px; background:transparent;
        font-size:13px; font-weight:500; color:{C_MUTED}; border-radius:0; }}
  .stTabs [aria-selected="true"] {{ color:{C_OXFORD} !important; font-weight:600;
        border-bottom:2px solid {C_BRASS} !important; }}
  /* ---------- panels, tables, inputs ---------- */
  [data-testid="stExpander"] {{ border:1px solid {C_LINE}; border-radius:3px;
        background:{C_PANEL}; box-shadow:0 1px 2px rgba(0,33,71,0.04); }}
  [data-testid="stExpander"] summary {{ font-size:13px; font-weight:600; color:{C_OXFORD}; }}
  [data-testid="stDataFrame"] {{ border:1px solid {C_LINE}; border-radius:3px; }}
  .stButton > button {{ background:{C_OXFORD}; color:#fff; border:0; border-radius:2px;
        font-weight:600; font-size:13px; padding:9px 20px; letter-spacing:0.3px;
        transition:background 150ms ease; }}
  .stButton > button:hover {{ background:{C_OXFORD_LT}; }}
  div[data-testid="stTextInput"] input {{ border-radius:2px; border:1px solid #cdd8e4;
        font-size:13.5px; padding:11px 13px; }}
  div[data-testid="stTextInput"] input:focus {{ border-color:{C_OXFORD};
        box-shadow:0 0 0 2px rgba(0,33,71,0.10); }}
  /* Callouts: a single Oxford keyline instead of Streamlit's saturated fills. */
  div[data-testid="stAlert"] {{ border-radius:2px; border-left:3px solid {C_OXFORD};
        background:#eef3f9; font-size:13px; }}
  [data-testid="stCaptionContainer"] {{ color:{C_MUTED}; font-size:11.5px; }}
  #MainMenu, footer {{ visibility:hidden; }}
</style>
""", unsafe_allow_html=True)

df = load_shortlist()
client = get_client()

def section(num, title, sub=""):
    """Numbered section heading.

    A dense single-page app needs visible structure; numbering the sections gives the
    reader a sense of position and lets the demo script refer to "section 3" out loud.
    """
    st.markdown(
        f"<div class='sec'><span class='n'>{num}</span><h3>{title}</h3>"
        f"<span class='sub'>{sub}</span></div>", unsafe_allow_html=True)


# Masthead. The emoji is dropped from the title -- at 35px in a serif face the
# wordmark carries itself, and the eyebrow line does the categorising instead.
st.markdown("""
<div class="hero">
  <p class="eyebrow">Directional target nomination &middot; CD4&#8314; T-cell Perturb-seq</p>
  <h1>GenePerturb2Target</h1>
  <p class="byline">Oishi Deb &nbsp;&middot; Yizhou Yu</p>
  <div class="rule"></div>
  <p><b>1,923 Directional Nominations </b> from a genome-scale CRISPRi screen &mdash; each one a
  call to <b>block a driver</b> or <b>activate a brake</b>, anchored in human genetics and
  filtered for druggability. Every figure and filter is computed on real dataset the
  language model only translates your words into filters and narrates a gene's own
  evidence.</p>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.header("Filters")
    spec = {}

    # --- Disease first: it is the filter most users reach for, and it is the one that
    # --- most sharply changes the result set, so it leads the panel.
    st.markdown("### 🩺 Disease")
    dis_groups = st.multiselect(
        "Disease anchor (Open Targets)", list(DISEASE_GROUPS.keys()),
        help="Filters on each target's Open Targets disease association. Every disease "
             "term present in the data belongs to exactly one group here, so nothing is "
             "unreachable. Groups collapse synonyms (Crohn / ulcerative colitis / colitis "
             "all sit under inflammatory bowel disease). Autoimmune and allergic groups "
             "are listed first because they are this project's focus; the remaining "
             "groups exist because the underlying annotation contains them.")
    if dis_groups: spec["disease_groups"] = dis_groups
    if MS_ANCHORED_GENES:
        if st.checkbox(f"MS-anchored genes only ({len(MS_ANCHORED_GENES)})",
                       help="Genes carrying multiple sclerosis evidence from the "
                            "patient-derived signature analysis. This is a separate, "
                            "non-overlapping line of evidence from the Open Targets "
                            "MS annotation above."):
            spec["ms_anchored_only"] = True
    disease = st.text_input("…or type a disease",
                            help="Accepts abbreviations and categories: MS, RA, T1D, IBD, "
                                 "SLE, T2D, neurodegenerative, leukaemia, cancer.")
    if disease.strip(): spec["disease_contains"] = disease.strip()

    st.markdown("---")
    st.markdown("### Direction & biology")
    # Options read from DIRECTION_LABEL so the filter, the table and the figure legends
    # cannot drift apart. Matching is by dict lookup, not string prefix.
    _dopts = {DIRECTION_LABEL[k]: k for k in ("driver_antagonize", "brake_agonize")}
    dsel = st.selectbox("Direction", ["(any)"] + list(_dopts))
    if dsel in _dopts: spec["direction"] = _dopts[dsel]
    classes = st.multiselect("Protein class", sorted(df.primary_class.dropna().unique()))
    if classes: spec["primary_class"] = classes
    tiers = st.multiselect("Genetics tier", ["strong", "moderate", "weak"])
    if tiers: spec["genetics_tier"] = tiers
    cond = st.selectbox("Peak condition", ["(any)", "Rest", "Stim8hr", "Stim48hr"])
    if cond != "(any)": spec["condition"] = cond
    if st.checkbox("Druggable only"): spec["druggable_only"] = True
    if st.checkbox("Novel / undrugged only"): spec["novel_only"] = True
    maxrank = st.number_input("Max rank (0 = no limit)", min_value=0, max_value=1923, value=0)
    if maxrank > 0: spec["max_rank"] = int(maxrank)

    st.markdown("---")
    st.caption("Structures load from the bundled set (offline); missing ones fall back to AlphaFold DB. "
               "Natural-language box needs an ANTHROPIC_API_KEY")

# natural-language query box
section("I", "Ask in plain language", "natural-language filter translation")
if client is None:
    with st.expander("💬 Natural-language search is off — how to enable it", expanded=False):
        st.markdown(
            "Every filter and figure in this app works without an API key. "
            "The key only powers the optional plain-language search box and the "
            "per-gene explanations.\n\n"
            "**On Streamlit Community Cloud** (this deployment): open your app at "
            "share.streamlit.io → **⋮ → Settings → Secrets**, then add\n"
            "```toml\n"
            'ANTHROPIC_API_KEY = "sk-ant-..."\n'
            "```\n"
            "Save; the app reboots automatically with the key available.\n\n"
            "**Running locally** — either export it in the shell first:\n"
            "```bash\n"
            'export ANTHROPIC_API_KEY="sk-ant-..."\n'
            "streamlit run app.py\n"
            "```\n"
            "or create `.streamlit/secrets.toml` next to `app.py` with the same "
            "`ANTHROPIC_API_KEY = \"sk-ant-...\"` line. Keep that file out of git.\n\n"
            "Get a key from console.anthropic.com → API keys."
        )
# Worked examples. Each was run through the parser and the filter, and the count
# shown is the number of targets it actually returns -- including the zero.
EXAMPLE_QUERIES = [
    ("Novel undrugged targets for type 1 diabetes", "4 targets"),
    ("Which multiple sclerosis targets should I block?", "7 targets"),
    ("IBD drivers with strong genetics that are druggable", "8 targets"),
    ("What should I activate in psoriasis?", "3 targets"),
    ("Brake-agonize kinases with strong asthma genetics and a clean patent space", "1 target"),
    ("Show me transcription factors that act as brakes at rest", "13 targets"),
    ("Top 10 druggable rheumatoid arthritis targets", "2 targets"),
    ("Neurodegenerative disease targets in this screen", "38 targets"),
    ("Kinases and GPCRs I could inhibit in lupus", "0 — no SLE-anchored kinase/GPCR"),
]

nlq = st.text_input(
    "e.g. 'brake-agonize kinases with strong asthma genetics and a clean patent space'",
    disabled=(client is None), label_visibility="collapsed")

with st.expander("💡 Example questions you can ask", expanded=False):
    st.caption(
        "The parser understands four things at once: **direction** (block a driver / "
        "activate a brake), **disease** (19 groups, abbreviations fine), **protein class** "
        "(kinase, GPCR, transcription factor, cytokine receptor, enzyme, transporter, ion "
        "channel, nuclear receptor), and **context** (Rest / Stim8hr / Stim48hr) — plus "
        "'novel', 'druggable', 'strong genetics' and 'top N'. Counts below are what each "
        "query actually returns on this shortlist."
    )
    for _q, _n in EXAMPLE_QUERIES:
        _c1, _c2 = st.columns([5, 1])
        _c1.markdown(f"“{_q}”")
        _c2.caption(_n)
    st.caption(
        "The last one returns nothing on purpose: no lupus-anchored target in this "
        "shortlist is a kinase or GPCR. An empty result is a real answer, not a failure."
    )
if nlq and client is not None:
    try:
        nlspec = nl_to_spec(client, nlq)
        spec = nlspec  # NL query overrides sidebar when present
        chips = " ".join(f"<span class='pill' style='background:#2a4a6a'>{k}={v}</span>"
                         for k, v in nlspec.items())
        st.markdown("**Parsed filters:** " + (chips or "<i>none</i>"), unsafe_allow_html=True)
        # #5 query trace: matches shown highlighted on load; animation is an optional replay
        match_res = apply_filters(df, nlspec)
        if 0 < len(match_res) <= df.shape[0]:
            st.plotly_chart(query_trace_map(df, set(match_res.target_gene)),
                            width="stretch", config={"displayModeBar": False})
            st.caption("Your matches are highlighted against the full landscape; everything else is "
                       "faded back. Press ↻ Replay reveal to watch them light up, or ◌ Show all "
                       "to see the whole shortlist undimmed.")
    except Exception as e:
        _emsg = str(e)
        if "not_found_error" in _emsg or "model:" in _emsg:
            st.error(
                f"The configured model is not available to this API key "
                f"(`{MODEL}`, {MODEL_SOURCE}). Model ids are retired over time. "
                "Either remove `ANTHROPIC_MODEL` from your secrets so the app "
                "auto-selects a valid model, or set it to one your key supports. "
                "Sidebar filters are unaffected and still give the full result set."
            )
        elif "authentication" in _emsg.lower() or "401" in _emsg:
            st.error("The API key was rejected. Check `ANTHROPIC_API_KEY` in "
                     "Settings → Secrets. Sidebar filters still work.")
        elif "rate_limit" in _emsg or "429" in _emsg:
            st.warning("Rate-limited by the API — try again shortly. "
                       "Sidebar filters are unaffected.")
        else:
            st.error(f"Parse failed ({e}); using sidebar filters.")

# results results 
res = apply_filters(df, spec)

# Disease filters are informative only where an anchor exists, so say so explicitly
# rather than letting an empty table look like "no such targets".
_dgroups = spec.get("disease_groups") or ([] if not spec.get("disease_contains") else None)
if spec.get("disease_groups") or spec.get("disease_contains") or spec.get("ms_anchored_only"):
    _n_anchored = int((df.ot_top_disease.notna() | df.ot_autoimmune_disease.notna()).sum())
    _bits = []
    if spec.get("disease_groups"):
        _bits.append("disease anchor: " + ", ".join(spec["disease_groups"]))
    if spec.get("ms_anchored_only"):
        _bits.append(f"MS-anchored set ({len(MS_ANCHORED_GENES)} genes)")
    if spec.get("disease_contains"):
        _g = resolve_disease_text(spec["disease_contains"])
        _bits.append(f'text "{spec["disease_contains"]}"' + (f" -> {_g}" if _g else " (no group match)"))
    st.info(
        f"**Filtered by disease** — {' · '.join(_bits)}. "
        f"{len(res)} target(s). Note only {_n_anchored} of {len(df)} shortlist entries carry any "
        f"Open Targets disease anchor, so a disease filter narrows the list sharply; an empty "
        f"result means no *anchored* target matched, not that no target is relevant."
    )

section("II", "The shortlist", "filtered nominations and how they distribute")

# KPI cards KPI CARDS
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
(tab_map, tab_ctx, tab_funnel, tab_land, tab_ms, tab_overlay,
 tab_p2d, tab_weights, tab_bench) = st.tabs(
    ["🗺 Directional map", "⏱ Context dynamics", "🔻 Method funnel", "🌅 Landscape",
     "🧭 MS generalization", "🧬 Patient manifold", "🕸 Path2Drug",
     "⚖️ Re-weight the score", "📊 Benchmark"])
with tab_map:
    left, right = st.columns([1.15, 1])
    with left:
        st.plotly_chart(directional_map(res), width="stretch", config={"displayModeBar": False})
    with right:
        st.markdown(f"**{n} matching nominations**")
        cols = ["rank", "target_gene", "direction", "integrated_score", "genetics_tier",
                "ot_top_disease", "suggested_modality", "novelty_class"]
        st.dataframe(for_display(res, cols).reset_index(drop=True),
                     width="stretch", height=420)
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
    st.caption("The genome-scale screen narrows from 11,526 genes to six deep-dive leads. Each stage "
               "is a strict subset of the one above it. Press ↻ Replay to watch the filters apply "
               "one at a time.")
with tab_land:
    st.plotly_chart(disease_sunburst(res), width="stretch", config={"displayModeBar": False})
    st.caption("Click a wedge to zoom in: direction → protein class → disease anchor.")

# per-gene deep dive: evidence bar + live 3D structure + grounded explanation
if n:
    st.markdown("---")
    section("III", "Nomination deep-dive", "per-target evidence, structure and mechanism")
    gene = st.selectbox("Select a gene", res.target_gene.tolist())
    row = res[res.target_gene == gene].iloc[0]
    dircol = C_DRIVER if row.direction == "driver_antagonize" else C_BRAKE
    action = "BLOCK (antagonize a driver)" if row.direction == "driver_antagonize" else "ACTIVATE (agonize a brake)"
    st.markdown(f"### {gene} &nbsp; <span class='pill' style='background:{dircol}'>{action}</span>",
                unsafe_allow_html=True)

    # gauges row (#6)
    st.plotly_chart(gauge_row(row), width="stretch", config={"displayModeBar": False})

    # score buildup (#4) — shown fully populated; the animation is an optional replay
    st.plotly_chart(score_buildup(row), width="stretch", config={"displayModeBar": False})
    st.caption("The four weighted components (causal 0.34, genetics 0.30, druggability 0.22, "
               "novelty 0.14) sum to the integrated score. Press ↻ Replay build to watch them "
               "stack up one at a time.")

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
                    st.warning(f"Structure not available for {gene}.")
        else:
            st.info("3D viewer hidden — tick 'Show 3D structure' to display it.")
    with cC:
        st.caption("Plain-language rationale — generated only from this gene's evidence row")
        if client is None:
            st.info("Add an ANTHROPIC_API_KEY (see the note under the search box) for a grounded explanation.")
        elif st.button("💬 Explain with Claude", key="explain"):
            with st.spinner("Generating from the gene's evidence…"):
                st.markdown(explain_gene(client, row))
        with st.expander("Raw evidence row"):
            st.json({"integrated_score": round(float(row.integrated_score), 3),
                     "peak_emp_z": round(float(row.peak_emp_z), 2),
                     "peak_emp_fdr": float(row.peak_emp_fdr),
                     "genetics_tier": row.genetics_tier, "top_disease": row.ot_top_disease,
                     "suggested_modality": row.suggested_modality, "novelty_class": row.novelty_class})

    # agentic reasoning layer — shown only for the 92 genes actually triaged
    _ag = load_agentic_calls()
    if _ag is not None and gene in set(_ag.gene):
        _call = _ag[_ag.gene == gene].iloc[0]
        _agree = str(_call["agrees_with_screen"])
        _badge = {"yes": ("agrees with the screen", C_AGREE),
                  "no": ("disagrees with the screen", C_DISAGREE),
                  "uncertain": ("uncertain", "#8b9aa8")}.get(_agree, ("uncertain", "#8b9aa8"))
        st.markdown(
            "<div class='sec' style='margin-top:22px'><span class='n'>&mdash;</span>"
            "<h3 style='font-size:15px'>Agentic reasoning layer</h3>"
            "<span class='sub'>an independent mechanistic read of this gene</span></div>",
            unsafe_allow_html=True)
        _r1, _r2 = st.columns([1.35, 1])
        with _r1:
            st.markdown(
                f"<span class='pill' style='background:{_badge[1]}'>{_badge[0]}</span> "
                f"<span class='pill' style='background:#5a6b7b'>{_call['confidence']} confidence</span>",
                unsafe_allow_html=True)
            _act = _call["therapeutic_action"]
            st.markdown(
                f"- **screen says** {_call['screen_direction']}\n"
                f"- **reasoning layer says** "
                f"{_act if isinstance(_act, str) and _act.strip() else '_no call_'}\n"
                f"- **reference rank** {int(_call['ref_rank'])}")
            st.markdown(f"**Mechanistic rationale.** {_call['mechanistic_rationale']}")
            _inc = _call["primary_inconsistency"]
            if isinstance(_inc, str) and _inc.strip():
                st.warning(f"**Flagged inconsistency.** {_inc}")
            _rep = _call["repurposing_note"]
            if isinstance(_rep, str) and _rep.strip():
                st.caption(f"Repurposing note: {_rep}")
        with _r2:
            st.plotly_chart(agentic_layers_fig(_call), width="stretch",
                            config={"displayModeBar": False})
            st.caption(
                "The reasoning layer reads each gene's evidence independently of the "
                "ranking. Where it disagrees with the screen's directional call, that "
                "disagreement is the signal worth attention — not an error to hide."
            )
        with st.expander("How this gene sits among all 92 triaged genes"):
            st.plotly_chart(agentic_overview_fig(_ag), width="stretch",
                            config={"displayModeBar": False})
            _nd = int((_ag.agrees_with_screen == "no").sum())
            st.caption(
                f"{_nd} of 92 genes disagree with the screen's call, and every one carries an "
                "explicit inconsistency flag rather than a silent override. Only 92 of the "
                "1,923 shortlisted genes were triaged, so absence of a panel means a gene "
                "was not assessed — not that it passed."
            )

    # signature-response heatmap (#1) — the mechanism reveal
    sigdf = load_signature_response()
    if sigdf is not None and gene in set(sigdf.target_gene):
        st.markdown("<div class='sec' style='margin-top:22px'><span class='n'>&mdash;</span><h3 style='font-size:15px'>Why this call? Signature-gene response</h3></div>", unsafe_allow_html=True)
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
            st.dataframe(for_display(show.round(2), cols).reset_index(drop=True),
                     width="stretch", height=340)
            st.caption("Druggable genetic-anchored leads: **IL2RB, TYK2, IL2RA, IL7R, TNFRSF1A** — "
                       "TYK2 is a validated approved-drug-class target.")
    with cR:
        st.markdown("**Drug-mechanism concordance** — direction vs how the approved drug acts")
        if ms_cc is not None:
            ccols = [c for c in ["gene", "drug", "drug_action", "our_call", "concordant"]
                     if c in ms_cc.columns]
            st.dataframe(for_display(ms_cc[ms_cc.concordant.notna()], ccols).reset_index(drop=True)
                         if "concordant" in ms_cc else for_display(ms_cc, ccols),
                         width="stretch", height=340)
            st.caption("Natalizumab/ITGA4 · fingolimod/S1PR1 · alemtuzumab/CD52 · TYK2 inhibitors · "
                       "anti-IL2RA — the corrected direction calls all five as *antagonize*.")

    st.info("**The thesis in miniature on new data:** undirected, uncorrected reversal recovers "
            "pleiotropic housekeeping hubs; directional scoring intersected with genetics recovers "
            "specific, druggable, correctly-directed targets. The framework applies to any "
            "CD4-T-cell-mediated disease with a definable signature.", icon="🧭")

with tab_overlay:
    _proj = load_ms_projection()
    _bg = load_ms_background()
    if _proj is None:
        st.info("Projection data (ms_projection_app.parquet) not found next to app.py app.py.")
    else:
        _c1, _c2 = st.columns([1, 3])
        with _c1:
            _ctx = st.selectbox("Activation context", ["All", "Rest", "Stim8hr", "Stim48hr"],
                                key="ovl_ctx",
                                help="The static figure pooled all three contexts. Selecting one "
                                     "shows that a knockdown's position on the manifold moves as "
                                     "the T cell activates.")
            _topn = st.slider("Label top N reversing", 5, 40, 20, 5, key="ovl_top")
        st.plotly_chart(ms_manifold_overlay(_proj, _bg, _ctx, _topn),
                        width="stretch", config={"displayModeBar": False})
        st.caption(
            "Each coloured point is one knockdown in one activation context, placed where its "
            "transcriptional effect lands in the integrated patient CD4 manifold (188,422 cells, "
            "25,000 drawn as the grey cloud). Red = the knockdown moves cells against the MS "
            "disease direction (candidate therapeutic); blue = with it. Ringed points are the "
            "top reversing knockdowns for the selected context."
        )
        if _bg is None:
            st.caption("Patient-cell background not found — showing knockdowns only.")
        st.info(
            "**Read this panel with the caveat from the MS generalization tab.** Raw reversal "
            "score is confounded by how many genes a knockdown perturbs: broad knockdowns score "
            "well simply by moving cells a long way. The footprint-matched correction reduces the "
            "confound correlation from −0.33 to −0.07, and the directional result quoted in the "
            "paper uses the corrected score, not the raw one shown here."
        )

with tab_p2d:
    _p2d_tfs = p2d_available()
    _p2d_sum = load_p2d_summary()
    if not _p2d_tfs or _p2d_sum is None:
        st.info("Path2Drug outputs not found. Copy `path2drug/` and the "
                "`path2drug_*.csv` tables into the app folder.")
    else:
        st.markdown(
            "**Thirteen of the top-100 nominations are transcription factors with no "
            "small-molecule pocket.** A directional call on an undruggable gene is not yet "
            "actionable. `Path2Drug` walks confidence-weighted STRING paths from the TF to "
            "the program's signature genes and asks which node *on that path* is druggable "
            "— turning an undruggable nomination into a tractable intervention point. "
            "Ten of those thirteen have a released backbone below; ANKZF1, ZNF236 and "
            "ZNF438 have too little STRING support to path from.")
        _c1, _c2 = st.columns([1, 1.6])
        with _c1:
            _tf = st.selectbox("Transcription factor", _p2d_tfs,
                               index=_p2d_tfs.index("STAT4") if "STAT4" in _p2d_tfs else 0)
        _row = _p2d_sum[_p2d_sum.tf == _tf]
        with _c2:
            if len(_row):
                _r = _row.iloc[0]
                st.markdown(
                    f"<div style='padding-top:26px'>"
                    f"<span class='pill' style='background:{C_OXFORD}'>{_r.direction}</span> "
                    f"&nbsp;{int(_r.n_intermediates)} intermediates &middot; "
                    f"<b>{int(_r.n_druggable)} druggable</b> &middot; "
                    f"{int(_r.n_dir_consistent)} direction-consistent &middot; "
                    f"{int(_r.n_edge_confirmed)} edge-confirmed</div>",
                    unsafe_allow_html=True)
        _d = load_p2d_target(_tf)
        if _d is None:
            st.info(f"No released backbone for {_tf}.")
        else:
            _fig = p2d_backbone_fig(_d)
            if _fig is not None:
                st.plotly_chart(_fig, width="stretch", config={"displayModeBar": False})
            _drug = [n for n in _d.get("intermediate_nodes", []) if n.get("druggable")]
            if _drug:
                _bits = []
                for n in _drug:
                    _dr = n["druggable"]
                    _eff = n.get("peak_signed_score")
                    _bits.append(
                        f"**{n['gene']}** (ChEMBL phase {_dr.get('max_phase','?')}, "
                        f"{str(_dr.get('action','?')).lower()})"
                        + (f", own program effect {_eff:+.2f}"
                           f"{', direction-consistent' if n.get('direction_consistent') else ''}"
                           if _eff is not None else ""))
                st.markdown("**Druggable intermediates on this backbone:** " + "; ".join(_bits) + ".")
            else:
                st.caption(f"No druggable intermediate on {_tf}'s backbone — "
                           "one of the ten TFs where the method returns nothing actionable.")
        st.markdown("---")
        _cA, _cB = st.columns(2)
        with _cA:
            _bl = load_p2d_baseline()
            if _bl is not None:
                st.plotly_chart(p2d_baseline_fig(_bl), width="stretch",
                                config={"displayModeBar": False})
                st.caption(
                    "The control that makes this a result rather than a diagram: a naive "
                    "highest-degree-druggable-neighbour heuristic collapses to the generic "
                    "hub CD4 for almost every TF. Path2Drug returns pathway-specific nodes. "
                    "Of the seven selected nodes with an observed program effect, 5 (71%)"
                    "move the program in the direction the TF's own call implies.")
        with _cB:
            _ec = load_p2d_edges()
            if _ec is not None:
                st.plotly_chart(p2d_edge_fig(_ec), width="stretch",
                                config={"displayModeBar": False})
                st.caption(
                    "STRING edges are literature co-occurrence, so each target→intermediate "
                    "edge is re-tested against the atlas: does knocking the TF down actually "
                    "move the intermediate's transcript? Of the 30 edges, 24 are testable "
                    "(both genes measured) and only 2 clear |z|≥1.96. The "
                    "backbone is a hypothesis generator, not a validated circuit.")
        st.markdown("---")
        st.markdown("**Where this sits in the argument.** The flagship case is "
                    "STAT4 → IL12RB2 → **JAK2** → IFNG: JAK2 has approved inhibitors, and its "
                    "own knockdown lowers the program (−0.42). That recovers — rather than "
                    "assumes — the established logic of treating STAT-driven inflammation "
                    "with JAK inhibitors, which is the point of the control: a method that "
                    "re-derives known pharmacology from perturbation data alone is one you "
                    "can trust a little further on the cases where no drug exists yet.")
        st.caption("Scope: the 10 hard-to-drug TFs above, none of which is among the six "
                   "deep-dive leads — this is a separate, complementary output. Any "
                   "language-model narration in the pipeline is constrained to the extracted "
                   "subnetwork and cannot introduce genes or edges; the deterministic "
                   "backbone, not the prose, is the deliverable.")

with tab_weights:
    _ws = load_weight_sensitivity()
    _intro = ("**The published integration weights are a choice, not a result.** Move them "
              "and watch the ranking respond.")
    if _ws:
        _intro += (" The panel is meant to be hard to break: under ±40% random perturbation "
                   f"the ranking holds at Spearman {_ws['spearman_median']:.3f}.")
    st.markdown(_intro)
    if _ws:
        st.caption(
            f"({_ws['perturbation']}; 95% CI {_ws['spearman_ci'][0]:.3f}–{_ws['spearman_ci'][1]:.3f}). "
            f"Top-20 membership is less stable — median Jaccard {_ws['top20_jaccard_median']:.2f} — "
            "so individual positions move even though the overall order does not.")
    _w1, _w2 = st.columns([1, 2.6])
    with _w1:
        st.markdown("**Weights**")
        _wc = st.slider("Causal", 0.0, 1.0, 0.34, 0.02, key="w_causal")
        _wg = st.slider("Genetics", 0.0, 1.0, 0.30, 0.02, key="w_gen")
        _wd = st.slider("Druggability", 0.0, 1.0, 0.22, 0.02, key="w_drug")
        _wn = st.slider("Novelty", 0.0, 1.0, 0.14, 0.02, key="w_nov")
        _wsum = _wc + _wg + _wd + _wn
        st.caption(f"Sum {_wsum:.2f} — normalised to 1 before scoring, so the score stays "
                   "on its published scale.")
        _wtop = st.slider("Show top N", 5, 30, 15, 5, key="w_topn")
        if st.button("↺ Reset to published", key="w_reset"):
            for _k, _v in [("w_causal", 0.34), ("w_gen", 0.30),
                           ("w_drug", 0.22), ("w_nov", 0.14)]:
                st.session_state[_k] = _v
            st.rerun()
    with _w2:
        _wfig, _wd_df = weight_reranking_fig(df, _wc, _wg, _wd, _wn, _wtop)
        st.plotly_chart(_wfig, width="stretch", config={"displayModeBar": False})
    if _ws:
        _dor = _ws["drop_one_retention"]
        st.markdown(
            "**Which layer is load-bearing?** Dropping each component and measuring how much "
            "of the top-50 survives: "
            + " · ".join(f"**{k}** {v:.0%}" for k, v in _dor.items())
            + ". Genetics matters most, novelty least — so the ranking is not simply a "
              "novelty filter in disguise.")

with tab_bench:
    _mb = load_multibaseline()
    if _mb is None:
        st.info("Baseline comparison data (multibaseline_comparison_results.csv) not found.")
    else:
        st.plotly_chart(multibaseline_fig(_mb), width="stretch",
                        config={"displayModeBar": False})
        st.markdown(
            "**This is the result we did not want.** Asked to recover the 19 known "
            "immunomodulatory drug targets in the shortlist, **network centrality (0.950) "
            "and Open Targets genetics (0.909) both beat our full integrated score "
            "(0.820)**. We report it because a target-prioritisation method that only "
            "shows the benchmarks it wins is not evidence of anything."
        )
        st.info(
            "**Why the method earns its place.** Recovering already-known "
            "targets is a test of how well a score reproduces existing knowledge — and a "
            "centrality measure trained on the literature's own citation structure should "
            "win it. It says nothing about *direction*: none of these baselines states "
            "whether to block or activate a target, which is the call this method makes. "
            "The CIs also overlap heavily (ours 0.695–0.926 vs centrality 0.928–0.971), so "
            "with 19 positives this ordering is not firmly established either way."
        )
        with st.expander("Full numbers"):
            st.dataframe(
                for_display(_mb, ["axis", "group", "auroc", "auroc_lo", "auroc_hi", "ap",
                                  "n_positives", "n_universe"]).reset_index(drop=True),
                width="stretch")

st.markdown("---")
st.caption("All nominations are computational hypotheses; see PROSPECTIVE_VALIDATION.md for the "
           "pre-registered falsification protocol. Directional map axes are the pipeline's causal and "
           "genetics score components; bubble size is the integrated score.")
