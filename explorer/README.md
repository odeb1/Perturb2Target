# Perturb2Target — interactive directional target explorer

A self-contained web app over the **1,923-target directional shortlist** from the
genome-scale CD4⁺ T-cell CRISPRi Perturb-seq nomination pipeline. Built for a live
demo: a judge can type a plain-language query and watch the real shortlist filter
in front of them.

## What makes it trustworthy (the design that matters)

### Visuals (what a judge sees)

- **KPI cards** — live counts (nominations / novel-undrugged / strong-genetics / druggable) that update with every filter.
- **Disease filter (top of the sidebar)** — a pick-list of **19 disease groups** covering **every** disease term present in the annotation, so nothing is unreachable: all 200 shortlist entries that carry an Open Targets disease association are selectable. Autoimmune and allergic groups lead the list (MS, RA, T1D, IBD, asthma/allergy, SLE, psoriasis, autoimmune thyroid, ankylosing spondylitis, vitiligo/alopecia, other autoimmune); the remaining groups exist because the annotation contains them — primary immunodeficiency, neurodegenerative disease (38 genes, the largest single category), neuropsychiatric, haematological malignancy, solid tumour, metabolic, cardiovascular, and a rare/Mendelian catch-all. Groups collapse synonymous terms, so *Crohn disease*, *ulcerative colitis* and *colitis* all sit under inflammatory bowel disease.
  A separate checkbox restricts to the 11 MS genes from the patient-derived signature analysis — a non-overlapping line of evidence from the Open Targets MS annotation. The free-text box resolves abbreviations and categories (MS, RA, T1D, IBD, SLE, T2D, neurodegenerative, Alzheimer, leukaemia, cancer, immunodeficiency) to the same groups; unrecognised text falls back to a word-boundary match, so "RA" no longer matches inside *random*.
  Disease coverage is sparse by nature: only 200 of 1,923 shortlist entries carry any scored disease association, so selecting a disease narrows the list sharply. The app states this inline whenever a disease filter is active, so an empty table reads as "no *anchored* target matched" rather than "no target is relevant".
- **Directional map** (hero tab) — interactive Plotly scatter of every matching target: x = causal directional strength, y = human-genetics support, bubble size = integrated score, colour = direction (red = block a driver, teal = activate a brake). Hover for gene, rank, disease.
- **Context dynamics** (hero tab) — the map *animated* across Rest → Stim8hr → Stim48hr with a ▶ Play button, so targets visibly move as the T cell activates — the project's context-specificity claim, live.
- **Method funnel** (hero tab) — an *animated* cascade: 11,526 genes screened → 1,923 with directional signal → 374 druggable → 150 genetics-anchored → 86 novel → 6 leads. Press ▶ to watch the genome narrow to the shortlist, one filter at a time (each stage a strict subset of the one above).
- **Landscape** (hero tab) — a clickable sunburst: direction → protein class → disease anchor.
- **MS generalization** (hero tab) — the whole pipeline re-run, unchanged, on an independent patient-derived multiple-sclerosis signature (1.9M CD4⁺ cells, two cohorts): KPI cards (86% of MS genetic anchors correctly-directed vs 58% background, p=0.007; 5/5 approved-MS-drug mechanism concordance; pleiotropy confound removed ρ −0.33→−0.07; 0 knockdowns survive genome-wide FDR — the honest headline), the four-panel result figure, and the MS nomination + drug-concordance tables. Demonstrates the method is dataset-agnostic, not tuned to the curated inflammatory program.
- **Query trace** (natural-language box) — when a plain-language query returns matches, an *animated* trace lights up the matching targets while the rest of the directional landscape fades back, so the filter feels physical.
- **Signature-response heatmap** (deep-dive) — the *mechanism reveal*: the 48 signature genes' z-scores under the selected knockdown, split into pro-inflammatory and regulatory arms. Shows *why* the directional call was made. A **▶ Animate across conditions** toggle replays the response Rest → 8 h → 48 h. Bundled for the 49-gene demo set × 3 conditions.
- **3D structure** (deep-dive) — two modes: **ligand-in-pocket** (Boltz-2 co-folded complex, docked small molecule as green sticks, for the 11 pocket-bearing nominations) and **AlphaFold** (coloured by pLDDT). The structure **renders automatically** when a nomination is selected (no button click; untick *Show 3D structure* to hide it). Both views **auto-rotate slowly on load**. The ligand view uses a genuine **staged reveal**: the protein cartoon renders first, then after ~0.9 s the pre-computed ligand *appears* while the camera animates a 1.2 s zoom into the pocket. (The ligand is shown/hidden in its Boltz-predicted position — no binding motion is invented; it becomes visible where it was co-folded.) 49 AlphaFold models + 11 Boltz complexes pre-bundled → fully offline; genes outside the bundle fall back to a live AlphaFold fetch (stdlib `urllib`).
- **Score buildup** (deep-dive) — an *animated* stacked bar: the four weighted evidence components (causal 0.34 / genetics 0.30 / druggability 0.22 / novelty 0.14) build up to the integrated score.
- **Gauges** (deep-dive) — integrated / genetics / tractability dials.
- **Compare mode** (deep-dive) — 2–3 nominations side by side.
- **Evidence bar** — the four score components (causal / genetics / druggability / novelty) for the selected gene.
- **Grounded explanation** — one-click Claude rationale generated only from that gene's evidence row.

### Bundled data files (must sit next to `app.py`)
- `signature_response.parquet` — per-signature-gene z-scores for the heatmap
- `map_by_condition.parquet` — per-condition directional scores for the animation
- `structures/` (49 AlphaFold `.pdb`) · `structures_ligand/` (11 Boltz `.cif`)

### The two-layer trust design

The app has two layers, and the split is deliberate:

1. **Deterministic core (always on, no network).** Sidebar filters and the natural-language
   box both emit the *same* structured filter spec, and all filtering runs on the real
   `target_shortlist.csv`. This layer works with no API key and never changes the data.

2. **Optional Claude layer (only if `ANTHROPIC_API_KEY` is set).**
   - **Query parsing:** Claude translates your free-form question into the structured
     filter spec via a **tool schema** (constrained to the real column value spaces).
     It maps *language → filter parameters*; it never selects genes itself.
   - **Per-gene explanation:** Claude narrates a selected gene **strictly from that
     gene's own computed evidence row** (a JSON of its scores, genetics tier, disease
     anchor, modality). The system prompt forbids adding any fact, number, disease, or
     mechanism not present in that row.

So the LLM makes the tool *conversational* without ever inventing a gene or a number —
the same honesty discipline the whole project rests on.

## Design system

The interface uses **Oxford Blue `#002147`** as its single brand colour, with a muted
brass as the only accent. Tokens live at the top of `app.py`:

| Token | Value | Use |
|---|---|---|
| `C_OXFORD` | `#002147` | masthead, sidebar, headings, metric figures |
| `C_OXFORD_LT` | `#01305f` | gradient stop, hover states, links |
| `C_BRASS` | `#b08d57` | accents **on Oxford Blue** only (5.2:1 there) |
| `C_BRASS_DK` | `#8a6a38` | brass text on light backgrounds (`#b08d57` is 2.9:1 there — under the floor) |
| `C_DRIVER` / `C_BRAKE` | `#c1272d` / `#1a7f9e` | **data encoding** — block a driver / activate a brake |

`C_DRIVER` and `C_BRAKE` are deliberately outside the brand palette: they carry the
directional call, which is the paper's central claim, so they must stay distinguishable
from decoration and from each other. Do not restyle them for visual consistency.

Typography is Inter for UI and Source Serif 4 for display (the wordmark, metric figures
and section numerals), both loaded from Google Fonts with a full system fallback so the
app degrades cleanly offline.

Every foreground/background pair in the stylesheet was measured against WCAG AA: all 18
clear the 4.5:1 body-text threshold except the section numeral, which is large display
text and clears the 3:1 large-text floor at 4.7:1.

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then open the URL Streamlit prints (usually http://localhost:8501).

## API key (optional)

Every filter, figure and structure view works **without** a key. A key only enables
two optional extras: the plain-language search box and the per-gene explanations.

The app looks for `ANTHROPIC_API_KEY` in Streamlit secrets first, then in the
environment, so either of these works.

**Streamlit Community Cloud** — the hosted deployment. Open the app on
share.streamlit.io, then **⋮ → Settings → Secrets**, and add:

```toml
ANTHROPIC_API_KEY = "sk-ant-..."
```

Save. The app reboots automatically and picks the key up. Do not commit the key to
the repository; the Secrets panel is the only place it needs to exist.

**Locally** — either export it before launching:

```bash
export ANTHROPIC_API_KEY="sk-ant-..."
streamlit run app.py
```

or copy `.streamlit/secrets.toml.example` to `.streamlit/secrets.toml` and fill in
the key there. That path is gitignored, so a real key cannot be committed by accident.

**Model selection is automatic.** The app calls the models endpoint with your key and
picks the cheapest suitable model available to it (Haiku-class first, falling back to
Sonnet then Opus), so a retired model id cannot break the app. Set `ANTHROPIC_MODEL` in
secrets only if you want to pin a specific one — a pinned value always wins, and a stale
pin will produce a 404, so prefer leaving it unset.
Keys come from console.anthropic.com → API keys.

## Tabs

| Tab | Shows |
|---|---|
| 🗺 Directional map | causal effect vs genetics, bubble = integrated score |
| ⏱ Context dynamics | targets move as the T cell activates (Rest → 8 h → 48 h) |
| 🔻 Method funnel | 11,526 genes → 1,923 directional → 150 genetics-anchored |
| 🌅 Landscape | direction → protein class → disease sunburst |
| 🧭 MS generalization | confound-corrected directional test on MS patient data |
| 🧬 Patient manifold | candidate MS-reversing knockdowns overlaid on the patient CD4 manifold |
| ⚖️ Re-weight the score | move the four integration weights and watch the ranking respond |
| 📊 Honest benchmark | our scores vs three external baselines — including where we lose |

Plus a per-gene **Agentic reasoning layer** panel in the deep-dive section, for the 92
genes the reasoning layer triaged.

### Patient manifold tab

`ms_projection_app.parquet` (7,874 knockdown×context effects) placed on the integrated
CD4 atlas UMAP, over a grey cloud of patient cells from
`ms_manifold_background.parquet` — 25,000 of the atlas's 188,422 cells, downsampled
because the cloud only conveys the manifold's shape.

Colour is the **raw** MS reversal score on a diverging scale centred on zero: red moves
cells against the MS disease direction (candidate therapeutic), blue with it. The
activation context is selectable, which the static version of this figure could not do —
pooling contexts hid the fact that a knockdown's position moves as the T cell activates,
and the top-reversing set genuinely differs between them (EIF1/ATXN2L at rest,
KATNIP/MPI at 8 h, RPUSD3/MICOS13 at 48 h).

The panel carries an on-screen caveat: raw reversal score is confounded by knockdown
breadth, and the result quoted in the paper uses the footprint-matched corrected score
(confound correlation −0.33 → −0.07), not the raw score plotted here.

### Agentic reasoning layer panel

`agentic_triage_calls.csv` — 92 of the 1,923 shortlist genes (4.8%), including the
top-ranked ones. For each, the reasoning layer's own `therapeutic_action` call, a
confidence, `layers_concordant` (1–4), and whether it agrees with the screen.

**29 of 92 disagree with the screen's directional call, and all 29 carry an explicit
`primary_inconsistency`.** The panel shows disagreements rather than suppressing them —
IL4R is the sharp case: ranked #1 as a brake to agonize, while the reasoning layer notes
at high confidence that dupilumab works by *blocking* IL4R. A pipeline that surfaces its
own contradiction is more trustworthy than one that hides it.

The layer indicator shows filled/empty slots for `layers_concordant` and deliberately does
not name *which* layer dissented — the data is a count, and rendering per-layer detail
would imply resolution the file does not contain.

### Re-weight the score tab

Sliders for the four integration weights, recomputing `custom_score` from the stored
components and re-ranking live. Weights are normalised to sum to 1, so the score stays on
its published scale.

Verified: at the published weights (0.34/0.30/0.22/0.14) the panel reproduces the
published top-15 exactly (Pearson r = 1.0, identical ranks). The panel is responsive —
genetics-only weighting swaps 7 of the top 15.

The tab also reports the robustness analysis: Spearman 0.995 under ±40% perturbation, but
top-20 Jaccard only 0.74, and drop-one retention showing genetics most load-bearing (0.46)
and novelty least (0.84).

### Honest benchmark tab

`multibaseline_comparison_results.csv` — six methods recovering 19 known immunomodulatory
drug targets from the 1,923 shortlist, with bootstrap CIs.

**This tab shows a result where the method loses.** Network centrality (0.950) and Open
Targets genetics (0.909) both beat the full integrated score (0.820). The tab states why
we keep it: recovering known targets rewards reproducing existing knowledge, none of the
baselines makes a *directional* call, and the CIs overlap heavily (0.695–0.926 vs
0.928–0.971) so with 19 positives the ordering is not firmly established.

## What you can ask in plain language

The parser reads four dimensions at once — **direction**, **disease**, **protein class**
and **activation context** — plus `novel`, `druggable`, `strong genetics` and `top N`.
Every example below was run through the real parser and filter; the count is what it
actually returns on the 1,923-target shortlist.

| Question | Returns |
|---|---|
| "Novel undrugged targets for type 1 diabetes" | 4 — CTSH, ACAP1, CPEB3, SLC25A37 |
| "Which multiple sclerosis targets should I block?" | 7 — IL2RA, TBX21, STAT1, BATF … |
| "IBD drivers with strong genetics that are druggable" | 8 — IL10RB, ADAM17, PPP5C, DAGLB … |
| "What should I activate in psoriasis?" | 3 — ITGAL, TNFAIP3, PPIF |
| "Brake-agonize kinases with strong asthma genetics and a clean patent space" | 1 — SIK2 |
| "Show me transcription factors that act as brakes at rest" | 13 |
| "Top 10 druggable rheumatoid arthritis targets" | 2 — IL6R, CD2 |
| "Neurodegenerative disease targets in this screen" | 38 |
| "Targets that peak only after 48 hours of stimulation" | 603 |
| "Kinases and GPCRs I could inhibit in lupus" | **0** — no SLE-anchored kinase or GPCR exists |

The last row is deliberate: an empty result is a real answer. No lupus-anchored target in
this shortlist is a kinase or GPCR, and the app says so rather than silently showing nothing.

Vocabulary the parser maps: *antagonize / block / inhibit / driver* → block a driver;
*agonize / activate / potentiate / brake* → activate a brake; *novel / undrugged / clean
patent space* → novelty filter; disease abbreviations MS, RA, T1D, IBD, SLE, T2D and
category words (neurodegenerative, leukaemia, cancer, immunodeficiency) all resolve.

## Demo script (90 seconds)

1. **Show the whole shortlist** (1,923 rows) — "every nomination is directional: block a
   driver, or activate a brake."
2. **Type the flagship query:**
   `brake-agonize kinases with strong asthma genetics and a clean patent space`
   → returns **SIK2**, the lead nomination. Show the parsed filter chips.
3. **Point at the directional map** — the red/teal split IS the core idea; SIK2 sits in
   the high-genetics region and is ringed as the selection.
4. **Open the deep-dive**, click **"Load 3D structure"** — SIK2's AlphaFold model renders
   live, coloured by confidence. This is the visual wow moment.
5. **Click "Explain with Claude"** on SIK2 — the explanation is generated only from
   SIK2's evidence row; point out it cites the same numbers shown in the panel.
6. **Contrast honesty:** flip Direction to *driver_antagonize*, `max_rank ≤ 20` — the
   top rows are the positive controls (IL2RA, STAT3, IL10RB), i.e. the method rediscovers
   known biology at the top.

## Files
- `app.py` — the Streamlit app (deterministic core + optional Claude layer + visuals)
- `target_shortlist.csv` — the 1,923-target scored shortlist (the app's data)
- `structures/` — 49 pre-bundled AlphaFold `.pdb` models (deep-dive noms, Boltz pocket set, positive controls, top-40 by rank)
- `requirements.txt` — dependencies

## Troubleshooting

- **`segmentation fault` when loading a 3D model** — this is a native crash in some
  miniforge builds where `requests`/urllib3 hits the macOS system SSL library. The app
  now (a) reads bundled structures from `structures/` with no network call at all, and
  (b) uses stdlib `urllib` for any fallback fetch, so the crash path is avoided. If you
  still see it, the bundled genes (all deep-dive noms and controls) never touch the
  network — stick to those for the demo.
- **`numpy.dtype size changed` on startup** — broken numpy/pandas pair in the base env.
  Run in a clean venv: `python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt`.

All nominations are **computational hypotheses**; see `../PROSPECTIVE_VALIDATION.md` for
the pre-registered falsification protocol.
