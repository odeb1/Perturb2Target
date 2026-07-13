# Perturb2Target — interactive directional target explorer

A self-contained web app over the **1,923-target directional shortlist** from the
genome-scale CD4⁺ T-cell CRISPRi Perturb-seq nomination pipeline. Built for a live
demo: a judge can type a plain-language query and watch the real shortlist filter
in front of them.

## What makes it trustworthy (the design that matters)

### Visuals (what a judge sees)

- **KPI cards** — live counts (nominations / novel-undrugged / strong-genetics / druggable) that update with every filter.
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

## Run

```bash
pip install -r requirements.txt
export ANTHROPIC_API_KEY=sk-...     # OPTIONAL — app works fully without it
streamlit run app.py
```

Then open the URL Streamlit prints (usually http://localhost:8501).
Set `ANTHROPIC_MODEL` to override the default (`claude-3-5-haiku-latest`).

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
