"""
Segfault isolation test. Run OUTSIDE streamlit:
    cd /Users/oishideb/Downloads/Perturb2Target/explorer
    python diagnose.py
Whatever line prints LAST before the crash names the culprit.
"""
import sys, os
print("[0] python:", sys.version.split()[0], "| exe:", sys.executable); sys.stdout.flush()

print("[1] importing numpy..."); sys.stdout.flush()
import numpy as np
print("    numpy", np.__version__); sys.stdout.flush()

print("[2] importing pandas..."); sys.stdout.flush()
import pandas as pd
print("    pandas", pd.__version__); sys.stdout.flush()

print("[3] importing plotly..."); sys.stdout.flush()
import plotly; print("    plotly", plotly.__version__); sys.stdout.flush()

print("[4] importing py3Dmol..."); sys.stdout.flush()
import py3Dmol; print("    py3Dmol", getattr(py3Dmol, "__version__", "?")); sys.stdout.flush()

print("[5] reading bundled SIK2.pdb..."); sys.stdout.flush()
pdb = open(os.path.join(os.path.dirname(__file__), "structures", "SIK2.pdb")).read()
print("    read", len(pdb), "bytes"); sys.stdout.flush()

print("[6] building py3Dmol view + HTML..."); sys.stdout.flush()
v = py3Dmol.view(width=400, height=400)
v.addModel(pdb, "pdb")
v.setStyle({"cartoon": {"color": "spectrum"}})
v.zoomTo()
html = v._make_html()
print("    made HTML", len(html), "bytes"); sys.stdout.flush()

print("[7] reading shortlist CSV..."); sys.stdout.flush()
df = pd.read_csv(os.path.join(os.path.dirname(__file__), "target_shortlist.csv"))
print("    csv", df.shape); sys.stdout.flush()

print("\nALL STEPS PASSED — no segfault outside streamlit.")
print("If streamlit still crashes but this doesn't, the issue is streamlit's file-watcher/rerun,")
print("not the 3D code. Launch with:  streamlit run app.py --server.fileWatcherType none")
