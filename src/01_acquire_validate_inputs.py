#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1 — Acquire & validate inputs. Stream GWCD4i.DE_stats.h5ad from S3, extract the
readout z-score matrix + obs/var metadata (de_obs.parquet, de_var.parquet, zscore_f32.npy).

Conda environment: perturbseq
Produces artifact: de_obs.parquet
Data source: Zhu, Dann, ... Pritchard & Marson (2025), bioRxiv 10.64898/2025.12.23.696273.

NOTE: This file is the reproduction code captured from the artifact lineage for this
pipeline step. Steps are ordered; each consumes the checkpoints written by earlier steps.
"""

import h5py
import numpy as np
import pandas as pd
import io
import requests


class HTTPRangeFile(io.RawIOBase):
    """Minimal seekable read-only file over HTTP range requests (proxy-friendly via requests)."""
    def __init__(self, url, size=None):
        self.url = url; self._pos = 0
        self.sess = requests.Session()
        self.size = size or int(self.sess.head(url, timeout=60).headers["Content-Length"])
        self.nreads = 0; self.nbytes = 0
    def seekable(self): return True
    def readable(self): return True
    def seek(self, offset, whence=0):
        self._pos = offset if whence==0 else (self._pos+offset if whence==1 else self.size+offset)
        return self._pos
    def tell(self): return self._pos
    def read(self, n=-1):
        if n is None or n < 0: n = self.size - self._pos
        if n == 0: return b""
        end = min(self._pos + n - 1, self.size - 1)
        h = {"Range": f"bytes={self._pos}-{end}"}
        r = self.sess.get(self.url, headers=h, timeout=120); r.raise_for_status()
        data = r.content
        self._pos += len(data); self.nreads += 1; self.nbytes += len(data)
        return data
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)


url = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.DE_stats.h5ad"
rf = HTTPRangeFile(url)
h = h5py.File(rf, "r")


def cat(colname):
    g = h["obs"][colname]
    if isinstance(g, h5py.Group):  # categorical
        cats = [x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]
        codes = g["codes"][:]
        return np.array(cats)[codes]
    return g[:]


obs_gene = cat("target_contrast_gene_name")
obs_cond = cat("culture_condition")


def maybe_cat(col):
    g = h["obs"][col]
    if isinstance(g, h5py.Group):
        cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][:]]); return cats[g["codes"][:]]
    return g[:]


obs_df = pd.DataFrame({
    "row": np.arange(len(obs_gene)),
    "target_gene": obs_gene,
    "condition": obs_cond,
    "target_contrast": maybe_cat("target_contrast"),
    "n_downstream": h["obs"]["n_downstream"][:],
    "n_cells_target": h["obs"]["n_cells_target"][:],
    "ontarget_effect_size": h["obs"]["ontarget_effect_size"][:],
})
obs_df.to_parquet("de_obs.parquet")
