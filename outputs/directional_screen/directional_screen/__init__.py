"""
directional_screen
===================

A small, assay-agnostic toolkit for **directional** target nomination from
massively-multiplexed perturbation screens.

Where a standard screen analysis ranks perturbations by *effect magnitude*, this
toolkit assigns each perturbation a signed **direction** — driver (antagonize)
vs brake (agonize) — from the sign of its effect on a signed gene program, with
an empirical-null significance test and external-truth direction calibration.

Quickstart
----------
    import pandas as pd
    from directional_screen import score_screen, load_benchmark, evaluate_directions

    # zmatrix: perturbations x readout-genes (z-scores); signature: gene -> +/-1
    scored = score_screen(zmatrix, signature)          # signed score + emp FDR + direction
    calls  = scored["direction"]
    print(evaluate_directions(calls))                  # score vs packaged benchmark

Modules
-------
    scoring      : program_score, empirical_null, score_screen
    calibration  : concordance_labels, fit_direction_confidence, reliability, brier_skill
    benchmark    : load_benchmark, load_dictionary, evaluate_directions
"""
from .scoring import program_score, empirical_null, bh_fdr, score_screen
from .calibration import (
    concordance_labels,
    fit_direction_confidence,
    reliability,
    brier_skill,
)
from .benchmark import load_benchmark, load_dictionary, evaluate_directions

__version__ = "1.0.0"
__all__ = [
    "program_score", "empirical_null", "bh_fdr", "score_screen",
    "concordance_labels", "fit_direction_confidence", "reliability", "brier_skill",
    "load_benchmark", "load_dictionary", "evaluate_directions",
]
