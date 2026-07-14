"""
Worked example: directional nomination on your own screen.
==========================================================

Run:  python examples/worked_example.py

This builds a small synthetic screen with a planted driver perturbation,
scores it directionally, and evaluates the directional calls against the
packaged benchmark. Replace the synthetic `zmatrix` and `signature` with your
own screen to reuse the framework.
"""
import numpy as np
import pandas as pd
import directional_screen as ds


def make_demo_screen(n_pert=300, n_gene=200, seed=1):
    """A synthetic perturbation x readout-gene z-score matrix with one planted driver."""
    rng = np.random.default_rng(seed)
    genes = [f"READ{i}" for i in range(n_gene)]
    perts = [f"PERT{i}" for i in range(n_pert)]
    Z = pd.DataFrame(rng.standard_normal((n_pert, n_gene)), index=perts, columns=genes)
    # signed program: first 12 genes pro-program (+1), next 8 anti-program (-1)
    signature = pd.Series({**{g: 1 for g in genes[:12]}, **{g: -1 for g in genes[12:20]}})
    # plant PERT0 as a driver: knockdown lowers the program
    Z.loc["PERT0", genes[:12]] -= 3
    Z.loc["PERT0", genes[12:20]] += 3
    return Z, signature


def main():
    # 1. Load or build your screen: zmatrix (perturbations x readout genes) + signed signature
    zmatrix, signature = make_demo_screen()

    # 2. Directional scoring (signed program score + empirical-null FDR + direction call)
    scored = ds.score_screen(zmatrix, signature, n_null=200, seed=0)
    print("Top driver (antagonize) and brake (agonize) hits:")
    print(scored.sort_values("program_score").head(3)[["program_score", "emp_z", "emp_fdr", "direction"]])
    print(scored.sort_values("program_score").tail(3)[["program_score", "emp_z", "emp_fdr", "direction"]])

    # 3. Score your directional calls against the packaged benchmark
    #    (map your gene symbols onto the benchmark; here we just show the API)
    bench = ds.load_benchmark()
    print(f"\nPackaged benchmark: {len(bench)} directional gold targets "
          f"({(bench.true_direction=='driver_antagonize').sum()} driver / "
          f"{(bench.true_direction=='brake_agonize').sum()} brake).")
    result = ds.evaluate_directions(bench.set_index("gene")["framework_call"])
    print("Directional evaluation (honest metrics):")
    for k in ["n_evaluated", "raw_accuracy", "balanced_accuracy", "permutation_p", "majority_baseline"]:
        print(f"  {k}: {result[k]}")
    print("  ->", result["note"])


if __name__ == "__main__":
    main()
