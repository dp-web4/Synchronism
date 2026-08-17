#!/usr/bin/env python3
"""Phase-1 toy model for the Markov coherence/emergence exploration arc.

Purpose
-------
Test a deliberately minimal candidate for dynamic coherence:

    coherence_ratio = tau_exit / tau_internal

where tau_exit is the metastable residence-time proxy and tau_internal is an
internal relaxation-time proxy.  The model also reports a simple strong-
lumpability error for a two-block partition.

This is a formalism probe, not a physics simulation and not evidence for
Synchronism.
"""

from __future__ import annotations

import math
import numpy as np


DYNAMIC_WITHIN = np.array(
    [
        [0.62, 0.23, 0.15],
        [0.18, 0.64, 0.18],
        [0.20, 0.20, 0.60],
    ],
    dtype=float,
)

FROZEN_WITHIN = np.eye(3, dtype=float)


def build_chain(leakage: float, asymmetry: float = 0.0, frozen: bool = False) -> np.ndarray:
    """Build two 3-state clusters with tunable inter-cluster leakage.

    asymmetry=0 gives exact strong lumpability for partition {0,1,2}/{3,4,5}.
    Positive asymmetry makes the per-state leakage rates differ inside each
    cluster and therefore breaks exact lumpability.
    """
    if not (0.0 < leakage < 1.0):
        raise ValueError("leakage must be between 0 and 1")

    within = FROZEN_WITHIN if frozen else DYNAMIC_WITHIN
    offsets = np.array([-1.0, 0.0, 1.0])
    state_leak = leakage * (1.0 + asymmetry * offsets)
    if np.any(state_leak <= 0.0) or np.any(state_leak >= 1.0):
        raise ValueError("asymmetry produced invalid state leakage")

    p = np.zeros((6, 6), dtype=float)

    for block in (0, 1):
        here = 3 * block
        there = 3 * (1 - block)
        for i in range(3):
            e_i = state_leak[i]
            p[here + i, here : here + 3] = (1.0 - e_i) * within[i]
            p[here + i, there : there + 3] = e_i / 3.0

    if not np.allclose(p.sum(axis=1), 1.0):
        raise RuntimeError("transition matrix rows do not sum to one")
    return p


def stationary_distribution(p: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(p.T)
    idx = int(np.argmin(np.abs(vals - 1.0)))
    pi = np.real(vecs[:, idx])
    if pi.sum() < 0:
        pi = -pi
    return pi / pi.sum()


def cluster_metrics(p: np.ndarray) -> dict[str, float]:
    """Metrics for cluster {0,1,2}.

    tau_exit uses the Perron root rho of the substochastic within-cluster
    matrix Q: 1/(1-rho), a standard near-metastability residence-time proxy.

    tau_internal uses 1/(1-|lambda_2|) for the row-normalized conditional
    within-cluster chain.  It is a relaxation-time proxy, not a rigorous total-
    variation mixing time.
    """
    q = p[:3, :3]
    q_eigs = np.linalg.eigvals(q)
    rho = float(np.max(np.abs(q_eigs)))
    tau_exit = 1.0 / (1.0 - rho)

    conditional = q / q.sum(axis=1, keepdims=True)
    eig_abs = sorted((float(abs(v)) for v in np.linalg.eigvals(conditional)), reverse=True)
    lambda_2 = eig_abs[1]

    if lambda_2 >= 1.0 - 1e-12:
        tau_internal = math.inf
        coherence_ratio = 0.0
    else:
        tau_internal = 1.0 / (1.0 - lambda_2)
        coherence_ratio = tau_exit / tau_internal

    return {
        "rho": rho,
        "tau_exit": tau_exit,
        "lambda_2_internal": lambda_2,
        "tau_internal": tau_internal,
        "coherence_ratio": coherence_ratio,
    }


def lumpability_error(p: np.ndarray) -> float:
    """Maximum violation of strong lumpability for the two-block partition.

    For states i,j in the same source block, exact strong lumpability requires
    equal aggregate transition probability into each target block.
    """
    blocks = (range(0, 3), range(3, 6))
    worst = 0.0
    for source in blocks:
        for target in blocks:
            target_idx = list(target)
            totals = [float(p[i, target_idx].sum()) for i in source]
            worst = max(worst, max(totals) - min(totals))
    return worst


def report_arm(label: str, asymmetry: float) -> list[dict[str, float]]:
    leakages = (0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30)
    rows = []
    print(f"\n{label}")
    print("leakage  tau_exit  tau_internal  coherence_ratio  lump_error")
    for leakage in leakages:
        p = build_chain(leakage, asymmetry=asymmetry)
        m = cluster_metrics(p)
        err = lumpability_error(p)
        row = {"leakage": leakage, **m, "lumpability_error": err}
        rows.append(row)
        print(
            f"{leakage:7.3f}  {m['tau_exit']:8.3f}  {m['tau_internal']:12.3f}"
            f"  {m['coherence_ratio']:15.3f}  {err:10.6f}"
        )
    return rows


def main() -> None:
    exact = report_arm("EXACTLY LUMPABLE ARM", asymmetry=0.0)
    perturbed = report_arm("PERTURBED ARM (asymmetry=0.25)", asymmetry=0.25)

    frozen_p = build_chain(0.02, asymmetry=0.0, frozen=True)
    frozen = cluster_metrics(frozen_p)
    print("\nFROZEN CONTROL (leakage=0.02)")
    print(f"tau_exit        = {frozen['tau_exit']:.3f}")
    print(f"tau_internal    = {frozen['tau_internal']}")
    print(f"coherence_ratio = {frozen['coherence_ratio']:.3f}")

    # Sanity checks: these are the qualitative claims Phase 1 is intended to test.
    exact_scores = [r["coherence_ratio"] for r in exact]
    assert all(a > b for a, b in zip(exact_scores, exact_scores[1:])), (
        "coherence ratio should decrease as leakage rises in this construction"
    )
    assert max(r["lumpability_error"] for r in exact) < 1e-12
    assert all(r["lumpability_error"] > 0 for r in perturbed)
    assert frozen["coherence_ratio"] == 0.0

    p = build_chain(0.02)
    pi = stationary_distribution(p)
    print("\nSTATIONARY DISTRIBUTION (exact arm, leakage=0.02)")
    print(np.array2string(pi, precision=6, suppress_small=True))

    print("\nAll Phase-1 qualitative checks passed.")


if __name__ == "__main__":
    main()
