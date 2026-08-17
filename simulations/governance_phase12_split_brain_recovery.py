#!/usr/bin/env python3
"""Governance Phase 12: split-brain recovery and assurance-relative continuity.

A deliberately small Monte Carlo model of recovery after a crash.

Timeline per episode:
- A is the original live token/runtime and crashes at t=0.
- B is restored from backup after a random restore delay R.
- Network/registry connectivity returns after random partition delay P.
- With probability p_resurrect, old A later comes back after delay O.
- An authoritative registry can finalize/fence the recovery only after both B exists
  and connectivity is restored, plus a small finalization delay. This occurs at F.

Two relying-party profiles:
- A1 / low-stakes: accepts B when signed descent/recovery evidence exists at R,
  even though exclusivity remains unknown.
- A2 / high-stakes: waits until F for scoped registry/fencing evidence supporting
  exclusive continuation.

The model separates:
- actual physical/runtime state (whether A and B are both live),
- evidence state (what a relying party can justify),
- acceptance policy (what evidence the relying party requires).

This is an architecture toy, not a security proof, protocol proposal, or physics model.
"""

import numpy as np


DEFAULTS = dict(
    mean_restore=8.0,
    mean_partition=30.0,
    p_resurrect=0.60,
    mean_resurrection=18.0,
    mean_finalize=5.0,
)


def simulate(
    n=200_000,
    seed=12345,
    mean_restore=DEFAULTS["mean_restore"],
    mean_partition=DEFAULTS["mean_partition"],
    p_resurrect=DEFAULTS["p_resurrect"],
    mean_resurrection=DEFAULTS["mean_resurrection"],
    mean_finalize=DEFAULTS["mean_finalize"],
):
    rng = np.random.default_rng(seed)

    restore = rng.exponential(mean_restore, n)
    partition_end = rng.exponential(mean_partition, n)

    does_resurrect = rng.random(n) < p_resurrect
    resurrection = np.where(
        does_resurrect,
        rng.exponential(mean_resurrection, n),
        np.inf,
    )

    # Registry/fencing finality cannot occur before both the recovered candidate exists
    # and registry connectivity is back. The finalization delay stands in for the
    # witnessed fence / registry commitment / quorum work after those prerequisites.
    finalized = np.maximum(restore, partition_end) + rng.exponential(mean_finalize, n)

    # If the old runtime resurrects before finalization, there is an interval during
    # which A and B are simultaneously live. B does not exist before restore.
    dual_start = np.maximum(restore, resurrection)
    dual_live = np.where(dual_start < finalized, finalized - dual_start, 0.0)
    has_dual_live = dual_live > 0

    # A1 accepts at restore; A2 waits for finalization/fencing evidence.
    assurance_gap = finalized - restore

    return {
        "n": n,
        "mean_restore": float(restore.mean()),
        "mean_finalized": float(finalized.mean()),
        "mean_assurance_gap": float(assurance_gap.mean()),
        "median_assurance_gap": float(np.median(assurance_gap)),
        "p_dual_live_before_finality": float(has_dual_live.mean()),
        "mean_dual_live_all": float(dual_live.mean()),
        "mean_dual_live_conditional": float(
            dual_live[has_dual_live].mean() if has_dual_live.any() else 0.0
        ),
        "p_resurrect_before_finality": float((resurrection < finalized).mean()),
        "p_resurrect_after_finality": float(
            (does_resurrect & (resurrection >= finalized)).mean()
        ),
        "p_no_resurrection": float((~does_resurrect).mean()),
        "p_partition_dominates_restore": float((partition_end > restore).mean()),
    }


def sweep_partition_delay():
    print("\nPartition-delay sweep")
    print(
        "mean_partition  mean_A2_minus_A1  p_dual_live  "
        "mean_dual_if_present"
    )
    for mean_partition in (2, 5, 10, 20, 30, 60, 120):
        out = simulate(
            n=100_000,
            seed=1000 + int(mean_partition * 10),
            mean_partition=float(mean_partition),
        )
        print(
            f"{mean_partition:14.1f}  "
            f"{out['mean_assurance_gap']:16.3f}  "
            f"{out['p_dual_live_before_finality']:11.4f}  "
            f"{out['mean_dual_live_conditional']:20.3f}"
        )


def main():
    out = simulate()

    print("Baseline split-brain recovery model")
    print(f"episodes: {out['n']}")
    print(f"mean restore time / A1 acceptance: {out['mean_restore']:.3f}")
    print(f"mean fence-finality / A2 acceptance: {out['mean_finalized']:.3f}")
    print(f"mean A2-A1 assurance latency: {out['mean_assurance_gap']:.3f}")
    print(f"median A2-A1 assurance latency: {out['median_assurance_gap']:.3f}")
    print(
        "probability of actual dual-live interval before finality: "
        f"{out['p_dual_live_before_finality']:.4f}"
    )
    print(
        "mean dual-live duration conditional on occurrence: "
        f"{out['mean_dual_live_conditional']:.3f}"
    )
    print(
        "probability old runtime resurrects only after finality: "
        f"{out['p_resurrect_after_finality']:.4f}"
    )
    print(
        "probability old runtime never resurrects: "
        f"{out['p_no_resurrection']:.4f}"
    )
    print(
        "probability registry partition outlasts restore: "
        f"{out['p_partition_dominates_restore']:.4f}"
    )

    sweep_partition_delay()


if __name__ == "__main__":
    main()
