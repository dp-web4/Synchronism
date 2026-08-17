#!/usr/bin/env python3
"""Governance Phase 14: evidence composition and semantic freshness.

A proof bundle declares which event classes can invalidate its claim. The toy compares:

1. global invalidation: ANY subsequent event makes every proof stale;
2. dependency-aware invalidation: a proof remains current until an event touching one
   of its declared dependency classes occurs.

It also composes continuity, trust, policy, and authority proofs. The action manifest
inherits the UNION of their dependency classes, so composition naturally shortens the
expected validity horizon.

This is an architecture/formalism toy. It is not a protocol proposal, security proof,
or claim of new event-sourcing/provenance mathematics.
"""

from __future__ import annotations

import numpy as np

EVENT_CLASSES = ("noise", "trust", "policy", "authority", "continuity")
EVENT_PROB = {
    "noise": 0.90,
    "trust": 0.04,
    "policy": 0.03,
    "authority": 0.02,
    "continuity": 0.01,
}

PROOF_DEPS = {
    "continuity": {"continuity"},
    "trust": {"trust", "continuity"},
    "policy": {"policy", "continuity"},
    "authority": {"authority", "policy", "continuity"},
}
PROOF_DEPS["action"] = set().union(*PROOF_DEPS.values())


def hazard(dependencies: set[str]) -> float:
    return sum(EVENT_PROB[event] for event in dependencies)


def sample_lifetimes(n=200_000, seed=20260817):
    rng = np.random.default_rng(seed)
    rows = []
    for name in ("continuity", "trust", "policy", "authority", "action"):
        h = hazard(PROOF_DEPS[name])
        # Number of events through and including the first invalidating event.
        samples = rng.geometric(h, size=n)
        rows.append(
            (
                name,
                h,
                1.0 / h,
                float(samples.mean()),
                float(np.median(samples)),
                float(np.quantile(samples, 0.95)),
            )
        )
    return rows


def is_current(created_at: int, now: int, dependencies: set[str], events: list[str]):
    for event in events[created_at:now]:
        if event in dependencies:
            return False
    return True


def frankenproof_demo():
    """A short trace where wall-clock recency gives the wrong answer.

    At t=0, all four proofs are generated.
    Events t=1..6 then occur. All proofs are only six events old, but different event
    classes invalidate different claims.
    """
    events = [
        "noise",       # t=1
        "noise",       # t=2
        "policy",      # t=3
        "noise",       # t=4
        "trust",       # t=5
        "authority",   # t=6
    ]
    now = len(events)
    out = {}
    for name in ("continuity", "trust", "policy", "authority"):
        out[name] = is_current(0, now, PROOF_DEPS[name], events)
    out["action"] = all(out.values())
    return events, out


def unrelated_noise_sweep():
    """Show that adding only irrelevant events does not age a semantic proof."""
    rows = []
    for noise_count in (0, 10, 100, 1_000, 100_000):
        events = ["noise"] * noise_count
        rows.append(
            (
                noise_count,
                is_current(0, len(events), PROOF_DEPS["continuity"], events),
                is_current(0, len(events), PROOF_DEPS["trust"], events),
                is_current(0, len(events), PROOF_DEPS["action"], events),
            )
        )
    return rows


def main():
    print("Declared event distribution")
    for event in EVENT_CLASSES:
        print(f"{event:12s} p={EVENT_PROB[event]:.3f}")

    print("\nProof validity horizons")
    print(
        "proof         invalidation_hazard  exact_mean  simulated_mean  median  p95"
    )
    for name, h, exact_mean, sim_mean, median, p95 in sample_lifetimes():
        print(
            f"{name:12s} {h:19.3f}  {exact_mean:10.3f}  "
            f"{sim_mean:14.3f}  {median:6.1f}  {p95:5.1f}"
        )

    print("\nGlobal-state invalidation baseline")
    print("If every event invalidates every proof, expected lifetime = 1.000 event.")
    print(
        "Dependency-aware action manifest expected lifetime = "
        f"{1.0 / hazard(PROOF_DEPS['action']):.3f} events."
    )

    events, current = frankenproof_demo()
    print("\nShort mixed-event trace")
    print("events:", ", ".join(events))
    for name in ("continuity", "trust", "policy", "authority", "action"):
        print(f"{name:12s} current={current[name]}")
    print(
        "All component proofs are six events old; recency alone cannot tell which "
        "claims remain current."
    )

    print("\nUnrelated-noise sweep")
    print("noise_events  continuity_current  trust_current  action_current")
    for n, c, t, a in unrelated_noise_sweep():
        print(f"{n:12d}  {str(c):18s}  {str(t):13s}  {str(a):14s}")


if __name__ == "__main__":
    main()
