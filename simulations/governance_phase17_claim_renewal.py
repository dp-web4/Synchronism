#!/usr/bin/env python3
"""Phase 17: claim renewal across entity/evaluator evolution.

Model three separate histories:
- entity lineage: which governed token exists now;
- evaluator lineage: which authorized plan computes the claim;
- claim lineage: each derived result is a new historical claim object.

The toy asks when a new claim can be derived incrementally from an old complete-for-
claim proof plus relevant delta evidence, and when a new evaluator dependency forces
historical backfill.

Dependencies are tagged as either:
- path-sensitive: relevant events since the relying party's anchor matter;
- current-state: only a current authenticated snapshot/commitment is needed.

This is an architecture toy, not a protocol standard or complexity proof.
"""

import random
from collections import Counter

SEED = 170817
N_EVENTS = 100_000
ANCHOR = 10_000
OLD_BASIS = 90_000

RATES = (
    ("identity", 0.001),
    ("policy", 0.005),
    ("trust", 0.010),
    ("sanction", 0.002),
)


def make_history(n=N_EVENTS, seed=SEED):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        u = rng.random()
        acc = 0.0
        event_type = "noise"
        for name, probability in RATES:
            acc += probability
            if u < acc:
                event_type = name
                break
        out.append(event_type)
    return out


def count(history, lo, hi, event_types):
    return sum(1 for e in history[lo:hi] if e in event_types)


def main():
    history = make_history()
    old_path_dependencies = {"identity", "policy"}

    old_proof = count(history, ANCHOR, OLD_BASIS, old_path_dependencies)
    delta_same = count(history, OLD_BASIS, N_EVENTS, old_path_dependencies)
    contracted_delta = count(history, OLD_BASIS, N_EVENTS, {"identity"})

    # A newly required current-state dependency can be satisfied from current state.
    # 21 is a schematic authenticated-snapshot/commitment proof cost.
    current_snapshot_cost = 21
    expanded_current_new_fetch = delta_same + current_snapshot_cost

    # A newly required path-sensitive dependency was NOT preserved by the old proof.
    # It must be backfilled from the admission anchor to the old basis, then continued.
    sanction_backfill = count(history, ANCHOR, OLD_BASIS, {"sanction"})
    sanction_delta = count(history, OLD_BASIS, N_EVENTS, {"sanction"})
    expanded_path_new_fetch = delta_same + sanction_backfill + sanction_delta

    full_rebuild_same = count(
        history, ANCHOR, N_EVENTS, old_path_dependencies
    )
    full_rebuild_expanded_path = count(
        history,
        ANCHOR,
        N_EVENTS,
        old_path_dependencies | {"sanction"},
    )

    print("whole history event counts:", Counter(history))
    print(f"old complete-for-claim proof objects: {old_proof}")
    print()
    print("new material needed to renew at current head:")
    print(f"same evaluator dependencies:        {delta_same}")
    print(f"contracted dependency evaluator:    {contracted_delta}")
    print(
        f"expanded current-state dependency: {expanded_current_new_fetch} "
        f"({current_snapshot_cost} snapshot-proof objects)"
    )
    print(
        f"expanded path dependency:          {expanded_path_new_fetch} "
        f"({sanction_backfill} historical backfill + {sanction_delta} delta)"
    )
    print()
    print("full rebuild comparison:")
    print(f"same dependencies full rebuild:     {full_rebuild_same}")
    print(f"expanded-path full rebuild:         {full_rebuild_expanded_path}")
    print()
    print("new-fetch fraction of corresponding rebuild:")
    print(f"same:             {delta_same / full_rebuild_same:.4f}")
    print(f"contracted:       {contracted_delta / full_rebuild_same:.4f}")
    print(
        "expanded-current: "
        f"{expanded_current_new_fetch / (full_rebuild_same + current_snapshot_cost):.4f}"
    )
    print(f"expanded-path:    {expanded_path_new_fetch / full_rebuild_expanded_path:.4f}")

    print("\nExpected scaling with fixed 10,000-event delta:")
    p_old = 0.001 + 0.005
    p_new_path = 0.002
    delta = 10_000
    print("old_history_span  same-plan-new  new-path-backfill+delta")
    for historical_span in (10_000, 100_000, 1_000_000, 10_000_000):
        same = delta * p_old
        expanded_path = delta * p_old + (historical_span + delta) * p_new_path
        print(f"{historical_span:16d}  {same:13.1f}  {expanded_path:23.1f}")


if __name__ == "__main__":
    main()
