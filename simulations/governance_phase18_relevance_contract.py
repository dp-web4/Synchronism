#!/usr/bin/env python3
"""Phase 18: a generic relevance-contract composition toy.

A contract describes the basis of one derived claim:
- subject;
- evaluator/plan;
- authoritative state position;
- dependency classes;
- completeness-for-claim;
- assurance basis;
- prior-claim lineage.

Composition carries each claim forward to one current state only when no intervening
event intersects its dependencies. It preserves assurance as a vector rather than
collapsing it to a scalar. Incomplete projections fail closed.

This is an architecture sketch, not a protocol or security proof.
"""

from dataclasses import dataclass
from typing import FrozenSet


@dataclass(frozen=True)
class RelevanceContract:
    name: str
    subject: str
    plan: str
    basis_position: int
    dependencies: FrozenSet[str]
    complete_for_claim: bool
    assurance: str
    prior_claim: str | None = None


def dependency_touches(history, start_position, end_position, dependencies):
    """Return relevant intervening events using 1-based displayed positions."""
    return [
        (i, event_type)
        for i, event_type in enumerate(
            history[start_position:end_position], start=start_position + 1
        )
        if event_type in dependencies
    ]


def carryable(contract, history, current_position):
    return contract.complete_for_claim and not dependency_touches(
        history,
        contract.basis_position,
        current_position,
        contract.dependencies,
    )


def compose(contracts, history, current_position):
    reasons = []
    subjects = {c.subject for c in contracts}
    if len(subjects) != 1:
        reasons.append("subject mismatch")

    for contract in contracts:
        if not contract.complete_for_claim:
            reasons.append(f"{contract.name}: incomplete-for-claim")
        touches = dependency_touches(
            history,
            contract.basis_position,
            current_position,
            contract.dependencies,
        )
        if touches:
            pos, event_type = touches[0]
            reasons.append(
                f"{contract.name}: stale via {event_type}@{pos}"
            )

    dependencies = frozenset().union(
        *(c.dependencies for c in contracts)
    )

    return {
        "ok": not reasons,
        "reasons": reasons,
        "composed_dependencies": sorted(dependencies),
        "basis_positions": [c.basis_position for c in contracts],
        "plans": [c.plan for c in contracts],
        # Deliberately keep assurance typed/per-component rather than averaging it.
        "assurance_vector": [c.assurance for c in contracts],
    }


def renew(contract, current_position):
    """Issue a distinct renewed claim; do not relabel the old claim as current."""
    return RelevanceContract(
        name=contract.name,
        subject=contract.subject,
        plan=contract.plan,
        basis_position=current_position,
        dependencies=contract.dependencies,
        complete_for_claim=True,
        assurance=contract.assurance,
        prior_claim=f"renewed-from:{contract.prior_claim or contract.name}",
    )


def main():
    history = ["noise"] * 20
    history[3] = "policy"       # displayed position 4
    history[7] = "trust"        # position 8
    history[14] = "identity"    # position 15
    history[17] = "authority"   # position 18
    current = 20

    contracts = [
        # Issued after the identity event: later authority churn is irrelevant.
        RelevanceContract(
            "continuity", "society-A", "continuity-v2", 16,
            frozenset({"identity"}), True, "A2", "C-cont-1"
        ),
        # Older policy/trust proofs are invalidated by the later identity change.
        RelevanceContract(
            "policy", "society-A", "policy-v4", 5,
            frozenset({"policy", "identity"}), True, "A2", "C-pol-7"
        ),
        RelevanceContract(
            "trust", "society-A", "trust-v8", 9,
            frozenset({"trust", "identity"}), True, "A1", "C-trust-9"
        ),
        RelevanceContract(
            "authority", "society-A", "authority-v3", 18,
            frozenset({"authority", "policy", "identity"}), True, "A2", "C-auth-2"
        ),
    ]

    print(f"current position: {current}")
    for c in contracts:
        print(
            f"{c.name:10s} carryable={carryable(c, history, current)} "
            f"touches={dependency_touches(history, c.basis_position, current, c.dependencies)}"
        )

    print("\ncomposition before renewal:")
    print(compose(contracts, history, current))

    renewed = [
        renew(c, current) if not carryable(c, history, current) else c
        for c in contracts
    ]

    print("\ncomposition after renewing stale claims:")
    print(compose(renewed, history, current))

    cheap_incomplete = RelevanceContract(
        "trust-cheap",
        "society-A",
        "trust-v8",
        current,
        frozenset({"trust"}),
        False,
        "A1",
    )

    print("\nincomplete projected reading composed with continuity:")
    print(compose([renewed[0], cheap_incomplete], history, current))


if __name__ == "__main__":
    main()
