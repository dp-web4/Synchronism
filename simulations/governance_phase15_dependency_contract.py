#!/usr/bin/env python3
"""Governance Phase 15: dependency declarations as executable contracts.

A toy evaluator has a rare branch that reads a hidden field. We compare:

- documentation-only declared dependencies;
- single-execution runtime access tracing;
- random mutation/coverage testing;
- an access-enforced evidence view where undeclared reads fail closed.

The point is not to propose Python dict wrappers as a protocol. It is to illustrate a
load-bearing architecture rule: a dependency declaration should constrain evaluator
access, not merely describe it after the fact.

Formalism/governance toy only; no security proof or physics result.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


class UndeclaredDependency(RuntimeError):
    pass


@dataclass
class EvidenceView:
    data: dict
    allowed: set[str] | None = None

    def __post_init__(self):
        self.reads: set[str] = set()

    def get(self, key: str, default=None):
        if self.allowed is not None and key not in self.allowed:
            raise UndeclaredDependency(f"evaluator attempted undeclared read: {key}")
        self.reads.add(key)
        return self.data.get(key, default)


def evaluator(view: EvidenceView) -> float:
    """Toy deterministic evaluator with one rare control-flow dependency."""
    score = float(view.get("base_score", 0.0))
    score += 0.1 * float(view.get("witness_quality", 0.0))

    # The discriminator is common and declared. The branch payload is deliberately
    # omitted from the bad declaration below to model a missed indirect/rare dependency.
    if view.get("mode") == "rare":
        score += float(view.get("rare_override", 0.0))

    return score


BAD_DECLARATION = {"base_score", "witness_quality", "mode"}
GOOD_DECLARATION = {"base_score", "witness_quality", "mode", "rare_override"}


def runtime_trace(mode: str):
    view = EvidenceView(
        {
            "base_score": 0.5,
            "witness_quality": 0.8,
            "mode": mode,
            "rare_override": 0.3,
        }
    )
    result = evaluator(view)
    return result, view.reads


def enforced_run(mode: str, declaration: set[str]):
    view = EvidenceView(
        {
            "base_score": 0.5,
            "witness_quality": 0.8,
            "mode": mode,
            "rare_override": 0.3,
        },
        allowed=declaration,
    )
    return evaluator(view), view.reads


def detection_probability(activation_probability: float, trials: int) -> float:
    """Chance random testing exercises a rare dependency at least once."""
    p = activation_probability
    return 1.0 - (1.0 - p) ** trials


def trials_for_detection(activation_probability: float, target: float) -> int:
    p = activation_probability
    return math.ceil(math.log(1.0 - target) / math.log(1.0 - p))


def main():
    print("Runtime access trace")
    for mode in ("normal", "rare"):
        result, reads = runtime_trace(mode)
        print(f"mode={mode:6s} result={result:.3f} reads={sorted(reads)}")

    print("\nBad documentation-only declaration")
    print("declared:", sorted(BAD_DECLARATION))
    print("actual rare-path dependency includes: rare_override")

    print("\nAccess-enforced evaluator")
    result, reads = enforced_run("normal", BAD_DECLARATION)
    print(f"normal path: PASS result={result:.3f} reads={sorted(reads)}")
    try:
        enforced_run("rare", BAD_DECLARATION)
        print("rare path: UNEXPECTED PASS")
    except UndeclaredDependency as exc:
        print(f"rare path: FAIL-CLOSED ({exc})")

    result, reads = enforced_run("rare", GOOD_DECLARATION)
    print(f"corrected declaration rare path: PASS result={result:.3f}")

    print("\nRandom-test chance of detecting a hidden rare dependency")
    print("activation_p   trials=10   trials=100   trials=1000   trials=10000")
    for p in (0.01, 0.001, 0.0001):
        probs = [detection_probability(p, n) for n in (10, 100, 1000, 10000)]
        print(
            f"{p:12.4f}   "
            + "   ".join(f"{v:10.6f}" for v in probs)
        )

    print("\nTrials needed to exercise hidden dependency at least once")
    print("activation_p   95% detection   99% detection")
    for p in (0.1, 0.01, 0.001, 0.0001):
        print(
            f"{p:12.4f}   {trials_for_detection(p, 0.95):13d}   "
            f"{trials_for_detection(p, 0.99):13d}"
        )


if __name__ == "__main__":
    main()
