#!/usr/bin/env python3
"""Phase 16: evaluator replacement as another component-turnover problem.

Three evaluator versions are compared against one baseline:

- f1_refactor: exact same behavior and same evidence dependencies;
- f2_narrower_mrh: same behavior only inside a restricted MRH where witness==base,
  while consuming a smaller evidence set;
- f3_hidden_expansion: identical on the currently observed region, but introduces
  extra dependencies and a rare branch outside that region.

The point is to separate:
1. behavioral compatibility;
2. evidence/proof-horizon compatibility;
3. evaluator lineage/authorization.

This is an architecture toy, not a security proof or protocol standard.
"""

from itertools import product
from math import log

FIELDS = ("base", "witness", "mode", "override")
STATES = [dict(zip(FIELDS, vals)) for vals in product((0.0, 1.0), repeat=4)]


def f0(s):
    return 0.6 * s["base"] + 0.4 * s["witness"]


def f1_refactor(s):
    return (3.0 * s["base"] + 2.0 * s["witness"]) / 5.0


def f2_narrower_mrh(s):
    return s["base"]


def f3_hidden_expansion(s):
    out = f0(s)
    if s["mode"] > 0.5 and s["override"] > 0.5:
        out += 0.15
    return min(out, 1.0)


EVALUATORS = {
    "f0": (f0, {"base", "witness"}),
    "f1_refactor": (f1_refactor, {"base", "witness"}),
    "f2_narrower_mrh": (f2_narrower_mrh, {"base"}),
    "f3_hidden_expansion": (
        f3_hidden_expansion,
        {"base", "witness", "mode", "override"},
    ),
}

# Two deliberately restricted domains.
OBSERVED = [s for s in STATES if s["mode"] == 0.0]
CORRELATED_MRH = [s for s in OBSERVED if s["witness"] == s["base"]]
ALL = STATES


def compare(fa, fb, domain):
    diffs = [abs(fa(s) - fb(s)) for s in domain]
    return max(diffs), sum(diffs) / len(diffs), sum(d > 1e-12 for d in diffs)


def classify(old_name, new_name, domain, authorized=True, eps=0.0):
    f_old, d_old = EVALUATORS[old_name]
    f_new, d_new = EVALUATORS[new_name]
    max_diff, _, _ = compare(f_old, f_new, domain)
    behavioral = max_diff <= eps + 1e-12
    evidence_compatible = d_new <= d_old

    if not authorized:
        return "clone/unrecognized evaluator lineage"
    if behavioral and evidence_compatible:
        if d_new < d_old:
            return "compatible + proof horizon contracts"
        return "drop-in compatible in this MRH"
    if behavioral and not evidence_compatible:
        return "behaviorally compatible but proof horizon expands"
    return "semantic change: re-evaluate reliance"


def trials_for_detection(p, confidence):
    return int(log(1.0 - confidence) / log(1.0 - p)) + 1


def main():
    print("Behavioral comparisons against f0")
    print("version              domain           max_diff  changed_states  dependencies")
    for name in ("f1_refactor", "f2_narrower_mrh", "f3_hidden_expansion"):
        f, deps = EVALUATORS[name]
        for domain_name, domain in (
            ("correlated-MRH", CORRELATED_MRH),
            ("observed-MRH", OBSERVED),
            ("all-states", ALL),
        ):
            max_diff, _, changed = compare(f0, f, domain)
            print(
                f"{name:20s} {domain_name:15s} {max_diff:8.3f} "
                f"{changed:14d}  {','.join(sorted(deps))}"
            )
        print()

    print("Replacement classification")
    for domain_name, domain in (
        ("correlated-MRH", CORRELATED_MRH),
        ("observed-MRH", OBSERVED),
        ("all-states", ALL),
    ):
        for name in ("f1_refactor", "f2_narrower_mrh", "f3_hidden_expansion"):
            print(f"{domain_name:15s} {name:20s} -> {classify('f0', name, domain)}")
        print()

    # f3 diverges when mode=1, override=1, and f0 is not already capped at 1.
    # With independent uniform binary base/witness/override, that is 3/8 of mode=1 states.
    print("Rare divergence sampling burden for f3")
    for p_mode in (0.10, 0.01, 0.001, 0.0001):
        p_diverge = 0.375 * p_mode
        n95 = trials_for_detection(p_diverge, 0.95)
        n99 = trials_for_detection(p_diverge, 0.99)
        print(
            f"p_mode={p_mode:0.4f}  p_diverge={p_diverge:0.6f}  "
            f"n95={n95:7d}  n99={n99:7d}"
        )


if __name__ == "__main__":
    main()
