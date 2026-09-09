#!/usr/bin/env python3
"""SA-3D: registered lifetime-mixture evidence gate after first stopping."""

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import mrh_sa3b_active_horizon as base
import mrh_sa3c_reopening_audit as previous

THRESHOLD = 20.0
MAIN_SEEDS = tuple(range(48, 72))
NULL_SEEDS = tuple(range(1024, 1536))


def gated_action(recommendation, stopped, evidence, enabled):
    veto = enabled and stopped and evidence < THRESHOLD and recommendation != "stop"
    return ("stop" if veto else recommendation), veto


class EvidenceObserver(base.Observer):
    def __init__(self, model, mode):
        super().__init__(model, depth=2)
        if mode not in ("none", "periodic", "gated"):
            raise ValueError("unknown policy mode")
        self.mode = mode
        self.null_index = model.family.index("noise")
        self.alternatives = tuple(i for i in range(len(model.family)) if i != self.null_index)
        self.ratios = [1.0] * len(model.family)
        self.evidence = self.max_evidence = 1.0
        self.round_number = 0
        self.first_stop = self.first_crossing = self.first_reopening = None
        self.audit = self.veto = self.post_stop_query = False
        self.recommendation = None
        self.pre_evidence = None

    def choose(self, rounds_left):
        self.recommendation = super().choose(rounds_left)
        self.round_number += 1
        self.pre_evidence = self.evidence
        if self.recommendation == "stop" and self.first_stop is None:
            self.first_stop = self.round_number
        proposed, self.veto = gated_action(self.recommendation, self.first_stop is not None,
                                          self.evidence, self.mode == "gated")
        self.action, self.audit = previous.audit_override(
            proposed, self.round_number, self.mode != "none")
        if base.CREDITS[base.ACTIONS.index(self.action)] > self.budget:
            raise ValueError("acquisition exceeds budget")
        self.post_stop_query = (self.first_stop is not None and self.action != "stop"
                                and not self.audit)
        if self.post_stop_query:
            if self.mode == "gated" and self.pre_evidence < THRESHOLD:
                raise AssertionError("post-stop unforced query without prior evidence")
            if self.first_reopening is None:
                self.first_reopening = self.round_number
        return self.action

    def feedback(self, target):
        # Validate phase/target via the original observer before changing evidence.
        pending = self.pending
        super().feedback(target)
        vector = pending[target]
        denominator = vector[self.null_index]
        for i in self.alternatives:
            self.ratios[i] *= vector[i] / denominator
        self.evidence = statistics.fmean(self.ratios[i] for i in self.alternatives)
        self.max_evidence = max(self.max_evidence, self.evidence)
        if self.evidence >= THRESHOLD and self.first_crossing is None:
            self.first_crossing = self.round_number

    def diagnostic(self):
        proposed, _ = self.model.choose(self.weights, depth=2)
        allowed, _ = gated_action(proposed, self.first_stop is not None,
                                  self.evidence, self.mode == "gated")
        return allowed, proposed


def run_seed(model, mode, rows):
    observer = EvidenceObserver(model, mode)
    errors, prices, logs, actions, audits, recommendations, vetoes = [], [], [], [], [], [], []
    violations = 0
    for t, row in enumerate(rows):
        action = observer.choose(len(rows) - t)
        if observer.post_stop_query and mode == "gated" and observer.pre_evidence < THRESHOLD:
            violations += 1
        probability, prediction = observer.predict(base.acquire(row, action, observer.budget))
        target = row[2]
        errors.append(float(prediction != target))
        prices.append(base.COSTS[base.ACTIONS.index(action)])
        logs.append(-math.log(probability if target else 1 - probability))
        actions.append(action)
        audits.append(observer.audit)
        vetoes.append(observer.veto)
        recommendations.append(observer.recommendation)
        observer.feedback(target)
    mean = statistics.fmean
    unforced, ungated = observer.diagnostic()
    return {"utility": mean(e + c for e, c in zip(errors, prices)),
            "error": mean(errors), "price": mean(prices), "log_loss": mean(logs),
            "tail_utility": mean(e + c for e, c in zip(errors[-16:], prices[-16:])),
            "action_counts": {a: actions.count(a) for a in base.ACTIONS},
            "tail_action_counts": {a: actions[-16:].count(a) for a in base.ACTIONS},
            "final_posterior": dict(zip(model.family, observer.weights)),
            "audits": sum(audits), "tail_audits": sum(audits[-16:]),
            "audit_rounds": [t + 1 for t, x in enumerate(audits) if x],
            "gate_veto_rounds": [t + 1 for t, x in enumerate(vetoes) if x],
            "tail_base_stop_count": recommendations[-16:].count("stop"),
            "final_unforced_recommendation": unforced,
            "final_ungated_recommendation": ungated,
            "first_stop": observer.first_stop, "first_crossing": observer.first_crossing,
            "first_reopening": observer.first_reopening,
            "ever_post_stop_query": observer.first_reopening is not None,
            "evidence_final": observer.evidence, "evidence_max": observer.max_evidence,
            "evidence_gate_violations": violations}


def rational_probability(world, action, obs):
    """Independent exact P(Y=1|observation,world), not model-table arithmetic."""
    if world == "memory" and action in ("history", "both"):
        return F(1, 10) + F(4, 5) * obs[0]
    if world == "sensor" and action in ("sensor", "both"):
        return F(9, 50) + F(16, 25) * obs[-1]
    if world == "parity" and action == "both":
        return F(9, 50) + F(16, 25) * (obs[0] ^ obs[1])
    return F(1, 2)


def controls():
    inherited = previous.controls()
    checks = 0

    def check(condition, label):
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    for family in (base.WORLDS[:3], base.WORLDS):
        model = base.Model(family)
        alternatives = [w for w in family if w != "noise"]
        for action, cells in model.events.items():
            for obs, v0, v1 in cells:
                for k, world in enumerate(family):
                    p1 = rational_probability(world, action, obs)
                    check(abs(v1[k] - float(p1 / len(cells))) < base.TOL, "rational likelihood 1")
                    check(abs(v0[k] - float((1 - p1) / len(cells))) < base.TOL, "rational likelihood 0")
            for wealth in (F(1), F(2, 7), F(19)):
                for world in alternatives:
                    expectation = F(0)
                    for obs, _, _ in cells:
                        p1 = rational_probability(world, action, obs)
                        for p in (p1, 1 - p1):
                            null_mass = F(1, 2 * len(cells))
                            ratio = p * 2
                            expectation += null_mass * wealth * ratio
                    check(expectation == wealth, "exact conditional likelihood-ratio martingale")
        for mode in ("none", "periodic"):
            rows = base.sample_world("sensor", 0)
            old = previous.run_seed(model, mode == "periodic", rows)
            new = run_seed(model, mode, rows)
            check(all(new[k] == v for k, v in old.items()), "SA-3C baseline regression")
        observer = EvidenceObserver(model, "gated")
        for t, row in enumerate(base.sample_world("noise", 100)):
            before = tuple(observer.ratios)
            before_evidence = observer.evidence
            action = observer.choose(64 - t)
            check(tuple(observer.ratios) == before, "choosing never resets evidence")
            observer.predict(base.acquire(row, action))
            check(observer.evidence == before_evidence, "current target unavailable at purchase")
            observer.feedback(row[2])
            if action == "stop":
                check(abs(observer.evidence - before_evidence) < base.TOL, "stop changes no evidence")
            odds_e = (1 - observer.weights[observer.null_index]) / (
                len(alternatives) * observer.weights[observer.null_index])
            check(abs(observer.evidence - odds_e) < 1e-10, "independent posterior-odds consistency")
        ratios = tuple(observer.ratios)
        try:
            observer.feedback(0)
        except ValueError:
            check(tuple(observer.ratios) == ratios, "double feedback rejects without evidence update")
        else:
            check(False, "double feedback must reject")
    for stopped in (False, True):
        for evidence in (0.0, 19.999, 20.0, 20.001):
            for action in base.ACTIONS:
                proposed, veto = gated_action(action, stopped, evidence, True)
                expected = stopped and evidence < 20 and action != "stop"
                check(veto == expected, "gate eligibility boundary")
                check(proposed == ("stop" if expected else action), "gate action boundary")
    action, veto = gated_action("sensor", True, 1, True)
    check(veto and previous.audit_override(action, 8, True) == ("both", True),
          "scheduled audit survives gate veto")
    return {**inherited, "sa3d_checks": checks}


def summarize(records):
    summary = base.summarize(records)
    for key in ("audits", "tail_audits"):
        summary["mean_" + key] = statistics.fmean(r[key] for r in records)
    summary["mean_gate_vetoes"] = statistics.fmean(len(r["gate_veto_rounds"]) for r in records)
    summary["ever_post_stop_query_count"] = sum(r["ever_post_stop_query"] for r in records)
    summary["evidence_crossing_count"] = sum(r["first_crossing"] is not None for r in records)
    summary["ever_stopped_count"] = sum(r["first_stop"] is not None for r in records)
    summary["final_unforced_action_counts"] = {
        a: sum(r["final_unforced_recommendation"] == a for r in records) for a in base.ACTIONS}
    summary["evidence_gate_violations"] = sum(r["evidence_gate_violations"] for r in records)
    return summary


def cohort(worlds, seeds, modes):
    results = {}
    families = {"restricted": base.Model(base.WORLDS[:3]), "expanded": base.Model(base.WORLDS)}
    for world in worlds:
        streams = [base.sample_world(world, seed) for seed in seeds]
        data = {"policies": {}, "paired_comparisons": {}}
        for family, model in families.items():
            for mode in modes:
                records = [{"seed": seed, **run_seed(model, mode, rows)}
                           for seed, rows in zip(seeds, streams)]
                data["policies"][family + "_" + mode] = {
                    "summary": summarize(records), "seeds": records}
            pairs = [("gated", "periodic")]
            if "none" in modes:
                pairs += [("periodic", "none"), ("gated", "none")]
            for first, second in pairs:
                a, b = (data["policies"][family + "_" + m]["seeds"] for m in (first, second))
                diffs = [x["utility"] - y["utility"] for x, y in zip(a, b)]
                data["paired_comparisons"][family + "_" + first + "_minus_" + second] = {
                    "mean": statistics.fmean(diffs),
                    "seed_standard_error": statistics.stdev(diffs) / math.sqrt(len(diffs)),
                    "paired_seed_differences": diffs}
        results[world] = data
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--output", type=Path, help="Create new artifact; refuse overwrite")
    args = parser.parse_args()
    result = {"protocol": "SA-3D", "registration_commit": "a323d705",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "sa3b_source_sha256": hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
              "sa3c_source_sha256": hashlib.sha256(Path(previous.__file__).read_bytes()).hexdigest(),
              "controls": controls()}
    if not args.controls_only:
        result["configuration"] = {"rounds": 64, "threshold": THRESHOLD, "alpha": 0.05,
                                   "main_seeds": MAIN_SEEDS, "null_seeds": NULL_SEEDS,
                                   "audit_rounds": list(range(8, 65, 8)), "per_round_budget": 3}
        result["main_cohort"] = cohort(base.WORLDS, MAIN_SEEDS, ("none", "periodic", "gated"))
        result["null_cohort"] = cohort(("noise",), NULL_SEEDS, ("periodic", "gated"))
    rendered = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as output:
            output.write(rendered)
        print(json.dumps({"output": str(args.output), "controls": result["controls"],
                          "source_sha256": result["source_sha256"]}, sort_keys=True))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
