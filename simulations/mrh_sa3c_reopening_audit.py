#!/usr/bin/env python3
"""SA-3C: registered periodic acquisition audits, reusing frozen SA-3B models."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import mrh_sa3b_active_horizon as base


SEEDS = tuple(range(24, 48))
OPTIMAL = dict(zip(base.WORLDS, ("history", "sensor", "stop", "both")))


def audit_override(recommendation, round_number, enabled):
    forced = enabled and round_number % 8 == 0 and recommendation == "stop"
    return ("both" if forced else recommendation), forced


class AuditedObserver(base.Observer):
    def __init__(self, model, enabled):
        super().__init__(model, depth=2)
        self.enabled = enabled
        self.round_number = 0
        self.audit = False
        self.recommendation = None

    def choose(self, rounds_left):
        self.recommendation = super().choose(rounds_left)
        self.round_number += 1
        self.action, self.audit = audit_override(
            self.recommendation, self.round_number, self.enabled)
        if base.CREDITS[base.ACTIONS.index(self.action)] > self.budget:
            raise ValueError("audit exceeds acquisition budget")
        return self.action


def run_seed(model, enabled, rows):
    observer = AuditedObserver(model, enabled)
    errors, prices, logs, actions, audits, recommendations = [], [], [], [], [], []
    for t, row in enumerate(rows):
        action = observer.choose(len(rows) - t)
        observation = base.acquire(row, action, observer.budget)
        probability, prediction = observer.predict(observation)
        target = row[2]
        errors.append(float(prediction != target))
        prices.append(base.COSTS[base.ACTIONS.index(action)])
        logs.append(-math.log(probability if target else 1 - probability))
        actions.append(action)
        audits.append(observer.audit)
        recommendations.append(observer.recommendation)
        observer.feedback(target)
    mean = statistics.fmean
    unforced, _ = model.choose(observer.weights, depth=2)
    return {"utility": mean(e + c for e, c in zip(errors, prices)),
            "error": mean(errors), "price": mean(prices), "log_loss": mean(logs),
            "tail_utility": mean(e + c for e, c in zip(errors[-16:], prices[-16:])),
            "action_counts": {a: actions.count(a) for a in base.ACTIONS},
            "tail_action_counts": {a: actions[-16:].count(a) for a in base.ACTIONS},
            "final_posterior": dict(zip(model.family, observer.weights)),
            "audits": sum(audits), "tail_audits": sum(audits[-16:]),
            "audit_rounds": [t + 1 for t, forced in enumerate(audits) if forced],
            "tail_base_stop_count": recommendations[-16:].count("stop"),
            "final_unforced_recommendation": unforced}


def controls():
    inherited = base.controls()
    checks = 0

    def check(condition, label):
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    for t in range(1, 65):
        for action in base.ACTIONS:
            proposed, forced = audit_override(action, t, True)
            expected = t % 8 == 0 and action == "stop"
            check(forced == expected, "exact audit eligibility")
            check(proposed == ("both" if expected else action), "audit action")
            check(audit_override(action, t, False) == (action, False), "disabled audit")
    model = base.Model(base.WORLDS)
    known_noise = AuditedObserver(model, True)
    known_noise.weights = (0.0, 0.0, 1.0, 0.0)
    spend, audits = 0.0, 0
    for t, row in enumerate(base.sample_world("noise", 100)):
        action = known_noise.choose(64 - t)
        probability, prediction = known_noise.predict(base.acquire(row, action))
        check(abs(probability - 0.5) < base.TOL and prediction == 0,
              "known noise audit has no prediction benefit")
        known_noise.feedback(row[2])
        spend += base.COSTS[base.ACTIONS.index(action)]
        audits += known_noise.audit
    check(audits == 8, "known noise audit count")
    check(abs(spend / 64 - 0.015) < base.TOL, "known noise wasted cost")
    left, right = AuditedObserver(model, True), AuditedObserver(model, True)
    for t, row in enumerate(base.sample_world("parity", 100)):
        action = left.choose(64 - t)
        check(action == right.choose(64 - t), "public replay action")
        obs = base.acquire(row, action)
        check(left.predict(obs) == right.predict(obs), "public replay prediction")
        left.feedback(row[2])
        right.feedback(row[2])
        check(left.weights == right.weights, "public replay posterior")
    for family in (base.WORLDS[:3], base.WORLDS):
        model = base.Model(family)
        rows = base.sample_world("sensor", 0)
        previous = base.run_seed(model, 2, None, rows)
        current = run_seed(model, False, rows)
        check(all(current[k] == value for k, value in previous.items()),
              "unaudited regression against SA-3B on old seed")
    return {"inherited_sa3b_checks": inherited["checks_passed"],
            "sa3c_checks": checks}


def experiment():
    families = {"restricted": base.Model(base.WORLDS[:3]),
                "expanded": base.Model(base.WORLDS)}
    results = {}
    for world in base.WORLDS:
        streams = [base.sample_world(world, seed) for seed in SEEDS]
        data = {"oracle_expected_utility": base.ORACLE[world], "policies": {},
                "paired_comparisons": {}}
        for family, model in families.items():
            for enabled in (False, True):
                name = family + ("_audit" if enabled else "_baseline")
                records = [{"seed": seed, **run_seed(model, enabled, rows)}
                           for seed, rows in zip(SEEDS, streams)]
                summary = base.summarize(records)
                for key in ("audits", "tail_audits", "tail_base_stop_count"):
                    summary["mean_" + key] = statistics.fmean(r[key] for r in records)
                summary["final_unforced_action_counts"] = {
                    a: sum(r["final_unforced_recommendation"] == a for r in records)
                    for a in base.ACTIONS}
                data["policies"][name] = {"summary": summary, "seeds": records}
            a, b = (data["policies"][family + suffix]["seeds"]
                    for suffix in ("_audit", "_baseline"))
            diffs = [x["utility"] - y["utility"] for x, y in zip(a, b)]
            correct = OPTIMAL[world]
            comparison = {"mean": statistics.fmean(diffs),
                          "seed_standard_error": statistics.stdev(diffs) / math.sqrt(len(diffs)),
                          "paired_seed_differences": diffs,
                          "more_costly_utility_seeds": sum(d > base.TOL for d in diffs),
                          "less_costly_utility_seeds": sum(d < -base.TOL for d in diffs)}
            comparison["autonomous_rescue_seeds"] = [
                x["seed"] for x, y in zip(a, b)
                if correct != "stop" and y["final_unforced_recommendation"] == "stop"
                and x["final_unforced_recommendation"] == correct]
            comparison["lost_correct_recommendation_seeds"] = [
                x["seed"] for x, y in zip(a, b)
                if y["final_unforced_recommendation"] == correct
                and x["final_unforced_recommendation"] != correct]
            data["paired_comparisons"][family + "_audit_minus_baseline"] = comparison
        results[world] = data
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--output", type=Path, help="Create new result file; refuse overwrite")
    args = parser.parse_args()
    result = {"protocol": "SA-3C", "registration_commit": "77c79f3e",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "sa3b_source_sha256": hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
              "controls": controls()}
    if not args.controls_only:
        result["configuration"] = {"rounds": 64, "seeds": SEEDS,
                                   "audit_rounds": list(range(8, 65, 8)),
                                   "per_round_budget": 3, "base_protocol": "SA-3B"}
        result["results"] = experiment()
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
