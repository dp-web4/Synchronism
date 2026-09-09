#!/usr/bin/env python3
"""SA-3F: bounded bursts with exact reachability and named-model power reports."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import mrh_sa3d_evidence_gate as gate

base = gate.base
THRESHOLD = F(20)
POWER_TARGET = F(4, 5)
POLICIES = {"periodic24": 24, "burst8": 8, "burst24": 24, "feasible24": 24}
MAIN_SEEDS = tuple(range(72, 96))
NULL_SEEDS = tuple(range(2048, 2176))


def agreement_probability(world):
    return F(9, 10) if world == "memory" else F(41, 50)


def factor(world, action, observation, target):
    p = gate.rational_probability(world, action, observation)
    return 2 * (p if target else 1 - p)


def envelope(ratios, audits):
    return sum((value * (2 * agreement_probability(world)) ** audits
                for world, value in ratios.items()), F(0)) / len(ratios)


@lru_cache(maxsize=4096)
def component_power(start, probability, audits, components):
    """Exact P(component LR ever >= K*20); sufficient for mixture crossing."""
    boundary = components * THRESHOLD
    if start >= boundary:
        return F(1)
    up, down = 2 * probability, 2 * (1 - probability)
    if start * up ** audits < boundary:
        return F(0)
    live, crossed = {0: F(1)}, F(0)
    for n in range(1, audits + 1):
        next_live = {}
        for successes, mass in live.items():
            for outcome, chance in ((1, probability), (0, 1 - probability)):
                s = successes + outcome
                amount = mass * chance
                if start * up ** s * down ** (n - s) >= boundary:
                    crossed += amount
                else:
                    next_live[s] = next_live.get(s, F(0)) + amount
        live = next_live
    return crossed


def report(ratios, audits, rounds_left, allowance_left):
    current = sum(ratios.values(), F(0)) / len(ratios)
    maximum = envelope(ratios, audits)
    powers = {world: component_power(value, agreement_probability(world), audits, len(ratios))
              for world, value in ratios.items()}
    supported = [world for world, probability in powers.items() if probability >= POWER_TARGET]
    if current >= THRESHOLD:
        status = "evidence_met"
    elif maximum < THRESHOLD:
        status = "unreachable_with_remaining_budget"
    elif supported:
        status = "power_supported_for_named_models"
    else:
        status = "reachable_power_not_certified"
    return {"status": status, "evidence": float(current), "evidence_exact": str(current),
            "maximum_evidence": float(maximum), "maximum_evidence_exact": str(maximum),
            "available_audits": audits, "rounds_left": rounds_left,
            "audit_allowance_left": allowance_left,
            "component_likelihood_ratios": {w: str(v) for w, v in ratios.items()},
            "power_lower_bounds": {w: float(p) for w, p in powers.items()},
            "power_lower_bounds_exact": {w: str(p) for w, p in powers.items()},
            "supported_models": supported, "power_target": float(POWER_TARGET),
            "scope": "supplied alternatives; sufficient component-crossing certificates",
            "unrepresented_alternatives": "unknown", "absence_of_signal_certified": False}


class FeasibleObserver(base.Observer):
    def __init__(self, model, policy):
        super().__init__(model, depth=2)
        if policy not in POLICIES:
            raise ValueError("unknown audit policy")
        self.policy, self.audit_cap = policy, POLICIES[policy]
        self.ratios = {w: F(1) for w in model.family if w != "noise"}
        self.audit_used = self.round_number = 0
        self.first_stop = self.first_crossing = self.first_reopening = None
        self.max_evidence = F(1)
        self.audit = self.veto = False
        self.recommendation = self.observation = None
        self.snapshots = {}

    @property
    def evidence(self):
        return sum(self.ratios.values(), F(0)) / len(self.ratios)

    def available_audits(self, rounds_left):
        slots = rounds_left
        if self.policy == "periodic24":
            slots = sum(t % 8 == 0 for t in range(self.round_number,
                                                 self.round_number + rounds_left))
        return min(self.audit_cap - self.audit_used, slots)

    def state_report(self, rounds_left):
        return {"round": self.round_number, **report(self.ratios,
                self.available_audits(rounds_left), rounds_left, self.audit_cap - self.audit_used)}

    def choose(self, rounds_left):
        self.recommendation = super().choose(rounds_left)
        self.round_number += 1
        if self.recommendation == "stop" and self.first_stop is None:
            self.first_stop = self.round_number
            self.snapshots["first_stop"] = self.state_report(rounds_left)
        proposed, self.veto = gate.gated_action(self.recommendation, self.first_stop is not None,
                                               self.evidence, True)
        self.audit = False
        scheduled = self.policy != "periodic24" or self.round_number % 8 == 0
        candidate = (proposed == "stop" and self.evidence < THRESHOLD
                     and self.audit_used < self.audit_cap and scheduled)
        if candidate:
            reachable = envelope(self.ratios, self.available_audits(rounds_left)) >= THRESHOLD
            if self.policy == "feasible24" and not reachable:
                if "first_unreachable_refusal" not in self.snapshots:
                    self.snapshots["first_unreachable_refusal"] = self.state_report(rounds_left)
            else:
                proposed, self.audit = "both", True
                self.audit_used += 1
        self.action = proposed
        if base.CREDITS[base.ACTIONS.index(proposed)] > self.budget:
            raise ValueError("per-round budget exceeded")
        if self.audit_used > self.audit_cap:
            raise AssertionError("lifetime audit allowance exceeded")
        if self.first_stop is not None and proposed != "stop" and not self.audit:
            if self.evidence < THRESHOLD:
                raise AssertionError("post-stop unforced query lacks preceding evidence")
            if self.first_reopening is None:
                self.first_reopening = self.round_number
        return proposed

    def predict(self, observation):
        result = super().predict(observation)
        self.observation = observation
        return result

    def feedback(self, target):
        action, observation = self.action, self.observation
        super().feedback(target)
        for world in self.ratios:
            self.ratios[world] *= factor(world, action, observation, target)
        self.observation = None
        self.max_evidence = max(self.max_evidence, self.evidence)
        if self.evidence >= THRESHOLD and self.first_crossing is None:
            self.first_crossing = self.round_number

    def diagnostic(self):
        action, _ = self.model.choose(self.weights, depth=2)
        return gate.gated_action(action, self.first_stop is not None, self.evidence, True)[0]


def run_seed(model, policy, rows):
    observer = FeasibleObserver(model, policy)
    errors, prices, logs, actions, audits, vetoes = [], [], [], [], [], []
    for t, row in enumerate(rows):
        action = observer.choose(len(rows) - t)
        probability, prediction = observer.predict(base.acquire(row, action, observer.budget))
        target = row[2]
        errors.append(float(prediction != target))
        prices.append(base.COSTS[base.ACTIONS.index(action)])
        logs.append(-math.log(probability if target else 1 - probability))
        actions.append(action)
        audits.append(observer.audit)
        vetoes.append(observer.veto)
        observer.feedback(target)
    observer.snapshots["final"] = observer.state_report(0)
    mean = statistics.fmean
    return {"utility": mean(e + c for e, c in zip(errors, prices)),
            "error": mean(errors), "price": mean(prices), "log_loss": mean(logs),
            "tail_utility": mean(e + c for e, c in zip(errors[-16:], prices[-16:])),
            "action_counts": {a: actions.count(a) for a in base.ACTIONS},
            "tail_action_counts": {a: actions[-16:].count(a) for a in base.ACTIONS},
            "final_posterior": dict(zip(model.family, observer.weights)),
            "audit_count": observer.audit_used, "audit_cap": observer.audit_cap,
            "audit_price_per_round": sum(audits) * 0.12 / len(rows),
            "total_audit_spend": sum(audits) * 0.12,
            "audit_rounds": [t + 1 for t, x in enumerate(audits) if x],
            "gate_veto_rounds": [t + 1 for t, x in enumerate(vetoes) if x],
            "first_stop": observer.first_stop, "first_crossing": observer.first_crossing,
            "first_reopening": observer.first_reopening,
            "ever_post_stop_query": observer.first_reopening is not None,
            "final_unforced_recommendation": observer.diagnostic(),
            "evidence_final": float(observer.evidence), "evidence_max": float(observer.max_evidence),
            "reports": observer.snapshots}


def controls():
    inherited = gate.controls()
    checks = 0

    def check(condition, label):
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    for family in (base.WORLDS[:3], base.WORLDS):
        model = base.Model(family)
        alternatives = [w for w in family if w != "noise"]
        k = len(alternatives)
        for action, cells in model.events.items():
            for obs, v0, v1 in cells:
                for i, world in enumerate(family):
                    for y, vector in enumerate((v0, v1)):
                        check(abs(float(factor(world, action, obs, y)) -
                                  vector[i] / vector[family.index("noise")]) < base.TOL,
                              "exact tracker matches independently checked likelihood table")
        # Full joint histories, not independent fictitious agreements for each model.
        ratios = {w: F(1) for w in alternatives}
        ratios[alternatives[-1]] = F(25)
        for world in alternatives:
            crossing = F(0)
            events = tuple(itertools.product((0, 1), repeat=3))
            for path in itertools.product(events, repeat=3):
                current, mass, hit = dict(ratios), F(1), False
                for h, s, y in path:
                    mass *= F(1, 8) * factor(world, "both", (h, s), y)
                    for candidate in current:
                        current[candidate] *= factor(candidate, "both", (h, s), y)
                    hit |= sum(current.values()) / k >= THRESHOLD
                crossing += mass * hit
            lower = component_power(ratios[world], agreement_probability(world), 3, k)
            check(crossing >= lower, "joint mixture power dominates component certificate")
        for audits in (0, 1, 3, 8, 24):
            current = dict(ratios)
            for _ in range(audits):
                for w in current:
                    current[w] *= factor(w, "both", (0, 0), 0)
            check(sum(current.values()) / k == envelope(ratios, audits), "attainable joint envelope")
        for w in alternatives:
            p, start = agreement_probability(w), F(10)
            for cap in range(9):
                brute = F(0)
                for path in itertools.product((0, 1), repeat=cap):
                    successes = sum(path)
                    mass = p ** successes * (1 - p) ** (cap - successes)
                    current, crossed = start, start >= k * THRESHOLD
                    for outcome in path:
                        current *= 2 * (p if outcome else 1 - p)
                        crossed |= current >= k * THRESHOLD
                    brute += mass * crossed
                check(brute == component_power(start, p, cap, k), "exhaustive binary component power")
        for world in ("sensor", "parity"):
            rows = base.sample_world(world, 0)
            old = gate.run_seed(model, "gated", rows)
            new = run_seed(model, "periodic24", rows)
            for key in ("utility", "error", "price", "action_counts", "audit_rounds",
                        "first_stop", "first_reopening", "final_unforced_recommendation"):
                check(new[key] == old[key], "periodic control-seed regression")
        left, right = FeasibleObserver(model, "feasible24"), FeasibleObserver(model, "feasible24")
        for t, row in enumerate(base.sample_world("noise", 100)):
            before = dict(left.ratios)
            action = left.choose(64 - t)
            check(action == right.choose(64 - t), "identical public-state action")
            check(left.ratios == before, "no reset or current-target evidence")
            obs = base.acquire(row, action)
            check(left.predict(obs) == right.predict(obs), "identical prediction")
            left.feedback(row[2]); right.feedback(row[2])
            check(left.ratios == right.ratios, "identical evidence update")
            check(left.audit_used <= left.audit_cap, "audit budget held")
        check(left.state_report(0) == right.state_report(0), "identical report replay")
        try:
            left.feedback(0)
        except ValueError:
            check(left.ratios == right.ratios, "invalid feedback cannot change evidence")
        else:
            check(False, "double feedback must fail")
    for ratios, cap, expected in (({"parity": F(1, 10)}, 8, "unreachable_with_remaining_budget"),
                                  ({"parity": F(1)}, 8, "reachable_power_not_certified"),
                                  ({"parity": F(1)}, 24, "power_supported_for_named_models"),
                                  ({"parity": F(20)}, 0, "evidence_met")):
        r = report(ratios, cap, cap, cap)
        check(r["status"] == expected, "report status boundary")
        check(not r["absence_of_signal_certified"], "no closure claim")
    for policy, expected in (("periodic24", 8), ("burst8", 8), ("burst24", 24)):
        model = base.Model(base.WORLDS)
        result = run_seed(model, policy, [(0, 0, 1)] * 64)
        check(result["audit_count"] == expected, "fixed adverse stream exhausts only allowed audits")
        if policy == "periodic24":
            check(result["audit_rounds"] == list(range(8, 65, 8)), "periodic calendar cap")
    observer = FeasibleObserver(base.Model(base.WORLDS), "periodic24")
    observer.round_number = 7
    check(observer.available_audits(2) == 1, "calendar includes current and final opportunities")
    check(observer.available_audits(0) == 0, "no future time means no audit capacity")
    return {**inherited, "sa3f_checks": checks}


def summarize(records):
    summary = base.summarize(records)
    for key in ("audit_count", "audit_price_per_round"):
        summary["mean_" + key] = statistics.fmean(r[key] for r in records)
    summary["ever_post_stop_query_count"] = sum(r["ever_post_stop_query"] for r in records)
    summary["final_unforced_action_counts"] = {
        a: sum(r["final_unforced_recommendation"] == a for r in records) for a in base.ACTIONS}
    for snapshot in ("first_stop", "first_unreachable_refusal", "final"):
        statuses = {}
        for r in records:
            status = r["reports"].get(snapshot, {}).get("status", "not_applicable")
            statuses[status] = statuses.get(status, 0) + 1
        summary[snapshot + "_status_counts"] = statuses
    return summary


def cohort(worlds, seeds):
    results = {}
    for world in worlds:
        streams = [base.sample_world(world, seed) for seed in seeds]
        data = {"policies": {}, "paired_comparisons": {}}
        for family, names in (("restricted", base.WORLDS[:3]), ("expanded", base.WORLDS)):
            model = base.Model(names)
            for policy in POLICIES:
                records = [{"seed": seed, **run_seed(model, policy, rows)}
                           for seed, rows in zip(seeds, streams)]
                data["policies"][family + "_" + policy] = {"summary": summarize(records), "seeds": records}
            for first, second in (("burst8", "periodic24"), ("burst24", "periodic24"),
                                  ("feasible24", "burst24")):
                a, b = (data["policies"][family + "_" + p]["seeds"] for p in (first, second))
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
    sources = {"sa3f": Path(__file__), "sa3d": Path(gate.__file__),
               "sa3b": Path(base.__file__), "sa3c": Path(gate.previous.__file__)}
    result = {"protocol": "SA-3F", "registration_commit": "cf65fd7a",
              "source_hashes": {k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in sources.items()},
              "controls": controls()}
    if not args.controls_only:
        result["configuration"] = {"rounds": 64, "main_seeds": MAIN_SEEDS, "null_seeds": NULL_SEEDS,
                                   "audit_allowances": POLICIES, "threshold": 20, "power_target": 0.8}
        result["main_cohort"] = cohort(base.WORLDS, MAIN_SEEDS)
        result["null_cohort"] = cohort(("noise",), NULL_SEEDS)
    rendered = json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as output:
            output.write(rendered)
        print(json.dumps({"output": str(args.output), "controls": result["controls"],
                          "source_hashes": result["source_hashes"]}, sort_keys=True))
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
