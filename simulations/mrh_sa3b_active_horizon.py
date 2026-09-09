#!/usr/bin/env python3
"""Registered SA-3B finite Bayesian acquisition experiment; stdlib only.

Policy code receives likelihood models and public feedback, never world identity
or latent rows. A history bit is a stored record, not a dynamical simulation.
See explorations/2026-09-08-codex-sa3b-active-horizon-charter.md.
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import statistics


WORLDS = ("memory", "sensor", "noise", "parity")
ACTIONS = ("stop", "history", "sensor", "both")
COSTS = (0.0, 0.04, 0.08, 0.12)
CREDITS = (0, 1, 2, 3)
ORACLE = dict(zip(WORLDS, (0.14, 0.26, 0.5, 0.3)))
ROUNDS = 64
SEEDS = tuple(range(24))
TOL = 1e-12


def acquire(row, action, budget=3):
    """Environment boundary: release only purchased fields, never the target."""
    if action not in ACTIONS or CREDITS[ACTIONS.index(action)] > budget:
        raise ValueError("invalid or over-budget acquisition")
    h, sensor, _ = row
    return {"stop": (), "history": (h,), "sensor": (sensor,),
            "both": (h, sensor)}[action]


def likelihood_tables():
    """Enumerate latent bits; return P(observation, Y | model, action)."""
    tables = {action: {} for action in ACTIONS}
    for h, s, n, m, u in itertools.product((0, 1), repeat=5):
        mass = 0.25 * (0.1 if n else 0.9) * (0.1 if m else 0.9) * 0.5
        targets = (h ^ n, s ^ n, u, h ^ s ^ n)
        for action in ACTIONS:
            obs = acquire((h, s ^ m, None), action)
            cell = tables[action].setdefault(obs, [[0.0] * 4, [0.0] * 4])
            for k, y in enumerate(targets):
                cell[y][k] += mass
    return tables


class Model:
    """Candidate family; no access to the evaluated world's identity."""

    def __init__(self, family):
        if not family or len(set(family)) != len(family):
            raise ValueError("empty or duplicate hypothesis family")
        indices = tuple(WORLDS.index(name) for name in family)
        self.family = tuple(family)
        self.events = {
            action: tuple((obs, tuple(vectors[0][k] for k in indices),
                           tuple(vectors[1][k] for k in indices))
                          for obs, vectors in sorted(cells.items()))
            for action, cells in likelihood_tables().items()
        }

    @staticmethod
    def mass(weights, vector):
        return sum(w * p for w, p in zip(weights, vector))

    @staticmethod
    def posterior(weights, vector):
        products = tuple(w * p for w, p in zip(weights, vector))
        total = sum(products)
        if total <= 0:
            raise ValueError("impossible feedback")
        return tuple(x / total for x in products)

    def stage_values(self, weights, budget=3, costs=COSTS):
        return {
            action: costs[i] + sum(min(self.mass(weights, v0),
                                      self.mass(weights, v1))
                                   for _, v0, v1 in self.events[action])
            for i, action in enumerate(ACTIONS) if CREDITS[i] <= budget
        }

    def choose(self, weights, depth=1, budget=3, costs=COSTS):
        if depth not in (1, 2):
            raise ValueError("only registered depths 1 and 2 are supported")
        values = self.stage_values(weights, budget, costs)
        if depth == 2:
            for action in values:
                for _, v0, v1 in self.events[action]:
                    for vector in (v0, v1):
                        mass = self.mass(weights, vector)
                        if mass:
                            post = self.posterior(weights, vector)
                            values[action] += mass * min(
                                self.stage_values(post, budget, costs).values())
        best = min(values.values())
        return next(a for a in ACTIONS if a in values and values[a] <= best + TOL), values


class Observer:
    """Stateful protocol enforces choose -> predict -> feedback sequencing."""

    def __init__(self, model, depth=1, fixed=None, budget=3):
        self.model = model
        self.weights = (1 / len(model.family),) * len(model.family)
        self.depth, self.fixed, self.budget = depth, fixed, budget
        self.phase = "choose"
        self.action = None
        self.pending = None

    def choose(self, rounds_left):
        if self.phase != "choose" or rounds_left < 1:
            raise ValueError("action requested out of sequence")
        action = self.fixed
        if action is None:
            action, _ = self.model.choose(self.weights, min(self.depth, rounds_left), self.budget)
        if action not in ACTIONS or CREDITS[ACTIONS.index(action)] > self.budget:
            raise ValueError("invalid or over-budget policy action")
        self.action, self.phase = action, "predict"
        return action

    def predict(self, observation):
        if self.phase != "predict":
            raise ValueError("prediction requested out of sequence")
        match = [cell for cell in self.model.events[self.action] if cell[0] == observation]
        if len(match) != 1:
            raise ValueError("observation does not match purchased channel")
        _, v0, v1 = match[0]
        m0, m1 = self.model.mass(self.weights, v0), self.model.mass(self.weights, v1)
        self.pending = (v0, v1)
        self.phase = "feedback"
        probability = m1 / (m0 + m1)
        return probability, int(probability > 0.5 + TOL)

    def feedback(self, target):
        if self.phase != "feedback" or target not in (0, 1):
            raise ValueError("feedback requested out of sequence or invalid target")
        self.weights = self.model.posterior(self.weights, self.pending[target])
        self.pending, self.action, self.phase = None, None, "choose"


def sample_world(world, seed):
    """Evaluator only; draws never depend on policy actions."""
    rng = random.Random(seed)
    rows = []
    for _ in range(ROUNDS):
        h, s = rng.randrange(2), rng.randrange(2)
        n, m = int(rng.random() < 0.1), int(rng.random() < 0.1)
        u = rng.randrange(2)
        y = {"memory": h ^ n, "sensor": s ^ n, "noise": u,
             "parity": h ^ s ^ n}[world]
        rows.append((h, s ^ m, y))
    return rows


def controls():
    checks = 0

    def check(condition, label):
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    def close(actual, expected, label):
        check(abs(actual - expected) < TOL, label)

    def rejects(function, label):
        try:
            function()
        except ValueError:
            check(True, label)
        else:
            check(False, label)

    model = Model(WORLDS)
    for action in ACTIONS:
        for k, world in enumerate(WORLDS):
            close(sum(v0[k] + v1[k] for _, v0, v1 in model.events[action]),
                  1, "likelihood normalization")
            for obs, v0, v1 in model.events[action]:
                close(v0[k] + v1[k], 1 / len(model.events[action]), "observation marginal")
                # Independent analytic formula, not a second call to enumeration.
                expected = 0.5
                if world == "memory" and action in ("history", "both"):
                    expected = 0.1 + 0.8 * obs[0]
                if world == "sensor" and action in ("sensor", "both"):
                    expected = 0.18 + 0.64 * obs[-1]
                if world == "parity" and action == "both":
                    expected = 0.18 + 0.64 * (obs[0] ^ obs[1])
                close(v1[k] / (v0[k] + v1[k]), expected, "conditional target law")
    for k, world in enumerate(WORLDS):
        weights = tuple(float(j == k) for j in range(4))
        for depth in (1, 2):
            action, values = model.choose(weights, depth)
            check(action == ("history", "sensor", "stop", "both")[k], "oracle action")
            close(values[action], depth * ORACLE[world], "oracle expected value")
        action, _ = model.choose(weights, budget=1)
        check(action == ("history" if world == "memory" else "stop"), "budget restriction")
        action, _ = model.choose(weights, costs=(0, 1, 1, 2))
        check(action == "stop", "priced-out channel")
    rejects(lambda: acquire((0, 1, 0), "both", 1), "over-budget environment")
    rejects(lambda: acquire((0, 1, 0), "unknown"), "invalid environment action")
    rejects(lambda: Observer(model, fixed="sensor", budget=1).choose(2), "over-budget observer")
    observer = Observer(model, fixed="stop")
    initial = observer.weights
    rejects(lambda: observer.feedback(1), "feedback before prediction")
    rejects(lambda: observer.predict(()), "prediction before action")
    for y in (0, 1, 1, 0):
        action = observer.choose(1)
        rejects(lambda: observer.feedback(y), "feedback before current prediction")
        rejects(lambda: observer.predict((1,)), "unbought observation rejected")
        obs = acquire((1, 0, y), action)
        check(obs == (), "stop reveals no fields")
        probability, prediction = observer.predict(obs)
        close(probability, 0.5, "stop predictive marginal")
        check(prediction == 0, "prediction tie")
        observer.feedback(y)
        for actual, expected in zip(observer.weights, initial):
            close(actual, expected, "labels without acquisition do not identify regime")
        rejects(lambda: observer.feedback(y), "double feedback")
    # The public transcript is identical; unpurchased rows may differ arbitrarily.
    left, right = Observer(model, depth=2), Observer(model, depth=2)
    for h, s, y in sample_world("memory", 100)[:16]:
        a = left.choose(16)
        check(a == right.choose(16), "replayed action")
        row = (h, s, y)
        alternative = (1 - h if a in ("stop", "sensor") else h,
                       1 - s if a in ("stop", "history") else s, 1 - y)
        obs = acquire(row, a)
        check(obs == acquire(alternative, a), "unbought values and targets hidden")
        check(left.predict(obs) == right.predict(obs), "replayed prediction")
        left.feedback(y)
        right.feedback(y)
        check(left.weights == right.weights, "replayed posterior")
    restricted = Model(WORLDS[:3])
    prior = (0.1, 0.1, 0.8)
    a1, v1 = restricted.choose(prior, 1)
    a2, v2 = restricted.choose(prior, 2)
    check(a1 == "stop", "registered break-even myopic probe")
    return {"checks_passed": checks,
            "planning_probe": {"prior": prior, "myopic_action": a1,
                               "myopic_values": v1, "two_step_action": a2,
                               "two_step_values": v2}}


def run_seed(model, depth, fixed, rows):
    observer = Observer(model, depth, fixed)
    errors, prices, logs, actions = [], [], [], []
    for t, row in enumerate(rows):
        action = observer.choose(len(rows) - t)
        observation = acquire(row, action, observer.budget)
        probability, prediction = observer.predict(observation)
        # Only here does the evaluator expose the target for scoring/feedback.
        y = row[2]
        errors.append(float(prediction != y))
        prices.append(COSTS[ACTIONS.index(action)])
        logs.append(-math.log(probability if y else 1 - probability))
        actions.append(action)
        observer.feedback(y)
    mean = statistics.fmean
    return {"utility": mean(e + c for e, c in zip(errors, prices)),
            "error": mean(errors), "price": mean(prices), "log_loss": mean(logs),
            "tail_utility": mean(e + c for e, c in zip(errors[-16:], prices[-16:])),
            "action_counts": {a: actions.count(a) for a in ACTIONS},
            "tail_action_counts": {a: actions[-16:].count(a) for a in ACTIONS},
            "final_posterior": dict(zip(model.family, observer.weights))}


def summarize(records):
    result = {}
    for key in ("utility", "error", "price", "log_loss", "tail_utility"):
        vals = [r[key] for r in records]
        result[key] = {"mean": statistics.fmean(vals),
                       "seed_standard_error": statistics.stdev(vals) / math.sqrt(len(vals))}
    for key, denominator in (("action_counts", ROUNDS), ("tail_action_counts", 16)):
        result[key.replace("counts", "frequencies")] = {
            a: statistics.fmean(r[key][a] / denominator for r in records) for a in ACTIONS}
    result["mean_final_posterior"] = {
        k: statistics.fmean(r["final_posterior"][k] for r in records)
        for k in records[0]["final_posterior"]}
    return result


def experiment():
    restricted, expanded = Model(WORLDS[:3]), Model(WORLDS)
    policies = {
        "restricted_myopic": (restricted, 1, None),
        "restricted_two_step": (restricted, 2, None),
        "expanded_myopic": (expanded, 1, None),
        "expanded_two_step": (expanded, 2, None),
        **{f"always_{a}": (expanded, 1, a) for a in ACTIONS},
    }
    results = {}
    for world in WORLDS:
        streams = [sample_world(world, seed) for seed in SEEDS]
        results[world] = {"oracle_expected_utility": ORACLE[world], "policies": {}}
        for name, (model, depth, fixed) in policies.items():
            records = [{"seed": seed, **run_seed(model, depth, fixed, rows)}
                       for seed, rows in zip(SEEDS, streams)]
            results[world]["policies"][name] = {"summary": summarize(records), "seeds": records}
        pairs = (("restricted_two_step", "restricted_myopic"),
                 ("expanded_two_step", "expanded_myopic"),
                 ("expanded_two_step", "restricted_two_step"))
        comparisons = {}
        for first, second in pairs:
            a, b = (results[world]["policies"][name]["seeds"] for name in (first, second))
            diffs = [x["utility"] - y["utility"] for x, y in zip(a, b)]
            comparisons[f"{first}_minus_{second}"] = {
                "mean": statistics.fmean(diffs),
                "seed_standard_error": statistics.stdev(diffs) / math.sqrt(len(diffs)),
                "paired_seed_differences": diffs}
        results[world]["paired_comparisons"] = comparisons
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--output", type=Path,
                        help="Create a new JSON result artifact; refuse to overwrite")
    args = parser.parse_args()
    result = {"protocol": "SA-3B", "registration_commit": "a09f253a",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "controls": controls()}
    if not args.controls_only:
        result["configuration"] = {"rounds": ROUNDS, "seeds": SEEDS,
                                   "costs": dict(zip(ACTIONS, COSTS)),
                                   "credits": dict(zip(ACTIONS, CREDITS)),
                                   "per_round_budget": 3}
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
