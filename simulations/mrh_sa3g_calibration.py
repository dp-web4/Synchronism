#!/usr/bin/env python3
"""SA-3G: exact separation of null validity, alternative power, and reachability."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import mrh_sa3e_evidence_budget as previous
import mrh_sa3f_feasible_bursts as observer

REGISTRATION = "d274f576"
THRESHOLD = F(20)
PRICE = F(3, 25)
STARTS = (F(1, 10), F(1))
ASSUMED = (F(13, 20), F(41, 50), F(9, 10))
TRUE = (F(1, 2), F(13, 25), F(11, 20), F(3, 5), F(13, 20),
        F(3, 4), F(41, 50), F(9, 10))
RULES = {"nominal": F(1, 2), "guarded": F(11, 20)}
CAP = 64
INTERVAL = (F(3, 4), F(9, 10))
DEPENDENCIES = {
    "mrh_sa3e_evidence_budget.py": "8b9627af221e83e6b123aa41720eaf5aa05f3c86e36cb872854b7dff2c68f5c0",
    "mrh_sa3f_feasible_bursts.py": "0d1f805d0f1f8af8d8a74b7ac8d74ae5f5f59c8a814b52b42c3b06815e9d1534",
}


def factors(q, u):
    if not 0 < u < q < 1:
        raise ValueError("require 0 < null ceiling < assumed alternative < 1")
    return q / u, (1 - q) / (1 - u)


def evidence(start, up, down, n, successes):
    return start * up ** successes * down ** (n - successes)


@lru_cache(maxsize=256)
def curve(start, q, u, p, cap=CAP):
    """Exact IID first passage; audit until crossing or cap, without pruning."""
    up, down = factors(q, u)
    if start <= 0 or not 0 <= p <= 1 or type(cap) is not int or cap < 0:
        raise ValueError("invalid starting evidence, true rate, or cap")
    live = {0: F(1)} if start < THRESHOLD else {}
    crossed = F(start >= THRESHOLD)
    stopped_evidence = start * crossed
    expected_count = F(0)
    records = []
    for n in range(cap + 1):
        survival = sum(live.values(), F(0))
        terminal = stopped_evidence + sum(
            (mass * evidence(start, up, down, n, s) for s, mass in live.items()), F(0))
        records.append({"cap": n, "crossing_probability": crossed,
                        "survival_probability": survival,
                        "expected_audits": expected_count,
                        "expected_acquisition_spend": PRICE * expected_count,
                        "expected_stopped_evidence": terminal,
                        "maximum_evidence": start * up ** n})
        if n == cap:
            break
        expected_count += survival
        next_live = {}
        for s, mass in live.items():
            for bit, chance in ((1, p), (0, 1 - p)):
                amount = mass * chance
                if not amount:
                    continue
                new_s = s + bit
                value = evidence(start, up, down, n + 1, new_s)
                if value >= THRESHOLD:
                    crossed += amount
                    stopped_evidence += amount * value
                else:
                    next_live[new_s] = next_live.get(new_s, F(0)) + amount
        live = next_live
    return records


def brute_force(start, q, u, p, cap):
    """Independent path traversal, including outcomes after stopping in path mass."""
    crossing = count = terminal = F(0)
    for path in itertools.product((0, 1), repeat=cap):
        probability = p ** sum(path) * (1 - p) ** (cap - sum(path))
        value, spent = start, 0
        for bit in path:
            if value >= THRESHOLD:
                break
            value *= q / u if bit else (1 - q) / (1 - u)
            spent += 1
        crossing += probability * (value >= THRESHOLD)
        count += probability * spent
        terminal += probability * value
    return crossing, count, terminal


def adversarial_null_crossing(start, q, u, cap):
    """Bellman maximum over conditional p in [0,u], attained at an endpoint."""
    @lru_cache(maxsize=None)
    def visit(n, s):
        current = start * (q / u) ** s * ((1 - q) / (1 - u)) ** (n - s)
        if current >= THRESHOLD:
            return F(1)
        if n == cap:
            return F(0)
        on_one, on_zero = visit(n + 1, s + 1), visit(n + 1, s)
        return max(on_zero, u * on_one + (1 - u) * on_zero)
    return visit(0, 0)


def controls():
    count = 0

    def check(condition, label):
        nonlocal count
        if not condition:
            raise AssertionError(label)
        count += 1

    for start, q, u in itertools.product(STARTS, ASSUMED, RULES.values()):
        up, down = factors(q, u)
        check(up > 1 > down > 0, "coordinatewise increasing positive factors")
        check(u * up + (1 - u) * down == 1, "null boundary mean one")
        for p in (F(0), u / 2, u):
            check(p * up + (1 - p) * down <= 1, "conditional null mean at most one")
        for cap in range(9):
            check(adversarial_null_crossing(start, q, u, cap)
                  == curve(start, q, u, u)[cap]["crossing_probability"],
                  "adaptive null adversary equals IID upper boundary")
        for p in (*TRUE, F(0), F(1)):
            records = curve(start, q, u, p)
            old_probability = old_count = F(0)
            for r in records:
                n, crossing = r["cap"], r["crossing_probability"]
                check(crossing + r["survival_probability"] == 1, "probability conserved")
                check(old_probability <= crossing <= 1, "monotone cap probability")
                check(old_count <= r["expected_audits"] <= n, "monotone bounded audit count")
                check(r["expected_acquisition_spend"] == PRICE * r["expected_audits"],
                      "every audit priced")
                if r["maximum_evidence"] < THRESHOLD:
                    check(crossing == 0, "zero power before reachability")
                if p <= u:
                    check(crossing <= start / THRESHOLD, "applicable null crossing bound")
                    check(r["expected_stopped_evidence"] <= start, "null stopped supermartingale")
                if p == u:
                    check(r["expected_stopped_evidence"] == start, "boundary stopped martingale")
                if p in (0, 1):
                    check(crossing == int(p == 1 and r["maximum_evidence"] >= THRESHOLD),
                          "deterministic crossing endpoints")
                if n <= 8:
                    actual = brute_force(start, q, u, p, n)
                    check(actual == (crossing, r["expected_audits"], r["expected_stopped_evidence"]),
                          "independent exhaustive first passage, counts, moment")
                old_probability, old_count = crossing, r["expected_audits"]
        for low, high in zip(TRUE, TRUE[1:]):
            for a, b in zip(curve(start, q, u, low), curve(start, q, u, high)):
                check(a["crossing_probability"] <= b["crossing_probability"],
                      "crossing monotone in true agreement")
        for n in range(9):
            maximum = max(start * up ** sum(path) * down ** (n - sum(path))
                          for path in itertools.product((0, 1), repeat=n))
            check(maximum == start * up ** n, "exhaustive exact envelope")
        check(curve(F(20), q, u, F(1, 2), 0)[0]["expected_audits"] == 0,
              "threshold already met costs zero")
        check(curve(F(20), q, u, F(1, 2), 0)[0]["crossing_probability"] == 1,
              "threshold boundary inclusive")

    for start, p in itertools.product(STARTS, (F(1, 2), F(41, 50))):
        for old, new in zip(previous.curve(start, p), curve(start, F(41, 50), F(1, 2), p)):
            check(all(old[k] == new[k] for k in ("cap", "crossing_probability",
                      "survival_probability", "expected_audits", "expected_acquisition_spend")),
                  "SA-3E exact record regression")
            check(old["expected_terminal_evidence"] == new["expected_stopped_evidence"],
                  "SA-3E stopped moment regression")
    for start, q, n in itertools.product(STARTS, ASSUMED, range(25)):
        check(observer.component_power(start, q, n, 1)
              == curve(start, q, F(1, 2), q)[n]["crossing_probability"],
              "SA-3F single-component power regression")
    for filename, expected in DEPENDENCIES.items():
        check(hashlib.sha256(Path(__file__).with_name(filename).read_bytes()).hexdigest()
              == expected, "frozen source dependency")
    return {"checks_passed": count}


def number(value):
    return {"exact": str(value), "decimal": float(value)}


def serialize_record(record):
    return {k: v if k == "cap" else number(v) for k, v in record.items()}


def minimum_caps(records):
    return {str(target): next((r["cap"] for r in records
                              if r["crossing_probability"] >= target), None)
            for target in (F(1, 2), F(4, 5), F(9, 10))}


def experiment():
    curves, reports = [], []
    overclaims = {"tested_p_cases": 0, "cases_at_caps_8_24_64": []}
    violations = []
    for start, q, (rule, u) in itertools.product(STARTS, ASSUMED, RULES.items()):
        declared = curve(start, q, u, q)
        floor = curve(start, q, u, INTERVAL[0])
        identifying = {"start": str(start), "assumed_alternative": str(q),
                       "rule": rule, "null_ceiling": str(u)}
        reports.append({**identifying, "iid_alternative_interval": [str(p) for p in INTERVAL],
                        "interval_provenance": "supplied, not estimated or certified from data",
                        "assumed_model_minimum_caps": minimum_caps(declared),
                        "interval_floor_minimum_caps": minimum_caps(floor),
                        "caps": [{"cap": n, "assumed_model_power": number(declared[n]["crossing_probability"]),
                                  "interval_floor_power": number(floor[n]["crossing_probability"]),
                                  "assumed_model_supports_80_percent": declared[n]["crossing_probability"] >= F(4, 5),
                                  "entire_interval_supports_80_percent": floor[n]["crossing_probability"] >= F(4, 5)}
                                 for n in range(CAP + 1)]})
        for p in TRUE:
            records = curve(start, q, u, p)
            applicable = p <= u
            curves.append({**identifying, "true_iid_agreement": str(p),
                           "null_bound_applicable": applicable,
                           "claimed_null_bound": number(start / THRESHOLD),
                           "one_step_mean_factor": number(p * q / u + (1 - p) * (1 - q) / (1 - u)),
                           "minimum_caps": minimum_caps(records),
                           "caps": [serialize_record(r) for r in records]})
            for n, r in enumerate(records):
                if declared[n]["crossing_probability"] >= F(4, 5) and r["crossing_probability"] < F(4, 5):
                    overclaims["tested_p_cases"] += 1
                    if n in (8, 24, 64):
                        overclaims["cases_at_caps_8_24_64"].append(
                            {**identifying, "true_iid_agreement": str(p), "cap": n,
                             "actual_power": number(r["crossing_probability"]),
                             "assumed_model_power": number(declared[n]["crossing_probability"])})
                if p <= RULES["guarded"] and r["crossing_probability"] > start / THRESHOLD:
                    violations.append({**identifying, "true_iid_agreement": str(p), "cap": n,
                                       "crossing_probability": number(r["crossing_probability"]),
                                       "bound": number(start / THRESHOLD)})
    assert len(curves) == 96 and sum(len(c["caps"]) for c in curves) == 6240
    return {"curves": curves, "power_reports": reports,
            "nominal_power_overclaims_if_true_p_ignored": overclaims,
            "broader_null_cases_exceeding_inapplicable_nominal_bound": violations}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--output", type=Path, help="Create new artifact; refuse overwrite")
    args = parser.parse_args()
    result = {"protocol": "SA-3G", "registration_commit": REGISTRATION,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "dependency_sha256": {filename: hashlib.sha256(Path(__file__).with_name(filename).read_bytes()).hexdigest()
                                    for filename in DEPENDENCIES},
              "controls": controls()}
    if not args.controls_only:
        result["configuration"] = {"threshold": str(THRESHOLD), "audit_price": str(PRICE),
                                   "starts": [str(v) for v in STARTS], "max_cap": CAP,
                                   "assumed_alternatives": [str(v) for v in ASSUMED],
                                   "true_rates": [str(v) for v in TRUE],
                                   "null_ceilings": {k: str(v) for k, v in RULES.items()}}
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
