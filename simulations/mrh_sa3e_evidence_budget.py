#!/usr/bin/env python3
"""SA-3E exact first-passage power/cost map for a supplied parity alternative."""

import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path

THRESHOLD = F(20)
STARTS = (F(1, 100), F(1, 10), F(1), F(10))
UP, DOWN = F(41, 25), F(9, 25)
WORLD_PROBABILITIES = {"noise": F(1, 2), "parity": F(41, 50)}
PRICE = F(3, 25)
MAX_CAP = 32


def evidence(start, n, successes):
    return start * UP ** successes * DOWN ** (n - successes)


def minimum_possible_audits(start):
    n, current = 0, start
    while current < THRESHOLD:
        current *= UP
        n += 1
    return n


def curve(start, probability, cap=MAX_CAP):
    if not 0 < start < THRESHOLD or not 0 <= probability <= 1 or cap < 0:
        raise ValueError("invalid start, probability, or cap")
    live = {0: F(1)}
    absorbed = stopped_evidence = expected_count = F(0)
    records = []
    for n in range(cap + 1):
        survival = sum(live.values(), F(0))
        terminal_evidence = stopped_evidence + sum(
            (mass * evidence(start, n, s) for s, mass in live.items()), F(0))
        records.append({"cap": n, "crossing_probability": absorbed,
                        "survival_probability": survival,
                        "expected_audits": expected_count,
                        "expected_acquisition_spend": PRICE * expected_count,
                        "expected_terminal_evidence": terminal_evidence})
        if n == cap:
            break
        expected_count += survival
        next_live = {}
        for successes, mass in live.items():
            for outcome, chance in ((1, probability), (0, 1 - probability)):
                new_mass = mass * chance
                if not new_mass:
                    continue
                new_s = successes + outcome
                new_e = evidence(start, n + 1, new_s)
                if new_e >= THRESHOLD:
                    absorbed += new_mass
                    stopped_evidence += new_mass * new_e
                else:
                    next_live[new_s] = next_live.get(new_s, F(0)) + new_mass
        live = next_live
    return records


def brute_force(start, probability, cap):
    crossing = audits = terminal = F(0)
    for path in itertools.product((0, 1), repeat=cap):
        successes = sum(path)
        mass = probability ** successes * (1 - probability) ** (cap - successes)
        current, spent, crossed = start, 0, False
        for outcome in path:
            current *= UP if outcome else DOWN
            spent += 1
            if current >= THRESHOLD:
                crossed = True
                break
        crossing += mass * crossed
        audits += mass * spent
        terminal += mass * current
    return crossing, audits, terminal


def controls():
    checks = 0

    def check(condition, label):
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    check((UP + DOWN) / 2 == 1, "null mean-one multiplier")
    for p in WORLD_PROBABILITIES.values():
        check(p + (1 - p) == 1, "outcome normalization")
    for start in STARTS:
        minimum = minimum_possible_audits(start)
        check(start * UP ** minimum >= THRESHOLD, "all-agreement reachability")
        check(start * UP ** (minimum - 1) < THRESHOLD, "reachability minimum")
        for world, p in WORLD_PROBABILITIES.items():
            records = curve(start, p)
            old_crossing = old_count = F(0)
            for r in records:
                check(r["crossing_probability"] + r["survival_probability"] == 1,
                      "probability conserved")
                check(0 <= r["expected_audits"] <= r["cap"], "audit count bounds")
                check(r["expected_audits"] >= old_count, "monotone expected count")
                check(r["crossing_probability"] >= old_crossing, "monotone crossing probability")
                if r["cap"] < minimum:
                    check(r["crossing_probability"] == 0, "unreachable cap has zero power")
                if world == "noise":
                    check(r["expected_terminal_evidence"] == start, "bounded optional stopping exact")
                    check(r["crossing_probability"] <= start / THRESHOLD, "conditional Ville bound")
                if r["cap"] <= 8:
                    crossing, audits, terminal = brute_force(start, p, r["cap"])
                    check(crossing == r["crossing_probability"], "brute-force first passage")
                    check(audits == r["expected_audits"], "brute-force audit count")
                    check(terminal == r["expected_terminal_evidence"], "brute-force terminal evidence")
                old_crossing, old_count = r["crossing_probability"], r["expected_audits"]
        for p in (F(0), F(1)):
            for r in curve(start, p):
                crosses = p == 1 and r["cap"] >= minimum
                count = min(r["cap"], minimum) if p == 1 else r["cap"]
                check(r["crossing_probability"] == int(crosses), "deterministic crossing endpoint")
                check(r["expected_audits"] == count, "deterministic cost endpoint")
    return {"checks_passed": checks}


def number(value):
    return {"exact": str(value), "decimal": float(value)}


def experiment():
    results = {}
    for start in STARTS:
        curves = {world: curve(start, p) for world, p in WORLD_PROBABILITIES.items()}
        results[str(start)] = {
            "starting_evidence": number(start),
            "conditional_null_bound": number(start / THRESHOLD),
            "minimum_possible_audits": minimum_possible_audits(start),
            "minimum_cap_for_power": {
                str(power): next((r["cap"] for r in curves["parity"]
                                  if r["crossing_probability"] >= power), None)
                for power in (F(1, 2), F(4, 5), F(9, 10))},
            "curves": {world: [{k: (v if k == "cap" else number(v)) for k, v in r.items()}
                               for r in records] for world, records in curves.items()},
        }
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--output", type=Path, help="Create a new artifact, refusing overwrite")
    args = parser.parse_args()
    result = {"protocol": "SA-3E", "registration_commit": "826005b9",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "controls": controls()}
    if not args.controls_only:
        result["configuration"] = {"threshold": str(THRESHOLD), "starts": [str(s) for s in STARTS],
                                   "max_cap": MAX_CAP, "audit_price": str(PRICE),
                                   "null_agreement_probability": "1/2",
                                   "alternative_agreement_probability": "41/50"}
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
