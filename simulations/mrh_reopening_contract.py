#!/usr/bin/env python3
"""Checked JSON handoff for SA-3F reports; stdin/stdout, no experiment execution.

Checks internal mathematical consistency, NOT the provenance of supplied
evidence, truth of candidate models, sensor calibration, or a real calendar.
"""

import argparse
from fractions import Fraction
import json
import re
import sys

sys.dont_write_bytecode = True
import mrh_sa3f_feasible_bursts as instrument

SCHEMA = "mrh.reopening-contract/1"
MODEL = "sa3b-known-channels/1"
STATE_KEYS = {"component_likelihood_ratios", "available_audits", "rounds_left",
              "audit_allowance_left"}
ASSUMPTIONS = [
    "Stationary supplied memory/sensor/parity alternatives against independent noise.",
    "Paired audits have the declared calibrated likelihoods and cost 0.12 each.",
    "Target feedback arrives after prediction; available audit opportunities are supplied.",
    "Likelihood ratios are inherited evidence, not permission to reset a test's error budget.",
    "Power bounds are conditional on named alternatives, not probabilities those alternatives are true.",
    "Validation does not attest evidence provenance, model completeness, or absence of signal.",
]


def _count(value, name, maximum):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [0, {maximum}]")
    return value


def _ratio(value):
    if type(value) is int:
        # Bound before string conversion as well as before exact arithmetic.
        if value.bit_length() > 850:
            raise ValueError("likelihood ratio exceeds the supported numeric size")
        value = str(value)
    elif isinstance(value, Fraction):
        if max(value.numerator.bit_length(), value.denominator.bit_length()) > 850:
            raise ValueError("likelihood ratio exceeds the supported numeric size")
        value = str(value)
    if not isinstance(value, str) or not re.fullmatch(r"[0-9]{1,256}(?:/[0-9]{1,256})?", value):
        raise ValueError("likelihood ratios must be positive integers or exact fraction strings")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("invalid likelihood ratio") from exc
    if result <= 0:
        raise ValueError("likelihood ratios must be strictly positive for these models")
    return result


def build_contract(state):
    """Compute a canonical report from supplied state; does not authenticate it."""
    if type(state) is not dict or set(state) != STATE_KEYS:
        raise ValueError("state must contain exactly the four documented fields")
    inputs = state["component_likelihood_ratios"]
    if type(inputs) is not dict or not inputs or set(inputs) - {"memory", "sensor", "parity"}:
        raise ValueError("supply one or more known non-noise hypotheses: memory, sensor, parity")
    ratios = {name: _ratio(value) for name, value in sorted(inputs.items())}
    audits = _count(state["available_audits"], "available_audits", 24)
    rounds = _count(state["rounds_left"], "rounds_left", 64)
    allowance = _count(state["audit_allowance_left"], "audit_allowance_left", 24)
    if audits > min(rounds, allowance):
        raise ValueError("available audits exceed remaining time or audit allowance")
    canonical_state = {"component_likelihood_ratios": {k: str(v) for k, v in ratios.items()},
                       "available_audits": audits, "rounds_left": rounds,
                       "audit_allowance_left": allowance}
    return {"schema": SCHEMA, "model": MODEL, "state": canonical_state,
            "assumptions": list(ASSUMPTIONS),
            "report": instrument.report(ratios, audits, rounds, allowance)}


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def verify_contract(contract):
    """Return True for consistent canonical JSON, otherwise raise ValueError.

    A fabricated but self-consistent state passes: provenance is out of scope.
    """
    if type(contract) is not dict or set(contract) != {"schema", "model", "state", "assumptions", "report"}:
        raise ValueError("invalid contract envelope")
    if contract["schema"] != SCHEMA or contract["model"] != MODEL:
        raise ValueError("unsupported schema or model")
    expected = build_contract(contract["state"])
    # Compare typed canonical JSON, not Python equality (where False == 0).
    if _canonical(contract) != _canonical(expected):
        raise ValueError("contract is noncanonical or inconsistent with its supplied state")
    return True


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON field: {key}")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("build", "verify"))
    args = parser.parse_args()
    try:
        payload = sys.stdin.read(65537)
        if len(payload) > 65536:
            raise ValueError("input exceeds 65536 characters")
        value = json.loads(payload, object_pairs_hook=_unique_object, parse_constant=_reject_constant)
        if args.operation == "build":
            result = build_contract(value)
        else:
            verify_contract(value)
            result = {"valid": True, "status": value["report"]["status"],
                      "verification_scope": "internal consistency, not evidence provenance"}
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))
    except (ValueError, TypeError, OverflowError, RecursionError) as exc:
        print(f"invalid reopening contract: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
