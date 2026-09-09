"""Contract API tests; no new acquisition policy or scientific evaluation."""

import copy
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import unittest

sys.dont_write_bytecode = True
import mrh_reopening_contract as contract


def state(ratio="1/10", audits=8):
    return {"component_likelihood_ratios": {"parity": ratio}, "available_audits": audits,
            "rounds_left": 32, "audit_allowance_left": 24}


class ReopeningContractTests(unittest.TestCase):
    def test_all_report_states(self):
        examples = (("1/10", 8, "unreachable_with_remaining_budget"),
                    ("1", 8, "reachable_power_not_certified"),
                    ("1", 24, "power_supported_for_named_models"),
                    ("20", 0, "evidence_met"))
        for ratio, audits, expected in examples:
            with self.subTest(expected=expected):
                result = contract.build_contract(state(ratio, audits))
                self.assertEqual(result["report"]["status"], expected)
                self.assertFalse(result["report"]["absence_of_signal_certified"])
                self.assertTrue(contract.verify_contract(result))

    def test_fraction_canonicalization(self):
        for ratio in ("2/20", Fraction(1, 10)):
            self.assertEqual(contract.build_contract(state(ratio)), contract.build_contract(state()))
        self.assertEqual(contract.build_contract(state(1))["state"]["component_likelihood_ratios"],
                         {"parity": "1"})

    def test_reject_invalid_ratios(self):
        for ratio in (0, -1, True, 0.1, "NaN", "Infinity", "1/0", "-1/2", "0",
                      "1e9", "1" * 257, 10 ** 1000, None, [], {}):
            with self.subTest(ratio_type=type(ratio).__name__):
                with self.assertRaises(ValueError):
                    contract.build_contract(state(ratio))

    def test_reject_invalid_budgets(self):
        for field, value in (("available_audits", 25), ("available_audits", -1),
                             ("available_audits", True), ("rounds_left", 65),
                             ("rounds_left", 7), ("audit_allowance_left", 7),
                             ("audit_allowance_left", 1.0)):
            item = state(); item[field] = value
            with self.subTest(field=field, value=value):
                with self.assertRaises(ValueError):
                    contract.build_contract(item)

    def test_reject_unknown_models_and_fields(self):
        for models in ({}, {"noise": "1"}, {"novel": "1"}, [], None):
            item = state(); item["component_likelihood_ratios"] = models
            with self.assertRaises(ValueError):
                contract.build_contract(item)
        item = state(); item["oracle_world"] = "parity"
        with self.assertRaises(ValueError):
            contract.build_contract(item)

    def test_reject_mutated_claims(self):
        original = contract.build_contract(state())
        mutations = (("status", "evidence_met"), ("absence_of_signal_certified", True),
                     ("absence_of_signal_certified", 0), ("supported_models", ["parity"]),
                     ("maximum_evidence_exact", "200"), ("available_audits", 24))
        for field, value in mutations:
            changed = copy.deepcopy(original); changed["report"][field] = value
            with self.subTest(field=field, value=value):
                with self.assertRaises(ValueError):
                    contract.verify_contract(changed)
        for field, value in (("schema", "unknown/2"), ("model", "unknown/2"),
                             ("assumptions", [])):
            changed = copy.deepcopy(original); changed[field] = value
            with self.assertRaises(ValueError):
                contract.verify_contract(changed)

    def test_state_changes_require_report_recomputation(self):
        item = contract.build_contract(state())
        item["state"]["component_likelihood_ratios"]["parity"] = "20"
        with self.assertRaises(ValueError):
            contract.verify_contract(item)

    def test_consistency_does_not_authenticate_evidence(self):
        # No observations establish this ratio. Recomputing a coherent report
        # still passes, documenting the interface's intentional trust boundary.
        fabricated = contract.build_contract(state("1000000", 0))
        self.assertTrue(contract.verify_contract(fabricated))
        self.assertEqual(fabricated["report"]["status"], "evidence_met")
        self.assertFalse(fabricated["report"]["absence_of_signal_certified"])

    def test_observer_report_compatibility(self):
        observer = contract.instrument.FeasibleObserver(
            contract.instrument.base.Model(contract.instrument.base.WORLDS), "feasible24")
        for t, row in enumerate(contract.instrument.base.sample_world("noise", 100)[:16]):
            action = observer.choose(64 - t)
            observer.predict(contract.instrument.base.acquire(row, action))
            observer.feedback(row[2])
        for snapshot in (observer.state_report(48), observer.state_report(0)):
            result = contract.build_contract({k: snapshot[k] for k in contract.STATE_KEYS})
            self.assertEqual(result["report"], {k: v for k, v in snapshot.items() if k != "round"})

    def cli(self, operation, payload):
        return subprocess.run([sys.executable, "-B", str(Path(contract.__file__)), operation],
                              input=payload, text=True, capture_output=True)

    def test_cli_round_trip(self):
        built = self.cli("build", json.dumps(state()))
        self.assertEqual(built.returncode, 0, built.stderr)
        verified = self.cli("verify", built.stdout)
        self.assertEqual(verified.returncode, 0, verified.stderr)
        self.assertTrue(json.loads(verified.stdout)["valid"])

    def test_cli_rejects_bad_wire_inputs(self):
        for payload in ('{"x": 1, "x": 2}', '{"x": NaN}', '[]', '{', ' ' * 65537):
            result = self.cli("build", payload)
            with self.subTest(payload_length=len(payload)):
                self.assertEqual(result.returncode, 2)
                self.assertEqual(result.stdout, "")
                self.assertIn("invalid reopening contract", result.stderr)


if __name__ == "__main__":
    unittest.main()
