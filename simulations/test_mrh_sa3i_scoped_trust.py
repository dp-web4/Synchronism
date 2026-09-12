"""Independent task gold, counterfactual exposure, and failure-accounting pins."""

from collections import Counter
import copy
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import mrh_sa3i_scoped_trust as instrument


def independent_gold(task):
    if task.domain == "checksum":
        a, b, c, d = task.values
        return int((a + b + b + c + c + c + d + d + d + d) % 5 == task.residue)
    covered = [set(range(lo, hi)) & set(range(8)) for lo, hi in task.values]
    exactly_two = ((covered[0] & covered[1]) | (covered[1] & covered[2]) |
                   (covered[0] & covered[2])) - (covered[0] & covered[1] & covered[2])
    return int(len(exactly_two) >= 2)


class ScopedTrustTests(unittest.TestCase):
    def setUp(self):
        self.episodes = instrument.episodes()
        self.baselines = {task.id: {"answer": task.gold()} for task in instrument.TASKS}

    def decision(self, episode, answer=None, verify=False):
        return {"action": "verify" if verify else "answer", "answer": answer,
                "receipt": instrument.expected_receipt(episode)}

    def test_base_answers_independent(self):
        self.assertEqual([t.gold() for t in instrument.TASKS], [1, 0, 1, 0])
        for task in instrument.TASKS:
            self.assertEqual(task.gold(), independent_gold(task))

    def test_full_small_input_gold_enumeration(self):
        for digits in itertools.product(range(3), repeat=4):
            for residue in range(5):
                task = instrument.Task("test", "checksum", digits, residue)
                self.assertEqual(task.gold(), independent_gold(task))
        intervals = [(lo, hi) for lo in range(4) for hi in range(lo + 1, 5)]
        for values in itertools.product(intervals, repeat=3):
            task = instrument.Task("test", "interval", values)
            self.assertEqual(task.gold(), independent_gold(task))

    def test_history_gold_no_task_overlap(self):
        base_questions = {t.question() for t in instrument.TASKS}
        for tasks in instrument.histories().values():
            self.assertEqual(len(tasks), 8)
            self.assertEqual(len({t.question() for t in tasks}), 8)
            self.assertTrue(base_questions.isdisjoint(t.question() for t in tasks))
            for task in tasks:
                self.assertEqual(task.gold(), independent_gold(task))
                self.assertIn(f"{task.id}: {task.question()} Adviser answer: {task.gold()}; "
                              f"verified answer: {task.gold()}.", instrument.history_text())

    def test_fully_crossed_and_balanced(self):
        self.assertEqual(len(self.episodes), 64)
        self.assertEqual(len({e.id for e in self.episodes}), 64)
        grid = Counter((e.task_id, e.condition, e.alias, e.advice) for e in self.episodes)
        self.assertEqual(set(grid), set(itertools.product(instrument.TASK_BY_ID,
                         instrument.CONDITIONS, instrument.ALIASES, (0, 1))))
        self.assertEqual(set(grid.values()), {1})
        self.assertEqual(set(Counter((e.condition, e.advice == instrument.TASK_BY_ID[e.task_id].gold())
                                     for e in self.episodes).values()), {8})

    def test_pair_content_identical_except_attribution_and_id(self):
        index = {(e.task_id, e.alias, e.advice, e.condition): e for e in self.episodes}
        for task, alias, advice, baseline in itertools.product(instrument.TASKS,
                                                             instrument.ALIASES, (0, 1), (0, 1)):
            normalized = []
            for condition in instrument.CONDITIONS:
                prompt = instrument.advice_prompt(index[task.id, alias, advice, condition], baseline)
                normalized.append("\n".join(line for line in prompt.splitlines()
                                             if not line.startswith(("Episode:", "Adviser:"))))
            self.assertEqual(len(set(normalized)), 1)

    def test_anonymous_alias_blocks_are_duplicates_not_new_content(self):
        for task, advice in itertools.product(instrument.TASKS, (0, 1)):
            rows = [e for e in self.episodes if e.task_id == task.id and e.advice == advice
                    and e.condition == "anonymous"]
            texts = [instrument.advice_prompt(e, 0).split("\n", 1)[1] for e in rows]
            self.assertEqual(texts[0], texts[1])

    def test_complete_histories_and_no_evaluator_labels_in_advice(self):
        for episode in self.episodes:
            for baseline in (0, 1):
                prompt = instrument.advice_prompt(episode, baseline)
                self.assertEqual(prompt.count(instrument.history_text()), 1)
                self.assertIn(f"Your previously recorded answer: {baseline}.", prompt)
                self.assertIn(f"Suggested answer: {episode.advice}.", prompt)
                for forbidden in ("advice_correct", "expected_receipt", '"gold"',
                                  "Perfect reference answer:", "named_unlinked", "mismatched"):
                    self.assertNotIn(forbidden, prompt)
                self.assertRegex(episode.id, r"^[0-9a-f]{20}$")
        for task in instrument.TASKS:
            self.assertNotIn("record", instrument.baseline_prompt(task))
            self.assertNotIn("Suggested", instrument.baseline_prompt(task))

    def test_attribution_receipts(self):
        for episode in self.episodes:
            receipt = instrument.expected_receipt(episode)
            domain = instrument.TASK_BY_ID[episode.task_id].domain
            if episode.condition == "matched":
                self.assertEqual(receipt["linked_domain"], domain)
            elif episode.condition == "mismatched":
                self.assertNotEqual(receipt["linked_domain"], domain)
                self.assertIn(receipt["linked_domain"], ("checksum", "interval"))
            else:
                self.assertEqual(receipt, {"linked_domain": "none", "successes": 0, "total": 0})

    def test_verification_delivers_gold_after_purchase(self):
        for episode in self.episodes:
            prompt = instrument.verification_prompt(episode, 0)
            self.assertTrue(prompt.startswith(instrument.advice_prompt(episode, 0)))
            self.assertIn(f"Perfect reference answer: {independent_gold(instrument.TASK_BY_ID[episode.task_id])}.",
                          prompt)

    def test_scripted_control_totals(self):
        expected = {"verify": (64, 64, 12.8), "follow": (32, 0, 32), "keep": (64, 0, 0),
                    "scoped": (56, 48, 17.6), "unscoped": (48, 32, 22.4)}
        for policy, (correct, verified, loss) in expected.items():
            report = instrument.score_run(self.baselines, instrument.scripted_records(policy, self.baselines))
            self.assertEqual(report["status"], "complete")
            self.assertEqual(sum(r["correct"] for r in report["rows"]), correct)
            self.assertEqual(sum(r["verified"] for r in report["rows"]), verified)
            self.assertAlmostEqual(sum(r["loss"] for r in report["rows"]), loss)
            self.assertEqual(len(report["cells"]), 8)
            self.assertEqual(len(report["paired_contrasts"]), 64)

    def test_scripted_scoped_vs_unscoped_contrast(self):
        for policy, expected in (("scoped", 0), ("unscoped", 1)):
            report = instrument.score_run(self.baselines, instrument.scripted_records(policy, self.baselines))
            pairs = [p for p in report["paired_contrasts"] if p["left"] == "mismatched"
                     and p["right"] == "named_unlinked" and not p["advice_correct"]]
            self.assertEqual(len(pairs), 8)
            self.assertEqual({p["differences"]["corruption"] for p in pairs}, {expected})

    def test_transition_eligibility_and_empty_denominators(self):
        for episode in self.episodes:
            gold = instrument.TASK_BY_ID[episode.task_id].gold()
            for baseline in (0, 1):
                row = instrument.score_episode(episode, baseline, self.decision(episode, episode.advice))
                self.assertEqual(row["switch"], 1 if baseline != episode.advice else None)
                self.assertEqual(row["rescue"], 1 if baseline != gold and episode.advice == gold else None)
                self.assertEqual(row["corruption"], 1 if baseline == gold and episode.advice != gold else None)
        rows = instrument.score_run(self.baselines, instrument.scripted_records("follow", self.baselines))["rows"]
        self.assertEqual(instrument.metric(rows, "rescue"), {"sum": 0, "n": 0, "mean": None})
        self.assertEqual(instrument.metric(rows, "corruption"), {"sum": 32, "n": 32, "mean": 1})

    def test_wrong_receipt_not_excluded(self):
        episode = next(e for e in self.episodes if e.condition == "matched")
        decision = self.decision(episode, instrument.TASK_BY_ID[episode.task_id].gold())
        decision["receipt"] = {"linked_domain": "none", "successes": 0, "total": 0}
        row = instrument.score_episode(episode, 0, decision)
        self.assertTrue(row["valid"])
        self.assertFalse(row["receipt_ok"])
        self.assertEqual(row["correct"], 1)

    def test_malformed_answers_pay_failure_loss(self):
        episode = self.episodes[0]
        decisions = [None, [], "```json\n{}\n```", {}, {"action": "guess"},
                     self.decision(episode, True), self.decision(episode, 1.0),
                     self.decision(episode, 2), self.decision(episode, None),
                     '{"action":"answer","action":"verify"}']
        for decision in decisions:
            row = instrument.score_episode(episode, 0, decision)
            self.assertFalse(row["valid"])
            self.assertEqual(row["correct"], 0)
            self.assertEqual(row["loss"], 1)

    def test_failed_verification_still_costs(self):
        episode = self.episodes[0]
        decision = self.decision(episode, verify=True)
        for followup in (None, {"answer": True}, {"answer": 0, "extra": 1}):
            row = instrument.score_episode(episode, 0, decision, followup)
            self.assertFalse(row["valid"])
            self.assertEqual(row["verified"], 1)
            self.assertEqual(row["loss"], 1.2)
        row = instrument.score_episode(episode, 0, decision,
                                       {"answer": 1 - instrument.TASK_BY_ID[episode.task_id].gold()})
        self.assertTrue(row["valid"])
        self.assertEqual(row["loss"], 1.2)
        with self.assertRaises(ValueError):
            instrument.score_episode(episode, 0, self.decision(episode, 0), {"answer": 0})

    def test_bad_receipt_schema_and_extra_decision_keys(self):
        episode = self.episodes[0]
        for receipt in ([], {}, {"linked_domain": "none", "successes": True, "total": 0},
                        {"linked_domain": "none", "successes": 8, "total": 0}):
            decision = self.decision(episode, 0)
            decision["receipt"] = receipt
            self.assertFalse(instrument.score_episode(episode, 0, decision)["valid"])
        decision = self.decision(episode, 0)
        decision["confidence"] = 1
        self.assertFalse(instrument.score_episode(episode, 0, decision)["valid"])

    def test_missing_duplicate_unknown_and_partial_records(self):
        records = instrument.scripted_records("follow", self.baselines)
        partial = instrument.score_run(self.baselines, records[:-1])
        self.assertEqual(partial["status"], "incomplete")
        self.assertEqual(partial["missing_episode_ids"], [records[-1]["episode_id"]])
        self.assertEqual(partial["observed"], 63)
        self.assertNotIn("paired_contrasts", partial)
        with self.assertRaises(ValueError):
            instrument.score_run(self.baselines, records + records[:1])
        bad = copy.deepcopy(records)
        bad[0]["episode_id"] = "unknown"
        with self.assertRaises(ValueError):
            instrument.score_run(self.baselines, bad)
        records[0]["decision"] = "invalid"
        report = instrument.score_run(self.baselines, records)
        self.assertEqual(report["observed"], 64)
        self.assertEqual(sum(r["valid"] for r in report["rows"]), 63)

    def test_invalid_baselines_rejected_no_replacement(self):
        for baselines in ({}, {**self.baselines, "C1": {"answer": True}},
                          {**self.baselines, "C1": '{"answer":0,"answer":1}'}):
            with self.assertRaises(ValueError):
                instrument.score_run(baselines, [])

    def test_replay_and_source_hash(self):
        first = instrument.canonical(instrument.validation_report())
        self.assertEqual(first, instrument.canonical(instrument.validation_report()))
        report = json.loads(first)
        self.assertEqual(report["source_sha256"], hashlib.sha256(Path(instrument.__file__).read_bytes()).hexdigest())
        self.assertEqual(report["fixture_sha256"], instrument.digest(instrument.canonical(instrument.fixture_artifact())))
        actual = subprocess.check_output([sys.executable, instrument.__file__], text=True)
        self.assertEqual(first, actual)

    def test_cli_never_overwrites_output(self):
        with tempfile.TemporaryDirectory(prefix="sa3i-test-") as directory:
            target = Path(directory) / "report.json"
            command = [sys.executable, instrument.__file__, "--output", str(target)]
            subprocess.run(command, check=True, capture_output=True)
            saved = target.read_bytes()
            retry = subprocess.run(command, capture_output=True)
            self.assertNotEqual(retry.returncode, 0)
            self.assertEqual(saved, target.read_bytes())

    def test_saved_validation_artifact(self):
        saved = Path(instrument.__file__).with_name("mrh_sa3i_scoped_trust_validation.json")
        self.assertEqual(saved.read_text(), instrument.canonical(instrument.validation_report()))

    def test_cli_scores_raw_records_and_exports_fixture(self):
        records = instrument.scripted_records("follow", self.baselines)
        # Exercise strict JSON responses as well as the programmatic dict interface.
        for record in records:
            record["decision"] = json.dumps(record["decision"])
        with tempfile.TemporaryDirectory(prefix="sa3i-raw-test-") as directory:
            source = Path(directory) / "raw.json"
            source.write_text(json.dumps({"baselines": self.baselines, "records": records}))
            output = subprocess.check_output([sys.executable, instrument.__file__,
                                              "--score-input", str(source)], text=True)
            self.assertEqual(json.loads(output), instrument.score_run(self.baselines, records))
        fixture = subprocess.check_output([sys.executable, instrument.__file__, "--fixture"], text=True)
        self.assertEqual(json.loads(fixture), json.loads(instrument.canonical(instrument.fixture_artifact())))

    def test_controls_detect_broken_instruments(self):
        with patch.object(instrument.Task, "gold", return_value=0):
            with self.assertRaises(AssertionError):
                self.test_base_answers_independent()
        with patch.object(instrument, "history_text", return_value=""):
            with self.assertRaises(AssertionError):
                self.test_complete_histories_and_no_evaluator_labels_in_advice()
        with patch.object(instrument, "PRICE", 0):
            with self.assertRaises(AssertionError):
                self.test_scripted_control_totals()


if __name__ == "__main__":
    unittest.main()
