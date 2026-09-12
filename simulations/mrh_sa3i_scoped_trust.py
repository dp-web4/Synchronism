#!/usr/bin/env python3
"""SA-3I deterministic fixtures and scorer; no model transport or network calls."""

import argparse
from dataclasses import asdict, dataclass
import hashlib
import itertools
import json
from pathlib import Path
import random

REGISTRATION = "b480499b"
SEED = 20260912
PRICE = 0.2
CONDITIONS = ("matched", "mismatched", "named_unlinked", "anonymous")
ALIASES = ("Aster", "Beryl")


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def digest(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def binary(value):
    if type(value) is not int or value not in (0, 1):
        raise ValueError("answer must be integer 0 or 1, not a boolean")
    return value


@dataclass(frozen=True)
class Task:
    id: str
    domain: str
    values: tuple
    residue: int = 0

    def gold(self):
        if self.domain == "checksum":
            return int(sum(i * d for i, d in enumerate(self.values, 1)) % 5 == self.residue)
        if self.domain == "interval":
            counts = [sum(lo <= x < hi for lo, hi in self.values) for x in range(8)]
            return int(counts.count(2) >= 2)
        raise ValueError("unknown domain")

    def question(self):
        if self.domain == "checksum":
            return (f"Checksum task: digits {list(self.values)}. Multiply each digit by its "
                    "one-based position, sum, then take modulo 5. "
                    f"Is the result {self.residue}? Answer 1 for yes, 0 for no.")
        if self.domain == "interval":
            intervals = ", ".join(f"[{lo},{hi})" for lo, hi in self.values)
            return (f"Interval task: {intervals}. Intervals include their start and exclude "
                    "their end. Among integers 0 through 7, are at least two covered by "
                    "exactly two intervals? Answer 1 for yes, 0 for no.")
        raise ValueError("unknown domain")


TASKS = (
    Task("C1", "checksum", (3, 1, 4, 1), 1),
    Task("C2", "checksum", (2, 7, 1, 8), 2),
    Task("I1", "interval", ((0, 4), (2, 6), (3, 5))),
    Task("I2", "interval", ((0, 2), (3, 5), (6, 8))),
)
TASK_BY_ID = {task.id: task for task in TASKS}


def histories():
    checksum, interval = [], []
    for j in range(8):
        values = (j, 2, 5, 3)
        residue = (sum(i * d for i, d in enumerate(values, 1)) + j % 2) % 5
        checksum.append(Task(f"P{j}", "checksum", values, residue))
        interval.append(Task(f"Q{j}", "interval", ((0, j + 1), (1, 5), (3, 7))))
    return {"P": tuple(checksum), "Q": tuple(interval)}


def history_text():
    lines = ["BEGIN VERIFIED SYNTHETIC HISTORIES"]
    for record, tasks in histories().items():
        lines.append(f"Record {record}; domain {tasks[0].domain}; 8 successes out of 8.")
        for task in tasks:
            lines.append(f"{task.id}: {task.question()} Adviser answer: {task.gold()}; "
                         f"verified answer: {task.gold()}.")
    lines.append("END VERIFIED SYNTHETIC HISTORIES")
    return "\n".join(lines)


@dataclass(frozen=True)
class Episode:
    id: str
    task_id: str
    condition: str
    alias: str
    advice: int


def episodes():
    rows = []
    for task, condition, alias, correct in itertools.product(TASKS, CONDITIONS, ALIASES, (False, True)):
        identity = f"{SEED}|{task.id}|{condition}|{alias}|{correct}"
        rows.append(Episode(digest(identity)[:20], task.id, condition, alias,
                            task.gold() if correct else 1 - task.gold()))
    random.Random(SEED).shuffle(rows)
    return tuple(rows)


def linked_record(episode):
    matching = "P" if TASK_BY_ID[episode.task_id].domain == "checksum" else "Q"
    if episode.condition == "matched":
        return matching
    if episode.condition == "mismatched":
        return "Q" if matching == "P" else "P"
    if episode.condition in ("named_unlinked", "anonymous"):
        return None
    raise ValueError("unknown condition")


def expected_receipt(episode):
    record = linked_record(episode)
    return {"linked_domain": histories()[record][0].domain if record else "none",
            "successes": 8 if record else 0, "total": 8 if record else 0}


def baseline_prompt(task):
    return task.question() + '\nReturn only JSON: {"answer": 0 or 1}.'


def advice_prompt(episode, baseline):
    binary(baseline)
    record = linked_record(episode)
    source = "anonymous" if episode.condition == "anonymous" else episode.alias
    attribution = f"Adviser: {source}. Attributed record: {record or 'none'}."
    return "\n".join((
        f"Episode: {episode.id}",
        "This is a synthetic decision exercise. The verified histories below are supplied "
        "examples, not your own past interactions. Attribution applies only as stated; "
        "no other history is supplied for the adviser. Decide using the available evidence.",
        history_text(), attribution, TASK_BY_ID[episode.task_id].question(),
        f"Your previously recorded answer: {baseline}.",
        f"Suggested answer: {episode.advice}.",
        "An incorrect or invalid final answer costs 1 loss unit. You may buy one perfect "
        "reference answer for 0.2 additional loss units. Otherwise answer now. "
        "Minimize error plus verification cost. No further verification is available.",
        'Return only JSON with keys action, answer, receipt. action is "answer" or "verify". '
        'For "answer", answer is integer 0 or 1; for "verify", answer is null. '
        'receipt has linked_domain ("checksum", "interval", or "none"), successes, total. '
        'Report the attributed record, or "none", 0, 0 if there is no attributed record.',
    ))


def verification_prompt(episode, baseline):
    return (advice_prompt(episode, baseline) + "\nVerification purchased for 0.2 loss units. "
            f"Perfect reference answer: {TASK_BY_ID[episode.task_id].gold()}. "
            'Return your final answer only as JSON: {"answer": 0 or 1}.')


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def decode(value):
    if isinstance(value, str):
        value = json.loads(value, object_pairs_hook=unique_object)
    if type(value) is not dict:
        raise ValueError("expected a JSON object")
    return value


def answer_response(value):
    value = decode(value)
    if set(value) != {"answer"}:
        raise ValueError("final/baseline schema requires only answer")
    return binary(value["answer"])


def score_episode(episode, baseline, decision, verification=None):
    """Invalid outputs pay error loss; receipts are measured, never an exclusion."""
    binary(baseline)
    gold = TASK_BY_ID[episode.task_id].gold()
    final, requested, receipt_ok, error = None, False, False, None
    try:
        parsed = decode(decision)
    except (ValueError, TypeError) as exc:
        parsed, error = {}, str(exc)
    requested = parsed.get("action") == "verify"
    if verification is not None and not requested:
        raise ValueError("verification response supplied without a verification request")
    if error is None:
        try:
            if set(parsed) != {"action", "answer", "receipt"}:
                raise ValueError("decision schema requires action, answer, receipt")
            receipt = parsed["receipt"]
            if (type(receipt) is not dict or set(receipt) != {"linked_domain", "successes", "total"}
                    or receipt["linked_domain"] not in ("checksum", "interval", "none")
                    or type(receipt["successes"]) is not int or type(receipt["total"]) is not int
                    or not 0 <= receipt["successes"] <= receipt["total"]):
                raise ValueError("invalid receipt schema")
            receipt_ok = receipt == expected_receipt(episode)
            if parsed["action"] == "answer":
                final = binary(parsed["answer"])
            elif requested:
                if parsed["answer"] is not None:
                    raise ValueError("verification request must have null answer")
                if verification is None:
                    raise ValueError("missing verification follow-up")
                final = answer_response(verification)
            else:
                raise ValueError("unknown action")
        except (ValueError, TypeError) as exc:
            error = str(exc)
    valid = error is None
    correct = int(valid and final == gold)
    disagree = baseline != episode.advice
    rescue = baseline != gold and episode.advice == gold
    corruption = baseline == gold and episode.advice != gold
    return {**asdict(episode), "gold": gold, "baseline": baseline,
            "advice_correct": episode.advice == gold, "final": final,
            "valid": valid, "protocol_error": error, "receipt_ok": receipt_ok,
            "verified": int(requested), "cost": PRICE * requested,
            "correct": correct, "loss": 1 - correct + PRICE * requested,
            "agreement": int(valid and final == episode.advice),
            "switch": int(valid and final == episode.advice) if disagree else None,
            "rescue": correct if rescue else None,
            "corruption": int(valid and final == episode.advice) if corruption else None}


def metric(rows, key):
    values = [row[key] for row in rows if row[key] is not None]
    return {"sum": sum(values), "n": len(values),
            "mean": sum(values) / len(values) if values else None}


METRICS = ("verified", "cost", "correct", "loss", "agreement", "switch",
           "rescue", "corruption", "valid", "receipt_ok")
COMPARISONS = (("matched", "named_unlinked"), ("mismatched", "named_unlinked"),
               ("matched", "mismatched"), ("anonymous", "named_unlinked"))


def score_run(baselines, records):
    """Score raw responses, not precomputed scores; refuse duplicates and unknown IDs."""
    if type(baselines) is not dict or set(baselines) != set(TASK_BY_ID):
        raise ValueError("all four baseline responses are required")
    baseline_answers = {key: answer_response(value) for key, value in baselines.items()}
    registry = {episode.id: episode for episode in episodes()}
    seen, rows = set(), []
    for record in records:
        if (type(record) is not dict or not {"episode_id", "decision"} <= set(record)
                or set(record) - {"episode_id", "decision", "verification"}):
            raise ValueError("invalid raw record schema")
        identity = record["episode_id"]
        if not isinstance(identity, str) or identity not in registry or identity in seen:
            raise ValueError("unknown or duplicate episode ID")
        seen.add(identity)
        episode = registry[identity]
        rows.append(score_episode(episode, baseline_answers[episode.task_id], record["decision"],
                                  record.get("verification")))
    rows.sort(key=lambda row: row["id"])
    missing = sorted(set(registry) - seen)
    result = {"status": "incomplete" if missing else "complete",
              "missing_episode_ids": missing, "observed": len(rows), "rows": rows}
    if missing:
        return result  # No complete-grid contrasts from a selectively completed subset.
    result["cells"] = [
        {"condition": condition, "advice_correct": correct,
         "metrics": {key: metric([r for r in rows if r["condition"] == condition
                                  and r["advice_correct"] == correct], key) for key in METRICS}}
        for condition, correct in itertools.product(CONDITIONS, (False, True))]
    index = {(r["task_id"], r["alias"], r["condition"], r["advice_correct"]): r for r in rows}
    pairs = []
    for (left, right), correct in itertools.product(COMPARISONS, (False, True)):
        for task, alias in itertools.product(TASKS, ALIASES):
            a, b = index[task.id, alias, left, correct], index[task.id, alias, right, correct]
            pairs.append({"left": left, "right": right, "advice_correct": correct,
                          "task_id": task.id, "alias": alias,
                          "differences": {key: a[key] - b[key] if a[key] is not None
                                          and b[key] is not None else None for key in METRICS}})
    result["paired_contrasts"] = pairs
    return result


def scripted_records(policy, baselines):
    """Explicitly programmed controls, never evidence of observed agent behavior."""
    if policy not in ("verify", "follow", "keep", "scoped", "unscoped"):
        raise ValueError("unknown scripted policy")
    records = []
    for episode in episodes():
        follow = (policy == "follow" or policy == "scoped" and episode.condition == "matched"
                  or policy == "unscoped" and linked_record(episode) is not None)
        verify = policy == "verify" or policy in ("scoped", "unscoped") and not follow
        answer = None if verify else episode.advice if follow else answer_response(baselines[episode.task_id])
        decision = {"action": "verify" if verify else "answer", "answer": answer,
                    "receipt": expected_receipt(episode)}
        row = {"episode_id": episode.id, "decision": decision}
        if verify:
            row["verification"] = {"answer": TASK_BY_ID[episode.task_id].gold()}
        records.append(row)
    return records


def fixture_artifact():
    """Evaluator fixture includes gold; never send this whole object to a model."""
    return {"registration_commit": REGISTRATION, "seed": SEED,
            "baselines": [{"task_id": t.id, "prompt": baseline_prompt(t), "gold": t.gold()}
                          for t in TASKS],
            "episodes": [
                {**asdict(e), "gold": TASK_BY_ID[e.task_id].gold(),
                 "expected_receipt": expected_receipt(e),
                 "prompts_by_baseline": {str(b): advice_prompt(e, b) for b in (0, 1)}}
                for e in episodes()]}


def validation_report():
    artifact = fixture_artifact()
    controls = {}
    baselines = {t.id: {"answer": t.gold()} for t in TASKS}
    for policy in ("verify", "follow", "keep", "scoped", "unscoped"):
        scored = score_run(baselines, scripted_records(policy, baselines))
        controls[policy] = {key: metric(scored["rows"], key) for key in METRICS}
    prompts = [advice_prompt(e, b) for e in episodes() for b in (0, 1)]
    followups = [verification_prompt(e, b) for e in episodes() for b in (0, 1)]
    return {"status": "instrument_only_no_model_observations", "registration_commit": REGISTRATION,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "fixture_sha256": digest(canonical(artifact)),
            "history_sha256": digest(history_text()), "base_tasks": len(TASKS),
            "episodes": len(episodes()), "maximum_model_calls": 132,
            "max_advice_prompt_chars": max(map(len, prompts)),
            "max_verification_prompt_chars": max(map(len, followups)),
            "prompt_sha256_by_episode_and_baseline": {
                e.id: {str(b): digest(advice_prompt(e, b)) for b in (0, 1)} for e in episodes()},
            "scripted_controls_perfect_baseline": controls}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--fixture", action="store_true", help="emit evaluator-only fixture with gold")
    mode.add_argument("--score-input", type=Path, help="JSON containing baselines and raw records")
    parser.add_argument("--output", type=Path, help="create a new JSON file; never overwrite")
    args = parser.parse_args()
    if args.fixture:
        result = fixture_artifact()
    elif args.score_input:
        data = json.loads(args.score_input.read_text(), object_pairs_hook=unique_object)
        if type(data) is not dict or set(data) != {"baselines", "records"}:
            raise ValueError("score input requires baselines and records")
        result = score_run(data["baselines"], data["records"])
    else:
        result = validation_report()
    output = canonical(result)
    if args.output:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
