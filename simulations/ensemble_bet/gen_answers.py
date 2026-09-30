#!/usr/bin/env python3
"""Answer-pool generator for the agent-ensemble compatibility bet.

For each item and each compatibility arm (persona-pool size K), generate POOL
independent answers; pool member i uses persona roster[i % K]. Same model,
same sampling for every arm — per-agent capability q is held fixed by
construction and MEASURED per persona (reported, not assumed).

Sequential is the default for fleet citizenship; inside a GPU courtesy window
(beats resting, no contention) --workers 2-3 halves the window the run holds.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import sys
import time
import urllib.request

OLLAMA = "http://127.0.0.1:11434/api/generate"

ROSTER = [
    ("plain", "Solve the problem."),
    ("hurried", "You are in a hurry. Read the problem once, compute immediately, answer at once."),
    ("careful", "You are meticulous. Name each quantity and its role before computing."),
    ("verbal", "You briefly restate what the question actually asks, then compute."),
    ("skeptic", "You distrust your first reading. Re-read the final sentence to see exactly what is asked, then compute."),
    ("arithfirst", "Go straight to the arithmetic; do not re-read or second-guess."),
    ("teacher", "Answer as if checking a student's work: spot the classic mistake in this kind of problem, then compute."),
    ("distractaware", "Some numbers in the problem may be irrelevant. Decide which numbers matter, then compute."),
    ("estimator", "Estimate the answer's magnitude first, then compute exactly and check against your estimate."),
    ("backwards", "Work backwards: ask what the final quantity must equal, then assemble it from the given numbers."),
    ("minimalist", "Use the fewest operations possible; do not engage with numbers you do not need."),
    ("accountant", "Keep a running ledger: after each operation, note what quantity you now hold."),
]


def ask(model: str, persona: str, item_text: str, seed: int, num_predict: int = 48) -> dict:
    prompt = (f"{persona}\n\n{item_text}\n\nRules: the last line of your answer must be only the "
              f"final integer, no units, no words. Keep any reasoning to one short line.")
    body = {"model": model, "prompt": prompt, "stream": False, "think": False,
            "options": {"temperature": 0.7, "top_p": 0.95, "num_predict": num_predict, "seed": seed}}
    req = urllib.request.Request(OLLAMA, data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=180) as r:
        d = json.loads(r.read().decode())
    return {"text": d.get("response", ""), "secs": round(time.time() - t0, 2),
            "done": d.get("done_reason")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--items", default="items.json")
    ap.add_argument("--model", default="qwen3.5:4b")
    ap.add_argument("--ks", default="1,2,3,4,6,8,12")
    ap.add_argument("--pool", type=int, default=12)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="cap items (pilot)")
    ap.add_argument("--out", default="answers.jsonl")
    a = ap.parse_args()

    items = json.load(open(a.items))
    if a.limit:
        items = items[: a.limit]
    ks = [int(x) for x in a.ks.split(",")]
    assert a.pool <= len(ROSTER), "pool members must have distinct personas available"

    jobs = [(item, k, i) for item in items for k in ks for i in range(a.pool)]
    n_gen = len(jobs)
    t0 = time.time()
    done = [0]
    lock_out = open(a.out, "a")

    def run(job):
        item, k, i = job
        pname, ptext = ROSTER[i % k]
        r = ask(a.model, ptext, item["text"], seed=hash((item["id"], k, i)) % (2**31))
        return {"item": item["id"], "template": item["template"], "answer": item["answer"],
                "k": k, "pool_i": i, "persona": pname, "raw": r["text"], "secs": r["secs"],
                "done": r["done"]}

    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for rec in ex.map(run, jobs):
            lock_out.write(json.dumps(rec) + "\n")
            lock_out.flush()
            done[0] += 1
            if done[0] % 100 == 0:
                rate = done[0] / (time.time() - t0)
                eta = (n_gen - done[0]) / rate / 60
                print(f"{done[0]}/{n_gen} ({rate:.2f}/s, ETA {eta:.0f} min)", file=sys.stderr)
    lock_out.close()


if __name__ == "__main__":
    main()
