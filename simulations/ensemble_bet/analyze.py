#!/usr/bin/env python3
"""Analysis for the agent-ensemble compatibility bet.

Implements exactly the measurements declared in
explorations/2026-09-30-prereg-exec-agent-ensemble-kimi.md. Reads answers.jsonl,
writes analysis.json + prints the verdict-relevant numbers.
"""
import json
import math
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit

POOL = 12
NS = [1, 3, 5, 7, 9, 12]
BOOT = 200
SLOPE_BOOT = 2000
RNG = np.random.default_rng(20260930)


def parse_int(text: str):
    """The last line's integer, per the prompt's rule; last integer token fallback."""
    lines = [ln.strip() for ln in text.strip().splitlines() if ln.strip()]
    for ln in reversed(lines):
        try:
            return int(ln.strip(".,* "))
        except ValueError:
            for tok in reversed(ln.replace(",", " ").split()):
                try:
                    return int(tok.strip(".,*"))
                except ValueError:
                    continue
    return None


def load(path: str):
    """k -> item -> pool_i -> (persona, parsed_answer, correct)."""
    arms = defaultdict(dict)
    persona_correct = defaultdict(list)
    unparsed = 0
    for line in open(path):
        r = json.loads(line)
        got = parse_int(r["raw"])
        if got is None:
            unparsed += 1
        correct = int(got is not None and got == r["answer"])
        arms[r["k"]].setdefault(r["item"], {})[r["pool_i"]] = (r["persona"], got, correct)
        persona_correct[r["persona"]].append(correct)
    return arms, persona_correct, unparsed


def arm_pack(arms):
    """k -> (answers array (items, pool) object, correctness array float, true answers array)."""
    out = {}
    for k, per_item in arms.items():
        full = sorted((it, d) for it, d in per_item.items() if len(d) == POOL)
        ans = np.array([[d[pi][1] if d[pi][1] is not None else -10**9 for pi in range(POOL)]
                        for _, d in full], dtype=object)
        cor = np.array([[d[pi][2] for pi in range(POOL)] for _, d in full], dtype=float)
        truth = None  # correctness already encodes truth
        out[k] = (ans, cor)
    return out


def compatibility(cor: np.ndarray) -> float:
    """1 - mean residual pairwise Pearson corr, residuals against pool-wide item difficulty."""
    difficulty = cor.mean(axis=1, keepdims=True)
    resid = cor - difficulty
    n = resid.shape[1]
    cors = []
    for i in range(n):
        for j in range(i + 1, n):
            a, b = resid[:, i], resid[:, j]
            if a.std() < 1e-9 or b.std() < 1e-9:
                continue
            cors.append(float(np.corrcoef(a, b)[0, 1]))
    return 1.0 - float(np.mean(cors)) if cors else float("nan")


def vote_correct(answers_row: np.ndarray, correct_row: np.ndarray) -> int:
    """Majority vote over parsed answers; ties break to numerically smallest candidate
    (declared). Returns 1 iff the winning answer is a correct one."""
    vals, counts = np.unique(answers_row.astype(np.int64), return_counts=True)
    top = counts.max()
    winners = sorted(int(v) for v, c in zip(vals, counts) if c == top)
    voted = winners[0]
    return int(any(int(a) == voted and c == 1 for a, c in zip(answers_row, correct_row)))


def ensemble_accuracy(ans: np.ndarray, cor: np.ndarray, n: int) -> float:
    n_items = ans.shape[0]
    if n_items == 0:
        return float("nan")
    accs = []
    for _ in range(BOOT):
        items = RNG.integers(0, n_items, n_items)
        hits = 0
        for it in items:
            members = RNG.integers(0, ans.shape[1], n)
            hits += vote_correct(ans[it, members], cor[it, members])
        accs.append(hits / n_items)
    return float(np.mean(accs))


def slope_vs_lnN(ans: np.ndarray, cor: np.ndarray) -> tuple[float, float, float]:
    """Bootstrap distribution of slope(accuracy ~ ln N) over item-resamples and member-draws."""
    n_items = ans.shape[0]
    slopes = []
    for _ in range(SLOPE_BOOT):
        items = RNG.integers(0, n_items, n_items)
        pts = []
        for n in NS:
            hits = 0
            for it in items:
                members = RNG.integers(0, ans.shape[1], n)
                hits += vote_correct(ans[it, members], cor[it, members])
            pts.append(hits / n_items)
        slopes.append(float(np.polyfit(np.log(NS), pts, 1)[0]))
    lo, hi = np.percentile(slopes, [2.5, 97.5])
    return float(np.mean(slopes)), float(lo), float(hi)


def hill(c, base, amp, c50):
    return base + amp * np.power(c, 4) / (np.power(c, 4) + c50**4)


def main(path: str = "answers.jsonl") -> None:
    arms, persona_correct, unparsed = load(path)
    packs = arm_pack(arms)
    ks = sorted(packs)
    print(f"arms: {ks}; unparsed answers: {unparsed}")
    print("per-persona accuracy: " + ", ".join(
        f"{p}={np.mean(v):.3f}(n={len(v)})" for p, v in sorted(persona_correct.items())))

    rows = {}
    for k in ks:
        ans, cor = packs[k]
        c_val = compatibility(cor)
        accs = {n: ensemble_accuracy(ans, cor, n) for n in NS}
        rows[k] = {"C": c_val, "q": float(cor.mean()), "n_items": ans.shape[0],
                   "acc": {str(n): a for n, a in accs.items()}}
        print(f"K={k:2d} q={cor.mean():.3f} C={c_val:+.3f} " +
              " ".join(f"N{n}:{a:.3f}" for n, a in sorted(accs.items())))

    # Kill 1: slope of accuracy on ln N at the lowest-C arm
    k_low = min(rows, key=lambda k: rows[k]["C"])
    mean_s, lo, hi = slope_vs_lnN(*packs[k_low])
    print(f"\nKILL-1 arm K={k_low}: slope acc~lnN = {mean_s:+.4f} 95% CI [{lo:+.4f}, {hi:+.4f}]"
          f" -> {'REFUTED (count compensates)' if lo > 0 else 'survives'}")

    # Kill 2: sharpness at N=9 — Hill vs linear on measured C
    xs = np.array([rows[k]["C"] for k in ks])
    ys = np.array([rows[k]["acc"]["9"] for k in ks])
    lin = np.polyfit(xs, ys, 1)
    resid_lin = ys - np.polyval(lin, xs)
    aic_lin = len(ys) * math.log(max(float(np.mean(resid_lin**2)), 1e-12)) + 2 * 2
    try:
        popt, _ = curve_fit(hill, xs, ys, p0=[ys.min(), max(ys.max() - ys.min(), 1e-3),
                                              max(np.median(xs), 1e-3)], maxfev=20000)
        resid_h = ys - hill(xs, *popt)
        aic_hill = len(ys) * math.log(max(float(np.mean(resid_h**2)), 1e-12)) + 2 * 3
    except Exception as e:
        aic_hill, popt = float("inf"), None
        print("hill fit failed:", e)
    print(f"KILL-2: AIC linear {aic_lin:.2f} vs Hill {aic_hill:.2f} "
          f"(ΔAIC Hill-lin = {aic_hill - aic_lin:+.2f}) -> "
          f"{'survives (Hill better)' if aic_hill <= aic_lin - 2 else 'REFUTED (no sharp knee)'}")

    # Kill 3: does C add beyond {q, N}? least squares on per-(k,N) accuracy
    X, Xn, Y = [], [], []
    for k in ks:
        for n in NS:
            Y.append(rows[k]["acc"][str(n)])
            X.append([1.0, math.log(n), rows[k]["q"], rows[k]["C"]])
            Xn.append([1.0, math.log(n), rows[k]["q"]])
    X, Xn, Y = np.array(X), np.array(Xn), np.array(Y)
    b_full, *_ = np.linalg.lstsq(X, Y, rcond=None)
    b_null, *_ = np.linalg.lstsq(Xn, Y, rcond=None)
    rss_full = float(np.mean((Y - X @ b_full) ** 2))
    rss_null = float(np.mean((Y - Xn @ b_null) ** 2))
    aic_full = len(Y) * math.log(max(rss_full, 1e-12)) + 2 * X.shape[1]
    aic_null = len(Y) * math.log(max(rss_null, 1e-12)) + 2 * Xn.shape[1]
    print(f"KILL-3: AIC with C {aic_full:.2f} vs without {aic_null:.2f} "
          f"(ΔAIC = {aic_full - aic_null:+.2f}; C coef {b_full[-1]:+.3f}) -> "
          f"{'survives (C adds)' if aic_full <= aic_null - 2 else 'REFUTED (C is a relabel)'}")

    json.dump({"per_arm": rows,
               "kill1": {"arm": k_low, "slope_mean": mean_s, "ci": [lo, hi]},
               "kill2": {"aic_linear": aic_lin, "aic_hill": aic_hill,
                          "hill_params": None if popt is None else [float(x) for x in popt]},
               "kill3": {"aic_with_C": aic_full, "aic_without_C": aic_null,
                          "C_coef": float(b_full[-1])},
               "per_persona_q": {p: float(np.mean(v)) for p, v in persona_correct.items()},
               "unparsed": unparsed},
              open("analysis.json", "w"), indent=1)
    print("\nwrote analysis.json")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "answers.jsonl")
