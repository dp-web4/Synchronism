"""Synthetic positive control for the interacting-ensemble compatibility design (CBP-Claude, 2026-10-16).

Question (pre-stated): with pre-interaction answers held LITERALLY IDENTICAL (same marginals, same votes),
can a pairwise property (who can verify whom) change post-interaction accuracy? And can a pairwise proxy
measured on separate calibration items recover it, where the 09-30 residual-correlation proxy cannot?

Model: I items, K categories, 4 options. Agents answer once (fixed matrix A). One critique-then-revise round:
agent a adopts peer b's differing answer with prob V[a,b,k] if that answer is correct, F if wrong
(false acceptance); among multiple accepted proposals it takes the most-proposed. Final = plurality.
Verification structure V varies across arms; A does not.
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ensemble_bet"))
from analyze import compatibility  # the 09-30 declared proxy, unchanged

rng = np.random.default_rng(20261016)
N, K, OPT = 9, 3, 4
I_TEST, I_CAL = 300, 120
F = 0.08  # false-acceptance rate of a wrong peer answer

# per-agent per-category skill: three specialist groups, each strong on one category (complementary coverage)
skill = np.full((N, K), 0.35)
for a in range(N): skill[a, a % K] = 0.80

def answer_matrix(I, cats):
    truth = rng.integers(0, OPT, I)
    A = np.empty((N, I), int)
    for a in range(N):
        for i in range(I):
            if rng.random() < skill[a, cats[i]]: A[a, i] = truth[i]
            else: A[a, i] = rng.choice([o for o in range(OPT) if o != truth[i]])
    return truth, A

def plurality(col):
    v, c = np.unique(col, return_counts=True); top = c.max()
    return min(int(x) for x, y in zip(v, c) if y == top)

def interact(A, truth, cats, V, rr):
    B = A.copy()
    for i in range(A.shape[1]):
        k = cats[i]
        for a in range(N):
            acc = []
            for b in range(N):
                if b == a or A[b, i] == A[a, i]: continue
                p = V[a, b, k] if A[b, i] == truth[i] else F
                if rr.random() < p: acc.append(A[b, i])
            if acc: B[a, i] = plurality(np.array(acc))
    return B

def accuracy(M, truth): return float(np.mean([plurality(M[:, i]) == truth[i] for i in range(M.shape[1])]))

# verification arms: same mean verification budget, different STRUCTURE
def V_arm(kind, level):
    V = np.zeros((N, N, K))
    for a in range(N):
        for b in range(N):
            for k in range(K):
                if kind == "uniform": V[a, b, k] = level
                elif kind == "expert-trusting":  # a can verify b where b is the category expert
                    V[a, b, k] = min(1.0, level * K) if b % K == k else 0.0
                elif kind == "anti":  # a only 'verifies' peers on categories they are weak in
                    V[a, b, k] = 0.0 if b % K == k else min(1.0, level * K / (K - 1))
    return V

def calib_proxy(V, rr):
    """Measured on SEPARATE calibration items: how often a accepts b's CORRECT differing answer minus how often
    it accepts b's WRONG one, weighted by how often b is right where a is wrong. Needs only accept/reject logs."""
    cats = rr.integers(0, K, I_CAL); truth, A = answer_matrix(I_CAL, cats)
    hits = fa = opp_r = opp_w = 0
    for i in range(I_CAL):
        k = cats[i]
        for a in range(N):
            for b in range(N):
                if b == a or A[b, i] == A[a, i]: continue
                if A[b, i] == truth[i]:
                    opp_r += 1; hits += rr.random() < V[a, b, k]
                else:
                    opp_w += 1; fa += rr.random() < F
    return hits / max(opp_r, 1) - fa / max(opp_w, 1)

cats = rng.integers(0, K, I_TEST)
truth, A = answer_matrix(I_TEST, cats)
base_vote = accuracy(A, truth)
decl = compatibility((A == truth).T.astype(float))
print(f"pre-interaction (IDENTICAL for every arm): q={np.mean(A == truth):.3f} vote_acc={base_vote:.3f} "
      f"declared09-30<C>={decl:.4f} (identity 1+1/(N-1)={1 + 1/(N-1):.4f})")
rows = []
for kind in ["uniform", "expert-trusting", "anti"]:
    for level in [0.05, 0.1, 0.2, 0.3, 0.5]:
        V = V_arm(kind, level); rr = np.random.default_rng(7)
        post = np.mean([accuracy(interact(A, truth, cats, V, rr), truth) for _ in range(5)])
        prox = calib_proxy(V, np.random.default_rng(11))
        rows.append(dict(kind=kind, level=level, mean_V=float(V.mean()), proxy=prox, post=float(post), gain=float(post - base_vote)))
        print(f"{kind:16s} level={level:.2f} meanV={V.mean():.3f} calib_proxy={prox:+.3f} post_acc={post:.3f} gain={post - base_vote:+.3f}")
mv = np.array([r["mean_V"] for r in rows]); px = np.array([r["proxy"] for r in rows]); g = np.array([r["gain"] for r in rows])
print(f"\ncorr(gain, mean_V)={np.corrcoef(g, mv)[0,1]:+.3f}   corr(gain, calib_proxy)={np.corrcoef(g, px)[0,1]:+.3f}")
json.dump(dict(base_vote=base_vote, declared_proxy=decl, rows=rows), open(os.path.join(HERE, "positive_control.json"), "w"), indent=1)
