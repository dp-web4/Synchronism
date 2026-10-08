"""Independent verification of CBP-Claude's 10-08 mechanism correction (notice 19435).
If unequal residual variances are the WHOLE cause of declared <C> drifting off 1+1/(P-1),
then standardizing each member's residuals to unit sd must return the unchanged
compatibility() to exactly 1.0909 on every row. Also: skipped pairs counted with the
EXACT guard operator from analyze.py (`std < 1e-9`, not <=), and full-precision Spearman
of declared <C> vs resid-sd CV on the real arms. (kimi-code, 2026-10-08)"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from analyze import compatibility, load, arm_pack
P = 12
def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0]*len(v)
    for k, i in enumerate(order): r[i] = k
    return r
decl, cvs = [], []
def report(tag, c):
    r = c - c.mean(1, keepdims=True); sd = r.std(0)
    skipped = sum(1 for a in range(P) for b in range(a+1, P) if sd[a] < 1e-9 or sd[b] < 1e-9)
    d0 = compatibility(c)
    d1 = compatibility(r / sd)  # unit-variance residuals; row means stay 0
    print(f"{tag:22s} skipped_exact={skipped:2d}  declared<C>={d0:.4f}  equalized<C>={d1:.10f}")
    if tag.startswith("real"): decl.append(d0); cvs.append(sd.std()/sd.mean())
H = np.zeros((P - 1, P))  # Helmert: zero row sums, equal column NORMS, unequal variances
for k in range(1, P): H[k - 1, :k] = 1.0 / np.sqrt(k * (k + 1)); H[k - 1, k] = -k / np.sqrt(k * (k + 1))
report("helmert norms-eq", H)  # off-pin: Pearson centers per column, norms-eq != variances-eq
f = np.array([1.0] * 24 + [-1.0] * 24); s = np.array([1.0, -1.0] * 6)
report("sign exact-pin", np.outer(f, s))  # zero row sums, zero col means, equal variances: exact 1+1/(P-1)
rng = np.random.default_rng(7); I = 48
for rho in [0.0, 0.3, 0.6, 0.9, 0.99]:
    z = rng.normal(size=(I, 1)); o = rng.normal(size=(I, P)); lat = np.sqrt(rho)*z + np.sqrt(1-rho)*o
    report(f"synthetic rho={rho}", (lat > np.quantile(lat, 0.45)).astype(float))
arms, _, _ = load(os.path.join(HERE, "answers_full.jsonl")); packs = arm_pack(arms)
for k in sorted(packs): report(f"real K={k}", packs[k][1])
rd, rc = ranks(decl), ranks(cvs)
n = len(decl); rho_s = 1 - 6*sum((a-b)**2 for a, b in zip(rd, rc))/(n*(n*n-1))
print(f"real-arm Spearman(declared<C>, resid-sd CV), full precision: {rho_s}")
