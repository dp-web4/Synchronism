"""Which mechanism moves the declared <C> off 1+1/(P-1): the degenerate-pair skip (analyze.py:75-76)
or unequal residual variances under Pearson normalisation? (CBP-Claude, 2026-10-08)
cov-ratio <C> = 1 - mean off-diagonal cov / mean variance is EXACTLY 1+1/(P-1) by the zero-sum identity."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from analyze import compatibility, load, arm_pack
P = 12
def report(tag, c):
    r = c - c.mean(1, keepdims=True); sd = r.std(0)
    skipped = sum(1 for a in range(P) for b in range(a+1, P) if sd[a] <= 1e-9 or sd[b] <= 1e-9)
    C = np.cov(r, rowvar=False, bias=True)
    covr = 1 - C[~np.eye(P, dtype=bool)].mean() / np.diag(C).mean()
    print(f"{tag:22s} skipped_pairs={skipped:2d}  declared<C>={compatibility(c):.4f}  "
          f"cov-ratio<C>={covr:.4f}  resid-sd CV={sd.std()/sd.mean():.2f}")
rng = np.random.default_rng(7); I = 48
for rho in [0.0, 0.3, 0.6, 0.9, 0.99]:
    z = rng.normal(size=(I, 1)); o = rng.normal(size=(I, P)); lat = np.sqrt(rho)*z + np.sqrt(1-rho)*o
    report(f"synthetic rho={rho}", (lat > np.quantile(lat, 0.45)).astype(float))
arms, _, _ = load(os.path.join(HERE, "answers_full.jsonl")); packs = arm_pack(arms)
for k in sorted(packs): report(f"real K={k}", packs[k][1])
