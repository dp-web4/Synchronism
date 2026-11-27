#!/usr/bin/env python3
"""Publisher 2026-11-09: how far under the bar is TEST-08, as an interval rather than a point?

The paper's count list prints "r² = 0.0001 against a registered > 0.20". PREDICTIONS.md 2026-11-08 (iv) says the
95% upper limit is r² ≈ 0.03. This recomputes that from the tracked per-galaxy table rather than relaying it.

Reads simulations/test08_per_galaxy_results.json (written by test08_sparc_environment_rar.py). Writes nothing.
Rebuilds the pre-declared primary statistic the same way that script does: per-galaxy RAR offset against
log10(1 + N_cyl[2 Mpc, 500 km/s]) with log D regressed out.

Expected before running: r = +0.012, N = 141, Fisher 95% upper limit on r² near 0.03, bootstrap within ~0.01 of that.
"""
import json
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
T = json.load(open(HERE / "test08_per_galaxy_results.json"))
names = sorted(T)
off = np.array([T[n]["offset_dex"] for n in names])
cyl = np.log10(1 + np.array([T[n]["cyl2"] for n in names], float))
logD = np.log10(np.array([T[n]["D"] for n in names], float))


def residual(y, x):
    A = np.stack([np.ones_like(x), x], axis=1)
    return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]


env = residual(cyl, logD)
r = float(np.corrcoef(off, env)[0, 1])
N = len(names)
print(f"N = {N}, r = {r:+.4f}, r^2 = {r**2:.5f}   (identity: registered run printed r=+0.012, N=141, r^2=0.0001)")

# Fisher z interval. The largest |r| in the interval sets the upper limit on r^2.
z, se = np.arctanh(r), 1 / np.sqrt(N - 3)
lo, hi = np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)
fisher_ul = max(lo**2, hi**2)
print(f"Fisher 95% r in [{lo:+.3f}, {hi:+.3f}]; upper limit on r^2 = {fisher_ul:.4f}")

# Bootstrap over galaxies, refitting the distance correction on each draw. Three seeds: a single seed's
# agreement with itself says nothing about whether the quantile is stable.
for seed in (1, 2, 3):
    rng = np.random.default_rng(seed)
    rs = np.empty(20000)
    for k in range(rs.size):
        i = rng.integers(0, N, N)
        rs[k] = np.corrcoef(off[i], residual(cyl[i], logD[i]))[0, 1]
    print(f"bootstrap seed {seed}: 95% r in [{np.quantile(rs, .025):+.3f}, {np.quantile(rs, .975):+.3f}]; "
          f"97.5% quantile of r^2 = {np.quantile(rs**2, .975):.4f}; 95% quantile of r^2 = {np.quantile(rs**2, .95):.4f}")

for label, bar in (("kill bar", 0.09), ("registered claim", 0.20)):
    print(f"{label} {bar}: {bar / fisher_ul:.1f}x the Fisher upper limit; {bar / r**2:.0f}x the point estimate")

print(f"\nper-galaxy offset sd = {off.std():.4f} dex (the variance r^2 is a share of)")
print("the registered sentence reads 'environment explains >20% of RAR scatter' and names no denominator;")
print("the pooled point scatter in the registered run is 0.133 dex, a different quantity from the per-galaxy offset sd.")
print(f"offset variance / pooled point variance = {off.var() / 0.133**2:.3f}")
