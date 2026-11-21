#!/usr/bin/env python3
"""Publisher control, 2026-11-11 (cbp-claude, interactive maintenance pass).

Question: the whitepaper's "Superconductivity | Tc ∝ exp(-γ/λ) | r = 0.948 | #62" row
(chemistry.md top-validated table; Appendix B.5). Where does 0.948 come from, and is the
γ in it independent of Tc?

Source: simulations/chemistry/superconductivity_coherence.py (Chemistry Session #62,
commit e4c903f5, 2026-01-17). This script executes that file's PARTS 1-4 unchanged (text up
to the "PART 5" banner; no file is written by that part), then recomputes the PART 5
correlation and three controls. Read-only.

Usage: python3 publisher_20261111_tc_gamma_control.py <path to superconductivity_coherence.py>
"""
import sys, io, contextlib
import numpy as np
from scipy import stats

src = open(sys.argv[1], encoding="utf-8").read()
cut = src.index('print("PART 5: Tc vs γ CORRELATION")')
cut = src.rfind("\n# ====", 0, cut)
ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(src[:cut], "superconductivity_coherence[:PART5]", "exec"), ns)

Tcs, Deltas, cats = ns["Tcs"], ns["Deltas"], ns["categories"]
gammas = ns["gammas"]
kB, bcs = ns["k_B_meV"], ns["BCS_ratio_theory"]
n = len(Tcs)
ratio = 2 * Deltas / (kB * Tcs)

def pr(label, x, y):
    r, p = stats.pearsonr(x, y)
    rs, ps = stats.spearmanr(x, y)
    print(f"{label:<58} Pearson r = {r:+.3f} (p = {p:.1e})   Spearman = {rs:+.3f}")
    return r

print(f"N = {n} superconductors; BCS_ratio_theory = {bcs}; categories: {sorted(set(cats))}")
print()
print("IDENTITY (the source's PART 5 number)")
r0 = pr("Tc vs 1/gamma, gamma as the source estimates it", 1 / gammas, Tcs)
print()
print("WHAT gamma IS: gamma = 2 / (ratio/BCS) x category factor, ratio = 2*Delta/(k_B*Tc)")
g_nocat = 2.0 / (ratio / bcs)
fac = gammas / g_nocat
for c in sorted(set(cats)):
    m = np.array([x == c for x in cats])
    print(f"   {c:<14} n = {m.sum():>2}  factor on gamma = {np.unique(np.round(fac[m],3))}  mean Tc = {Tcs[m].mean():7.1f} K")
print()
print("CONTROLS")
r1 = pr("C1  Tc vs 1/gamma, category factors removed (= ratio/2BCS)", 1 / g_nocat, Tcs)
r2 = pr("C2  Tc vs category factor alone (1/factor)", 1 / fac, Tcs)
rng = np.random.default_rng(20261111)
rs = []
for _ in range(20000):
    perm = rng.permutation(n)
    # shuffle the gap ratios among materials, keep each material's category factor
    rs.append(stats.pearsonr(1 / (2.0 / (ratio[perm] / bcs) * fac), Tcs)[0])
rs = np.array(rs)
print(f"C3  gap ratios shuffled among materials, category factors kept (20,000 draws, seed 20261111):")
print(f"      median r = {np.median(rs):+.3f}, 2.5-97.5% = [{np.percentile(rs,2.5):+.3f}, {np.percentile(rs,97.5):+.3f}], "
      f"fraction >= source r: {np.mean(rs >= r0):.4f}")
print()
print("WITHIN CATEGORY (factor constant, so this is Tc vs gap ratio)")
for c in sorted(set(cats)):
    m = np.array([x == c for x in cats])
    if m.sum() >= 4:
        pr(f"   {c} (n = {m.sum()})", 1 / gammas[m], Tcs[m])
    else:
        print(f"   {c} (n = {m.sum()}): too few")
