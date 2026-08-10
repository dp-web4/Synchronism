#!/usr/bin/env python3
"""
Publisher 2026-08-10 — audit of the *corpus* rationale for the chemistry sector's gamma ~ 1.

Context
-------
On 2026-08-09 the synchronism-site maintainer lane established that the site's stated
rationale for gamma ~ 1 --

    "At gamma ~ 1, the coherence function has maximum curvature."

-- is mathematically false, and concluded from a literal grep of "maximum curvature" that the
claim is *site-originated*: "Nothing in Research/ needs correcting for this."  That proposal
explicitly deferred one check: "Worth confirming against the chemistry session archive
(sessions 134-2660), which this grep did not cover."

This script runs that deferred check's consequence.  The phrase is indeed absent from the
chemistry archive -- but the *claim* is present, in different words, in two load-bearing
corpus documents:

    Research/Chemistry/MASTER_PREDICTIONS.md:574  (P15.4)
        "Mechanism: Correlation length diverges at gamma = 1 (N_corr = 4)"

    Research/Chemistry/Session20_Complexity_Emergence.md:33-43
        C_eff  ~  (2/gamma) * (gamma/2) * exp(-(gamma-1)^2 / sigma^2)
        "Peaks at gamma = 1.0 where N_corr = 4."

So there is a corpus-side rationale after all.  This script tests it.

Checks
------
A. Reproduce the site lane's three results (independent re-derivation, symbolic + numeric).
B. The Session-20 "balance of order and disorder" factor: is it doing any work?
C. Can the published Session-20 table be reproduced by the Gaussian factor ALONE?
D. P15.4's "correlation length diverges at gamma = 1" against the framework's own
   gamma = 2/sqrt(N_corr).

Data: none required -- these are claims about closed-form expressions.
"""

import numpy as np
import sympy as sp

print("=" * 78)
print("PUBLISHER 2026-08-10 -- gamma ~ 1 mechanism audit (corpus-side)")
print("=" * 78)

# ---------------------------------------------------------------------------
# A. Re-derive the site lane's three claims independently (sympy, not numeric diff).
#    My own standing lesson: a check that contradicts a proof is the suspect --
#    so use a CAS, and state which way the check cut.
# ---------------------------------------------------------------------------
print("\n[A] Independent re-derivation of the 2026-08-09 site-lane result")
print("-" * 78)

x, g = sp.symbols("x gamma", positive=True)
C = sp.tanh(g * sp.log(1 + x))

dC_dx = sp.simplify(sp.diff(C, x))
slope_at_0 = sp.simplify(sp.limit(dC_dx, x, 0))
print(f"  dC/dx at x->0            : {slope_at_0}")
print(f"  monotone increasing in g : {sp.simplify(sp.diff(slope_at_0, g))} > 0  -> no interior max")

d2C_dx2 = sp.simplify(sp.diff(C, x, 2))
# concavity: evaluate on a grid at gamma = 1 and check sign
neg = []
for gv in (0.25, 0.5, 1.0, 2.0, 4.0):
    vals = [float(d2C_dx2.subs({g: gv, x: xv})) for xv in (0.01, 0.1, 0.5, 2.0, 10.0, 100.0)]
    neg.append(all(v < 0 for v in vals))
print(f"  d2C/dx2 < 0 on all tested (gamma, x) : {all(neg)}   -> C concave, no inflection")

# max of dC/dlnx over x, as a function of gamma
def max_dC_dlnx(gv, n=200001):
    lx = np.linspace(-12, 12, n)
    xv = np.exp(lx)
    t = 1.0 + xv
    Cv = np.tanh(gv * np.log(t))
    return np.max(np.gradient(Cv, lx))

print("\n  max_x dC/dln x  (the log-density reading):")
prev = -np.inf
mono = True
for gv in (0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 50.0):
    m = max_dC_dlnx(gv)
    mono &= (m > prev)
    prev = m
    star = "  <- gamma = 1" if gv == 1.0 else ""
    print(f"    gamma={gv:<6g} max={m:.4f}{star}")
print(f"  strictly increasing in gamma: {mono}  (saturating, ceiling 0.5*tanh-limit)")
print("  VERDICT A: site lane CONFIRMED -- no feature of C singles out gamma ~ 1.")
print("             (This check cut TOWARD the site lane's finding, not against it.)")

# ---------------------------------------------------------------------------
# B. The Session-20 'balance' factor.
# ---------------------------------------------------------------------------
print("\n[B] Session 20's stated mechanism: 'Balance between order and disorder'")
print("-" * 78)
order_term = 2 / g
disorder_term = g / 2
balance = sp.simplify(order_term * disorder_term)
print(f"  order term      (2/gamma) : {order_term}")
print(f"  disorder term   (gamma/2) : {disorder_term}")
print(f"  product                   : {balance}")
print(f"  is identically 1          : {sp.simplify(balance - 1) == 0}")
print("  => the 'balance of order and disorder' contributes EXACTLY NOTHING.")
print("     The two terms are reciprocals; their product is 1 for every gamma.")

sigma = sp.symbols("sigma", positive=True)
C_eff_full = balance * sp.exp(-((g - 1) ** 2) / sigma ** 2)
C_eff_gauss = sp.exp(-((g - 1) ** 2) / sigma ** 2)
print(f"  C_eff reduces to          : {sp.simplify(C_eff_full)}")
print(f"  identical to Gaussian only: {sp.simplify(C_eff_full - C_eff_gauss) == 0}")
print("  => the peak at gamma = 1 is the hand-written '1' in (gamma - 1)^2.")
print("     The ansatz ASSUMES its conclusion. This is circular, not derived.")

# ---------------------------------------------------------------------------
# C. Reproduce the published Session-20 table from the Gaussian factor alone.
# ---------------------------------------------------------------------------
print("\n[C] Published Session-20 table vs Gaussian-only prediction")
print("-" * 78)
published = [(0.5, 16, 0.73), (0.85, 5.5, 0.97), (1.0, 4, 1.00),
             (1.5, 1.8, 0.73), (2.0, 1, 0.29)]

# Calibrate sigma on the single gamma = 0.5 row, then PREDICT the rest.
s2 = -0.25 / np.log(0.73)
print(f"  sigma^2 calibrated on the gamma=0.5 row alone: {s2:.4f}  (sigma = {np.sqrt(s2):.4f})")
print(f"  {'gamma':>6} {'N_corr(pub)':>12} {'C_eff(pub)':>11} {'Gaussian only':>14} {'|diff|':>8}")
worst = 0.0
for gv, nc, pub in published:
    pred = np.exp(-((gv - 1) ** 2) / s2)
    worst = max(worst, abs(pred - pub))
    print(f"  {gv:>6g} {nc:>12g} {pub:>11.2f} {pred:>14.4f} {abs(pred-pub):>8.4f}")
print(f"  worst absolute deviation over all 5 rows: {worst:.4f}")
print("  => every published row is reproduced by exp(-(gamma-1)^2/sigma^2) with ONE free")
print("     parameter. The order/disorder physics contributes zero explanatory content.")

# The symmetry tell, stated independently of the algebra:
print(f"\n  Symmetry tell: published C_eff(0.5) = 0.73 and C_eff(1.5) = 0.73 -- exactly equal.")
print("  A genuine (2/gamma)*(gamma/2) modulation could only preserve that symmetry by")
print("  being identically 1. The table itself testifies that the prefactor does nothing.")

# ---------------------------------------------------------------------------
# D. P15.4: "Correlation length diverges at gamma = 1 (N_corr = 4)"
# ---------------------------------------------------------------------------
print("\n[D] P15.4 mechanism vs the framework's own gamma = 2/sqrt(N_corr)")
print("-" * 78)
print("  Framework identity (CLAUDE.md, GAMMA_UNIFICATION.md):  gamma = 2/sqrt(N_corr)")
print("  Inverted:                                              N_corr = 4/gamma^2")
print(f"  {'gamma':>8} {'N_corr':>14}   regime label used by the corpus itself")
labels = {2.0: "Uncorrelated / Disorder", 1.0: "'Critical (peak)'",
          0.5: "Ordered (protein-like)", 0.1: "-", 0.01: "-"}
for gv in (2.0, 1.0, 0.5, 0.1, 0.01):
    print(f"  {gv:>8g} {4.0/gv**2:>14.1f}   {labels[gv]}")
print("\n  N_corr is finite and SMALL (= 4) at gamma = 1.")
print("  N_corr -> infinity as gamma -> 0, i.e. the corpus's 'Fully correlated / Rigid order'")
print("  end of its own axis -- NOT at gamma = 1.")
print("  Correlation length is monotone in N_corr under any reading, so under the framework's")
print("  OWN formula the divergence sits at gamma -> 0.")
print("  => P15.4's stated mechanism is not merely unsupported; it points the WRONG WAY.")
print("\n  Internal inconsistency, same corpus:")
print("    MASTER_PREDICTIONS.md:574  'Correlation length diverges at gamma = 1'")
print("    Session20 Part 4.1         'Correlation length comparable to system size'")
print("  'Diverges' and 'comparable to system size' are different claims, and the framework's")
print("  own N_corr = 4 supports neither.")

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
print("""
The site lane's 2026-08-09 finding is CONFIRMED and its provenance verdict is CORRECTED.

  * "Nothing in Research/ needs correcting" is false. The corpus carries its own gamma ~ 1
    rationale in two load-bearing documents; the literal grep for "maximum curvature" could
    not see it because the corpus states it in different words.

  * The corpus rationale is WEAKER than the one the site retracted:
      - Session 20's derivation is CIRCULAR -- its 'order x disorder balance' factor is
        identically 1, so the peak at gamma = 1 is inserted by hand in exp(-(gamma-1)^2).
      - P15.4's 'correlation length diverges at gamma = 1' is WRONG-SIGNED under the
        framework's own gamma = 2/sqrt(N_corr): divergence is at gamma -> 0.

  * This is a DEMOTION, not a refutation, and it should NOT move the refutation count.
    Voiding a rationale removes support for a claim; it does not refute the claim. The
    gamma ~ 1 clustering survives as an empirical regularity with no surviving mechanism
    on either the site side or the corpus side.
""")
