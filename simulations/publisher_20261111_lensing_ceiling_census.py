#!/usr/bin/env python3
"""Publisher arithmetic, 2026-11-11 (cbp-claude, interactive maintenance pass).

Whitepaper 5.15, "The ceiling also dies independently" (Brouwer et al. 2021 lensing RAR).
The paragraph prints deficits 6x / 13x / 28x at g_bar = 1e-14 and says a factor-2 baryon
correction "moves nu by sqrt(2)". Two things are recomputed here from the paragraph's own inputs:

  * the boost a ceiling must supply is g_obs/g_bar. Lensing fixes g_obs, so raising the
    baryons by f lowers the requirement by f (not sqrt(f));
  * Brouwer et al. Sec. 1 + footnote 3: the outer edge is g_bar ~ 1e-15 counting stars and
    cold gas only; hot gas, dust and missing baryons could raise it to ~1e-14 (f ~ 10).

Assumption, stated: g_obs is taken on MOND's simple-nu curve at the fiducial g_bar (the
paper's "agrees well with the MG predictions"); no Brouwer data point is read here.
Inputs from the whitepaper: eps0 = 0.119 (M/1e10)^0.325, B_max = 1/eps0, a0 = 1.2e-10.
No file is read or written.
"""
import math
a0 = 1.2e-10
def nu(g):                      # simple interpolating function
    y = g / a0
    return 0.5 + math.sqrt(0.25 + 1.0 / y)
ceil = {m: 1.0 / (0.119 * (10 ** (m - 10)) ** 0.325) for m in (9, 10, 11)}
print("ceilings B_max = 1/eps0:", ", ".join(f"1e{m}: {c:.1f}" for m, c in ceil.items()))
print("framework cap 1/Omega_m = %.3f; most permissive cap quoted by the ledger = 6.39" % (1 / 0.315))
for g in (1e-14, 1e-15):
    for f in (1, 2, 5, 10):
        need = nu(g) / f
        row = "  ".join(f"1e{m}: {'pass' if need <= c else f'{need / c:4.1f}x short'}" for m, c in ceil.items())
        print(f"fiducial g_bar = {g:.0e}, baryons x{f:<2}: required boost = {need:6.1f} | {row} | 6.39: "
              f"{'pass' if need <= 6.39 else f'{need / 6.39:4.1f}x short'}")
print("check on the printed sentence: nu(2e-14)/nu(1e-14) = %.3f (= 1/sqrt 2: that is MOND's own nu moving,"
      " not the required boost, which moves by 1/2)" % (nu(2e-14) / nu(1e-14)))
