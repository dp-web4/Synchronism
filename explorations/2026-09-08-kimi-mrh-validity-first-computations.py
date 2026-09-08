#!/usr/bin/env python3
"""MRH-validity program, first computations (kimi-code, 2026-09-08).

Three measurements for the exploration chartered in
2026-09-08-kimi-mrh-validity-charter.md:

  1. THE GAUSSIAN PRICE OF A HORIZON. A linear system (X included, Y excluded)
     with a coupling knob b. The observer who excludes Y pays a prediction-error
     inflation of exactly exp(2 * I(Y ; X' | X)) — an identity for Gaussians,
     not a bound. Blanket closure (b=0) makes exclusion free: the Ptolemy case.

  2. THE DEBYE BOUNDARY ZONE. Exact Debye C_V vs the two MRH equations (low-T
     T^3 law, high-T Dulong-Petit 3Nk). Where each effective equation's error
     crosses 1% — and the shape of the loss (exponential mode-truncation vs
     continuum-approximation), which differs per horizon.

  3. N_corr AS THE CORRELATION STRUCTURE'S SHADOW. The relative fluctuation of
     the mean of N unit variables with pairwise correlation rho is
     sqrt(rho + (1-rho)/N), not 1/sqrt(N): the effective count N_corr collapses
     to 1/rho. The recurrence of gamma = 2/sqrt(N_corr) is the CLT fixed point —
     and deviations of the fluctuation exponent from 1/2 are where non-Gaussian
     (i.e., potentially new) physics would live.

Stdlib + numpy only; deterministic; prints the numbers the results doc quotes.
"""
import numpy as np

print("=" * 78)
print("1. THE GAUSSIAN PRICE OF A HORIZON  (error inflation = e^{2*CMI})")
print("=" * 78)
print("X' = a*X + b*Y + e_x ;  Y' = d*Y + e_y   (Y autonomous, Var e = 1)")
print("Observer's MRH includes X only. Sweep the exclusion leak b.\n")
print(f"{'b':>6} {'CMI (nats)':>12} {'Var(X\'|X)':>12} {'Var(X\'|X,Y)':>13} "
      f"{'inflation':>10} {'e^2CMI':>10} {'match':>7}")
a, d = 0.5, 0.6
for b in (0.0, 0.1, 0.3, 0.7, 1.5):
    # stationary moments: vy = 1/(1-d^2); cxy = b*d*vy/(1-a*d); vx = (a*b*cxy + b^2*vy + 1)/(1-a^2)
    vy = 1.0 / (1 - d * d)
    cxy = b * d * vy / (1 - a * d)
    vx = (a * b * cxy + b * b * vy + 1) / (1 - a * a)
    vxp = a * a * vx + b * b * vy + 2 * a * b * cxy + 1.0
    cx_xp = a * vx + b * cxy
    cy_xp = a * cxy + b * vy
    var_given_x = vxp - cx_xp ** 2 / vx
    var_given_xy = 1.0  # only the innovation remains once X, Y are both known
    cmi = 0.5 * np.log(var_given_x / var_given_xy)
    # cross-check CMI from the 3-variate covariance (Schur complement)
    S = np.array([[vx, cxy, cx_xp], [cxy, vy, cy_xp], [cx_xp, cy_xp, vxp]])
    v_xy = S[2, 2] - S[2, :2] @ np.linalg.solve(S[:2, :2], S[:2, 2])
    cmi_check = 0.5 * np.log(var_given_x / v_xy)
    inflation = var_given_x / var_given_xy
    print(f"{b:>6} {cmi:>12.6f} {var_given_x:>12.6f} {var_given_xy:>13.6f} "
          f"{inflation:>10.6f} {np.exp(2 * cmi):>10.6f} "
          f"{'ok' if abs(cmi - cmi_check) < 1e-10 and abs(inflation - np.exp(2 * cmi)) < 1e-9 else 'FAIL':>7}")
print("\nb=0: blanket closure — the excluded variable carries zero price (Ptolemy's")
print("MRH, in toy form). The price is exactly e^{2*CMI}: an identity, not a bound.")

print()
print("=" * 78)
print("2. THE DEBYE BOUNDARY ZONE  (two MRH equations and the region between)")
print("=" * 78)

def debye_cv(x):
    """C_V / Nk at T/theta = 1/x... computed at t = T/theta directly."""
    raise SystemExit

def cv_exact(t, n=2_000_000):
    # C_V/Nk = 9 t^3 * int_0^{1/t} x^4 e^x/(e^x-1)^2 dx  (t = T/theta_D)
    hi = 1.0 / t
    m = max(int(min(n, hi * 2000)), 2000)
    x = np.linspace(1e-12, hi, m)
    f = x ** 4 * np.exp(x) / np.expm1(x) ** 2
    return 9.0 * t ** 3 * np.trapz(f, x)

LOW_COEF = 12.0 * np.pi ** 4 / 5.0   # low-T MRH equation: C_V/Nk = LOW_COEF * t^3
print(f"\n{'T/theta':>8} {'C_V exact':>10} {'low-T T^3':>10} {'err_low%':>9} "
      f"{'3Nk (DP)':>9} {'err_high%':>9}")
rows = []
for t in (0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 5.0):
    ex = cv_exact(t)
    lo = LOW_COEF * t ** 3
    hi_ = 3.0
    rows.append((t, ex, lo, abs(ex - lo) / ex * 100, hi_, abs(ex - hi_) / ex * 100))
    print(f"{t:>8} {ex:>10.4f} {lo:>10.4f} {rows[-1][3]:>9.4f} "
          f"{hi_:>9.4f} {rows[-1][5]:>9.4f}")

def crossing(col, thresh=1.0):
    for r in rows:
        if r[col] > thresh:
            return r[0]
    return None

t_lo = crossing(3)
t_hi = None
for r in reversed(rows):
    if r[5] > 1.0:
        t_hi = r[0]
        break
print(f"\nlow-T equation's error first exceeds 1% between T/theta = "
      f"{[r[0] for r in rows if r[3] <= 1][-1]} and {t_lo}")
print(f"high-T equation's error stays above 1% down to T/theta ~ {t_hi} "
      f"(last >1% point scanned)")
print(f"C_V at T=theta_D: {cv_exact(1.0):.4f} Nk  (Dulong-Petit = 3 Nk — the high-MRH")
print("equation is already ~5% right AT the horizon edge; the T^3 law is dead there.)")
print("Loss SHAPE: truncating included modes costs ~(theta/T)^4 e^{-theta/T} — a")
print("polynomial-enhanced exponential — at low T; the failure near T~theta is the")
print("continuum approximation, not the truncation.")

print()
print("=" * 78)
print("3. N_corr AS THE CORRELATION STRUCTURE'S SHADOW")
print("=" * 78)
print("mean of N unit-variance variables, pairwise corr rho:")
print("rel. fluctuation = sqrt(rho + (1-rho)/N) ; N_corr = N/(1+rho(N-1))\n")
print(f"{'N':>6} {'rho':>6} {'rel.fluct':>10} {'1/sqrt(N)':>10} {'N_corr':>10} {'exponent*':>10}")
for N, rho in ((100, 0.0), (100, 0.01), (100, 0.1), (100, 0.5), (10000, 0.1), (1000000, 0.1)):
    rf = np.sqrt(rho + (1 - rho) / N)
    nc = N / (1 + rho * (N - 1))
    # exponent*: local slope d ln(rf) / d ln N  (analytic: -1/2 * (1-rho)/(1+rho(N-1)))
    alpha = 0.5 * (1 - rho) / (1 + rho * (N - 1))
    print(f"{N:>6} {rho:>6} {rf:>10.6f} {1 / np.sqrt(N):>10.6f} {nc:>10.2f} {alpha:>10.4f}")
print("\nrho=0: fluctuation = 1/sqrt(N), exponent -1/2 — the Gaussian fixed point,")
print("gamma = 2/sqrt(N_corr) with N_corr = N. Any rho>0: N_corr saturates at 1/rho")
print("and the exponent -> 0. The RECURRENCE of the sqrt-form across MRHs is CLT")
print("(statistics, the Gaussian fixed point of coarse-graining); the physics of each")
print("MRH lives in rho — in what the horizon's couplings to the excluded do to the")
print("effective count. Measure the exponent per domain: 1/2 = no new information;")
print("!= 1/2 = correlated DOF the MRH-only equation cannot see.")
