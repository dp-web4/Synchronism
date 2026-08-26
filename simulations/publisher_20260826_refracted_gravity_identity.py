"""
Is the Synchronism galaxy sector's field equation Refracted Gravity (Matsakos & Diaferio 2016)?

Clean-room verification for the Publisher lane, 2026-08-26.  Written from the PUBLISHED functional
forms only (arXiv:1603.04943 Eqs. 2.3 / 4.1, read at source); no synchronism-site script was
imported, read, or copied.  The routing note this checks is
Research/proposals/refracted_gravity_prior_art_and_density_lever_is_log_size_20260825.md.

Legs
  1. C_Omega  ==  epsilon(rho)   -- exact, symbolic, both log bases
  2. B_max    ==  1/epsilon_0    -- the bounded boost IS the vacuum floor
  3. C_rho    ~=  epsilon(eps0=0) -- numeric, regulator free, at fitted gamma AND at derived gamma

Exit 0 on all assertions passing.
"""
import numpy as np
import sympy as sp

rho, rho_c, q, e0 = sp.symbols('rho rho_c q epsilon_0', positive=True)

# ---------------------------------------------------------------- leg 1
# M&D Eq. 4.1:  eps(rho) = e0 + (1-e0) * (1/2) * ( tanh[ log((rho/rho_c)^q) ] + 1 )
# Framework  :  C_Omega  = Om + (1-Om)  * x/(1+x)
# Claim      :  identical with x = (rho/rho_c)^p,  p = 2q/ln(base),  Om = e0.
print("=== leg 1: C_Omega == epsilon(rho), exact ===")
for base, name in [(10, 'log10'), (sp.E, 'ln')]:
    p = 2 * q / sp.log(base)
    u = sp.log((rho / rho_c) ** q, base)
    eps = e0 + (1 - e0) * sp.Rational(1, 2) * (sp.tanh(u).rewrite(sp.exp) + 1)
    x = (rho / rho_c) ** p
    C_Omega = e0 + (1 - e0) * x / (1 + x)
    diff = sp.simplify(sp.powsimp(sp.expand(sp.simplify(eps - C_Omega)), force=True))
    print(f"  M&D 'log' read as {name:6s}:  p = 2q/ln(base),  eps - C_Omega = {diff}")
    assert diff == 0, f"identity failed for {name}"
print("  -> EXACT, not 2.2e-16.  The bounded ceiling function IS the gravitational permittivity.")

# ---------------------------------------------------------------- leg 2
print("\n=== leg 2: boost ceiling is the reciprocal vacuum floor ===")
x_s = sp.symbols('x', positive=True)
C_O = e0 + (1 - e0) * x_s / (1 + x_s)
floor = sp.limit(C_O, x_s, 0)
assert sp.simplify(floor - e0) == 0
print(f"  C_Omega(rho->0) = {floor};  B_max = 1/C = 1/epsilon_0")
Om = 0.315
print(f"  framework floor Omega_m = {Om}  ->  B <= {1/Om:.3f}"
      f"   |  M&D fit epsilon_0 = 0.20-0.25 -> B <= {1/0.25:.2f}-{1/0.20:.2f}")
print("  -> 'bounded boost' is a fitted parameter of the 2016 construction, not a novelty.")

# ---------------------------------------------------------------- leg 3
print("\n=== leg 3: C_rho vs epsilon at epsilon_0 = 0, regulator free ===")
r = np.logspace(-3, 4, 4000)                       # rho / rho_crit, 7 decades
pgrid = np.linspace(0.05, 3.0, 600)
sgrid = np.logspace(-2, 2, 400)                    # rho_c / rho_crit
rows = []
for gam in (0.489, 0.498, 2.0):
    C_rho = np.tanh(gam * np.log(1.0 + r))
    best = None
    for pp in pgrid:
        y = (r[None, :] / sgrid[:, None]) ** pp
        E0 = y / (1.0 + y)
        d = C_rho[None, :] - E0
        rms = np.sqrt((d * d).mean(axis=1))
        i = int(np.argmin(rms))
        if best is None or rms[i] < best[0]:
            best = (float(rms[i]), float(pp), float(sgrid[i]), float(np.abs(d[i]).max()))
    rms, pp, s, mx = best
    rows.append((gam, rms, mx, pp, s))
    print(f"  gamma={gam:<6}: rms={rms:.2e}  max={mx:.2e}  p={pp:.3f}  "
          f"rho_c/rho_crit={s:.3g}  q=p*ln10/2={pp*np.log(10)/2:.3f}  "
          f"(q / M&D disc 0.75 = {pp*np.log(10)/2/0.75:.2f}x)")

fit = [x for x in rows if x[0] == 0.489][0]
der = [x for x in rows if x[0] == 2.0][0]
assert fit[1] < 1e-3, "fitted-gamma match should be ~1e-3 rms"
assert abs(der[1] - 0.014) < 0.002 and abs(der[2] - 0.034) < 0.004, \
    "gamma=2 residual should reproduce the routing note's 0.014 / 0.032"
print("  -> the note's quoted 'rms 0.014 / max 0.032' is the gamma=2 figure;")
print(f"     at the FITTED gamma the match is {0.0142/fit[1]:.0f}x tighter (rms {fit[1]:.1e}).")
print("  -> the '+1' regulator is absorbed by a ~2x shift of rho_c: RG's published form has no regulator.")

print("\nAll assertions passed.")
