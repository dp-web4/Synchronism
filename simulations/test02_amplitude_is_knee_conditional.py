"""
TEST-02 (Gaia wide binaries): is the published 0.05-0.4% velocity-deviation band
consistent with the knee placements the rest of the site uses?

Coupling as Tier 1 states it:   g_eff = g_N / C(rho),  C(rho) = tanh(gamma * ln(rho/rho_crit + 1))
Local disc density (Oort-scale): rho_local = 0.09 M_sun/pc^3.

The fractional velocity excess is sqrt(1/C) - 1; the fractional ACCELERATION excess is 1/C - 1.
Tier 1's band is quoted as a velocity deviation, so both are reported.

Run: python3 test02_amplitude_is_knee_conditional.py
"""
import math

RHO_LOCAL = 0.09  # M_sun/pc^3

def C(rho, rho_crit, gamma):
    return math.tanh(gamma * math.log(rho / rho_crit + 1.0))

KNEES = [
    (0.029 * 229.0**2, "SPARC/V_flat calibration, Milky Way V_flat=229 km/s (/galaxy-plotter, /key-claims)"),
    (0.029 * 47.0**2,  "SPARC/V_flat calibration, DDO 154 V_flat=47 km/s (/galaxy-plotter)"),
    (0.161,            "measured velocity-blind knee, 2026-08-27 (/key-claims, GC fork)"),
    (0.0083,           "Refracted Gravity E0 knee (explorer 2026-09-09)"),
    (3.2e-4,           "low end of the refuted rho_c grid (PREDICTIONS.md 2026-09-09)"),
]

def solve_knee(target_vel_pct, gamma):
    """rho_crit giving a given fractional VELOCITY excess, in percent."""
    lo, hi = 1e-9, 1e6
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        v = 100.0 * (1.0 / math.sqrt(C(RHO_LOCAL, mid, gamma)) - 1.0)
        if v < target_vel_pct:
            lo = mid          # excess increases with rho_crit
        else:
            hi = mid
    return math.sqrt(lo * hi)

print("TEST-02 local boost at rho_local = %.3f M_sun/pc^3\n" % RHO_LOCAL)
print("%-12s %-7s %-12s %-14s %-14s  %s" % ("rho_crit", "gamma", "C(rho_local)", "accel excess %", "vel excess %", "knee provenance"))
for rc, tag in KNEES:
    for g in (0.489, 2.0):
        c = C(RHO_LOCAL, rc, g)
        print("%-12.4g %-7.3f %-12.5g %-14.4g %-14.4g  %s"
              % (rc, g, c, 100.0 * (1.0 / c - 1.0), 100.0 * (1.0 / math.sqrt(c) - 1.0), tag))

print("\nKnee window REQUIRED to produce Tier 1's quoted 0.05-0.4% VELOCITY band:")
for g in (0.489, 2.0):
    lo = solve_knee(0.05, g)
    hi = solve_knee(0.40, g)
    print("  gamma = %.3f :  rho_crit in [%.4g, %.4g] M_sun/pc^3" % (g, lo, hi))

print("\nSame, read as an ACCELERATION band (if the 0.05-0.4% is 1/C - 1):")
def solve_knee_acc(target_pct, gamma):
    lo, hi = 1e-9, 1e6
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        v = 100.0 * (1.0 / C(RHO_LOCAL, mid, gamma) - 1.0)
        if v < target_pct:
            lo = mid          # excess increases with rho_crit
        else:
            hi = mid
    return math.sqrt(lo * hi)
for g in (0.489, 2.0):
    print("  gamma = %.3f :  rho_crit in [%.4g, %.4g] M_sun/pc^3"
          % (g, solve_knee_acc(0.05, g), solve_knee_acc(0.40, g)))

# Presentation-rounded values, printed so a published table can be checked line-for-line
# against this artifact rather than against a mental rounding (findings_lint R3).
print("\nAS PUBLISHED (rounded for the site / proposal tables):")
for rc, tag in KNEES:
    for g in (0.489, 2.0):
        c = C(RHO_LOCAL, rc, g)
        acc = 100.0 * (1.0 / c - 1.0)
        vel = 100.0 * (1.0 / math.sqrt(c) - 1.0)
        def eng(v):
            e = math.floor(math.log10(abs(v))) if v else 0
            m = v / (10.0 ** e)
            return "%.1fe%d" % (m, e) if (e < -1 or e > 3) else "%.4g" % v
        print("  rho_crit=%-9.4g gamma=%-6.3f C=%-10s accel=%-10s %% vel=%-10s %%   %s"
              % (rc, g, eng(c), eng(acc), eng(vel), tag))
print("  (read 1.8e4 as 1.8x10^4; 2.9e-5 as 2.9x10^-5; 9.1e3 as 9.1x10^3)")

print("""
READING
-------
The quoted 0.05-0.4% band is reachable only inside a narrow knee window, and that window
does not contain ANY knee the site uses in the galaxy sector. At the published
rho_crit = 0.029*V_flat^2 calibration the Milky Way's local boost is ~3.5e6 percent in
acceleration (a factor ~3.5e4), not 0.4%. So 'practically untestable / ~80x below Gaia
reach' is a property of an otherwise-unused knee placement, not of the framework's
published calibration. At the published calibration the density branch is excluded locally
by Oort / Solar-System ephemerides by orders of magnitude -- which is the same statement
/key-claims already makes as 'rotation curves do not fail to flatten, they blow up'.
""")
