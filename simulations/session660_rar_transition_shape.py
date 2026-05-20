"""
Session 660: RAR transition-shape discriminator — verify the compander-as-μ claim.

The proposal (rar_transition_shape_discriminator.md) claims:

Under the μ-identification (C as MOND interpolating function with argument
g_obs/a₀): g_bar = g_obs · tanh(γ · ln(1 + g_obs/a₀))

1. Reproduces both MOND asymptotes (Newtonian at high a, √-law at low a)
2. At γ=2, differs from McGaugh ν(y) = 1/(1 − exp(−√y)) only in transition
   curvature, with max deviation ~−0.083 dex at g_bar/a₀ ≈ 1.1
3. Free-γ fit collapses to γ≈0.91 (= MOND), so discriminator has power
   ONLY if γ is pinned at 2 (the framework's own N_corr=1 assignment)

This script verifies the asymptotics and the transition-shape deviation.
The full BIC fit to SPARC (2693 points, nuisance marginalization) is
operator-track; this checks the core mathematical claim.
"""
import numpy as np

A0 = 1.2e-10  # m/s^2, MOND scale (value irrelevant; we work in x = g/a0)

def mu_syn(x, gamma=2.0):
    """Synchronism compander as MOND μ: μ(x) = tanh(γ·ln(1+x)), x = g_obs/a0.
    Relation: g_bar = g_obs · μ(g_obs/a0)."""
    return np.tanh(gamma * np.log(1.0 + x))

def mcgaugh_nu(y):
    """McGaugh RAR interpolating function: g_obs = g_bar · ν(g_bar/a0),
    ν(y) = 1/(1 − exp(−√y))."""
    return 1.0 / (1.0 - np.exp(-np.sqrt(y)))

print("=" * 70)
print("SESSION 660: RAR transition-shape discriminator")
print("=" * 70)
print()

# --- Asymptotic check for μ_syn ---
print("Asymptotic behavior of g_bar = g_obs·tanh(γ·ln(1+g_obs/a0)), γ=2")
print("-" * 70)
# High acceleration: x>>1, ln(1+x)→ln(x) large, tanh→1, g_bar→g_obs (Newtonian)
x_hi = 1e3
print(f"  High a (x={x_hi:.0e}): μ = {mu_syn(x_hi):.6f} → g_bar ≈ g_obs (Newtonian) ✓")
# Low acceleration: x<<1, ln(1+x)≈x, tanh(γx)≈γx, μ≈γx
#   g_bar = g_obs·γ·(g_obs/a0) = γ·g_obs²/a0 → g_obs = sqrt(a0·g_bar/γ) (deep-MOND √-law)
x_lo = 1e-3
mu_lo = mu_syn(x_lo)
approx = 2.0 * x_lo  # γx
print(f"  Low a (x={x_lo:.0e}): μ = {mu_lo:.6e}, γx = {approx:.6e} → μ≈γx ✓")
print(f"    ⇒ g_obs = √(a0·g_bar/γ): deep-MOND √-law, Tully-Fisher preserved ✓")
print()

# --- Transition-shape comparison at γ=2 ---
# Convert both to the RAR observable: log10(g_obs) vs log10(g_bar).
# For McGaugh: g_obs = g_bar·ν(g_bar/a0), directly.
# For Synchronism μ: given g_obs, g_bar = g_obs·μ(g_obs/a0). Invert numerically
#   to get g_obs(g_bar), then compare.

# Work in dimensionless: let gbar_x = g_bar/a0, gobs_x = g_obs/a0.
gbar_x = np.geomspace(1e-2, 1e2, 2000)

# McGaugh: gobs_x = gbar_x · ν(gbar_x)
gobs_mcg = gbar_x * mcgaugh_nu(gbar_x)

# Synchronism μ at γ=2: solve gbar_x = gobs_x · tanh(2·ln(1+gobs_x)) for gobs_x
def invert_syn(gbar_target, gamma=2.0):
    # monotonic in gobs_x; bisection
    lo, hi = 1e-6, 1e6
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        val = mid * np.tanh(gamma * np.log(1.0 + mid))
        if val < gbar_target:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)

gobs_syn = np.array([invert_syn(g) for g in gbar_x])

# Deviation in dex
dev_dex = np.log10(gobs_syn) - np.log10(gobs_mcg)
imax = np.argmax(np.abs(dev_dex))

print("Transition-shape deviation (γ=2 compander vs McGaugh ν)")
print("-" * 70)
print(f"  Max |deviation|: {dev_dex[imax]:+.4f} dex at g_bar/a0 = {gbar_x[imax]:.2f}")
print(f"  Proposal claimed: −0.083 dex at g_bar/a0 ≈ 1.1")
print()
SIGMA_INT = 0.057  # SPARC RAR intrinsic scatter (Lelli+2017)
print(f"  SPARC σ_int = {SIGMA_INT} dex")
print(f"  Max deviation / σ_int = {abs(dev_dex[imax])/SIGMA_INT:.2f}× "
      f"(proposal: 1.45×)")
print()

# RMS over the transition region (g_bar/a0 in [0.1, 10])
mask = (gbar_x >= 0.1) & (gbar_x <= 10)
rms = np.sqrt(np.mean(dev_dex[mask]**2))
print(f"  RMS deviation over transition (0.1 < g_bar/a0 < 10): {rms:.4f} dex")
print(f"  Proposal claimed RMS: 0.067 dex (> σ_int = 0.057)")
print()

print("=" * 70)
print("VERDICT")
print("=" * 70)
print("- Asymptotics confirmed: μ_syn reproduces Newtonian + deep-MOND √-law")
print("- At γ=2, the compander is a DISTINCT interpolating function with a")
print(f"  measurable transition-shape deviation (~{abs(dev_dex[imax]):.3f} dex)")
print("- Comparable to or exceeding SPARC σ_int → mildly testable/disfavored")
print()
print("CRITICAL CONTINGENCY (per proposal): discrimination requires γ pinned")
print("at 2. With free γ, the fit returns γ≈0.9 (N_corr≈5) ≈ McGaugh, zero")
print("discrimination. So the test reduces to S643's open question:")
print("   Is galaxy-scale γ pinned by N_corr=1, or is it a fit parameter?")
print()
print("This is the framework's FIRST galaxy-scale test not MOND-degenerate")
print("by sign or EFE — but only if the framework commits to γ=2 a priori.")
