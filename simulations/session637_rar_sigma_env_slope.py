"""
Session 637: RAR σ_int(ρ_env) slope — derivation attempt

Visitor proposal asks: derive a numerical slope for the framework's claim
that σ_int(RAR) depends on local environmental density ρ_env.

Framework's tools:
  - C(ρ) = tanh(γ · log(ρ/ρ_crit + 1)) with γ = 2 (from 6D phase space)
  - V²_total = V²_baryon / C(ρ)  (Session #65, local form)
  - ρ_crit = A · V_flat²,   A = 0.029 (km/s)⁻² (Session #66, α=1.0 fiducial)

Question: how does ρ_env enter C(ρ)? The framework's derivation defines
ρ in C(ρ) as "local density at radius r" (Session #64). In the V_flat
regime (galaxy outskirts), the local density is dominated by the galaxy
itself; environmental density is a small additive perturbation:

  ρ_total(r) = ρ_galactic(r) + ρ_env

Propagating variance: σ_int(RAR bin) ≈ |d(log C)/d(log ρ)|_<ρ> × σ_(log ρ|bin).

For the "environment-dependent scatter" claim to be a real prediction,
the difference Δρ_env between bins must produce a detectable Δσ_int.
"""
import numpy as np

# Physical constants and typical scales
RHO_COSMIC_MEAN = 3e-30      # g/cm^3 (Ω_m * ρ_crit_cosmo)
RHO_ENV_CLUSTER = 200 * RHO_COSMIC_MEAN   # ~6e-28 g/cm^3 (within R_200)
RHO_ENV_GROUP   = 50  * RHO_COSMIC_MEAN
RHO_ENV_VOID    = 0.1 * RHO_COSMIC_MEAN

# Galactic outer-disk mean density (V_flat regime)
# M ~ V^2 R / G; V=200 km/s, R=20 kpc gives M~2e11 M_sun in sphere, ρ_avg~4e-25 g/cm^3
RHO_GALACTIC_OUTER = 4e-25   # g/cm^3, at galaxy outskirts

# Coherence parameters
GAMMA = 2.0  # derived from 3+3-4 = 2 (Session #64)


def C(rho, rho_crit, gamma=GAMMA):
    """Coherence function."""
    return np.tanh(gamma * np.log(rho / rho_crit + 1.0))


def dlogC_dlogrho(rho, rho_crit, gamma=GAMMA):
    """d(log C)/d(log ρ) — the slope that propagates variance."""
    x = rho / rho_crit
    y = np.log(x + 1.0)
    Cv = np.tanh(gamma * y)
    sech2 = 1.0 - Cv**2
    # dC/dy = gamma * sech^2(gamma y); dy/dlogx = x/(x+1); dlogC/dlogx = (1/C)*dC/dy*dy/dlogx
    return (gamma / Cv) * sech2 * (x / (x + 1.0))


def pretty(label, val, fmt="{:.4g}"):
    print(f"  {label:<48s} {fmt.format(val)}")


print("=" * 72)
print("SESSION 637: RAR σ_int(ρ_env) slope derivation")
print("=" * 72)
print()
print("Setup")
print("-" * 72)
pretty("γ (coherence sharpness)", GAMMA)
pretty("ρ_galactic_outer (V_flat regime, g/cm^3)", RHO_GALACTIC_OUTER)
pretty("ρ_env (cluster outskirts, g/cm^3)", RHO_ENV_CLUSTER)
pretty("ρ_env (group, g/cm^3)", RHO_ENV_GROUP)
pretty("ρ_env (void, g/cm^3)", RHO_ENV_VOID)
print()

# We need ρ_crit in g/cm^3 to compare. The framework gives ρ_crit = A*V^2 with
# A = 4π/(α²GR₀²) ≈ 0.029 (km/s)⁻². For V_flat = 200 km/s:
#   ρ_crit = 0.029 * 200^2 = 1160 (a number; needs a base scale to be dimensional)
# Session #66 treats ρ as dimensionless (ρ/I_max), so ρ_crit and ρ are in
# whatever units the local density is measured. Take the framework's
# self-consistent reading: ρ at galaxy outskirts is the unit-1 reference, and
# ρ_crit = 1 there. Then ρ_env / ρ_crit = ρ_env / ρ_galactic_outer.

print("Local density (galactic) >> environmental density at galaxy outskirts")
print("-" * 72)
for label, rho_env in [
    ("isolated/void", RHO_ENV_VOID),
    ("group",        RHO_ENV_GROUP),
    ("cluster",      RHO_ENV_CLUSTER),
]:
    ratio = rho_env / RHO_GALACTIC_OUTER
    pretty(f"ρ_env({label}) / ρ_galactic_outer", ratio, "{:.2e}")
print()

# Total local density at galaxy outskirts in different environments
# Take ρ_crit = ρ_galactic_outer (framework normalization).
RHO_CRIT = RHO_GALACTIC_OUTER

results = []
for label, rho_env in [
    ("isolated/void", RHO_ENV_VOID),
    ("group",        RHO_ENV_GROUP),
    ("cluster",      RHO_ENV_CLUSTER),
]:
    rho_total = RHO_GALACTIC_OUTER + rho_env
    Cval = C(rho_total, RHO_CRIT)
    slope = dlogC_dlogrho(rho_total, RHO_CRIT)
    log_ratio = np.log10(rho_total / RHO_GALACTIC_OUTER)
    results.append((label, rho_env, rho_total, Cval, slope, log_ratio))

print("Coherence and slope at galaxy outskirts in each environment")
print("-" * 72)
print(f"  {'environment':<16s} {'C(ρ)':>10s} {'dlogC/dlogρ':>14s} {'Δlog ρ':>12s}")
for label, rho_env, rho_total, Cval, slope, log_ratio in results:
    print(f"  {label:<16s} {Cval:>10.4f} {slope:>14.4f} {log_ratio:>12.2e}")
print()

# Effect on σ_int: σ_int(bin) ≈ |slope| · σ_(log ρ|bin) where σ_(log ρ|bin)
# is the standard deviation of log ρ across galaxies within that environment bin.
# Most of σ comes from galaxy-to-galaxy variation in ρ_galactic, which dominates
# ρ_env. A typical scatter in ρ_galactic_outer across SPARC is ~0.5 dex.

SIGMA_LOG_RHO_GAL = 0.5  # rough scatter in galactic outer density across SPARC

print("Predicted σ_int from environmental ρ variation alone")
print("-" * 72)
print("Per-bin scatter from environment is set by Δρ_env between bins,")
print("not σ within bins (which is dominated by ρ_galactic).")
print()
slopes = [s for _, _, _, _, s, _ in results]
mean_slope = np.mean(slopes)
log_ratio_cluster_void = np.log10(
    (RHO_GALACTIC_OUTER + RHO_ENV_CLUSTER) / (RHO_GALACTIC_OUTER + RHO_ENV_VOID)
)
delta_logC_cluster_void = mean_slope * log_ratio_cluster_void
pretty("mean d(log C)/d(log ρ)", mean_slope)
pretty("Δlog ρ from void to cluster (dex)", log_ratio_cluster_void)
pretty("Predicted Δlog C (cluster − void), dex", delta_logC_cluster_void)
pretty("|Δσ_int| (cluster − void), dex", abs(delta_logC_cluster_void))
print()

# Compare to baseline σ_int and SPARC measurement floor
SIGMA_INT_BASELINE = 0.13   # Lelli+2017
SPARC_FLOOR        = 0.02   # per-bin precision

print("Comparison")
print("-" * 72)
pretty("σ_int RAR baseline (Lelli+2017), dex", SIGMA_INT_BASELINE)
pretty("SPARC σ_int per-bin precision, dex",  SPARC_FLOOR)
pretty("Predicted Δσ_int as fraction of baseline",
       abs(delta_logC_cluster_void) / SIGMA_INT_BASELINE,
       "{:.1%}")
pretty("Predicted Δσ_int as fraction of SPARC floor",
       abs(delta_logC_cluster_void) / SPARC_FLOOR,
       "{:.1%}")
print()
print("=" * 72)
if abs(delta_logC_cluster_void) < SPARC_FLOOR:
    print("VERDICT: predicted environmental Δσ_int is BELOW SPARC measurement")
    print("         floor. The 'environment-dependent scatter' claim, when")
    print("         derived from C(ρ) with γ=2, is ~%.0f%% of the per-bin floor."
          % (100 * abs(delta_logC_cluster_void) / SPARC_FLOOR))
    print("         Cannot be discriminated from σ_int = const (MOND).")
else:
    print("VERDICT: predicted Δσ_int is potentially detectable.")
print("=" * 72)
