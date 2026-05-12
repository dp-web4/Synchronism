"""
Session 653: Compute C_galactic / C_cosmic using framework equations.

Visitor proposal asks: does the framework's own C(ρ) formula support the
suppressor-class mechanism Session 107 assumed?

Session 107's assumption: C_cosmic/C_galactic < 1 (cosmic coherence is
lower than galactic coherence) → effective G_local < G_global → growth
suppression at low z.

DR1 measured the opposite (enhancement). Two diagnoses:
  Branch 1: sign-flip recoverable (re-interpret coupling direction)
  Branch 2: suppressor class dead

This script computes C_galactic and C_cosmic in the framework's natural
units and reports the ratio, so the framework's own arithmetic can answer
which branch is dictated by its equations vs which requires interpretation.
"""
import numpy as np

# Physical scales (g/cm^3)
RHO_GALACTIC_OUTER = 4e-25   # V_flat regime, galaxy outskirts
RHO_COSMIC_MEAN    = 3e-30   # Ω_m × ρ_crit_cosmo

# Take ρ_crit at the galactic reference (S637 normalization)
RHO_CRIT = RHO_GALACTIC_OUTER

GAMMA = 2.0

def C(rho, rho_crit, gamma=GAMMA):
    return np.tanh(gamma * np.log(rho / rho_crit + 1.0))

C_gal = C(RHO_GALACTIC_OUTER, RHO_CRIT)
C_cos = C(RHO_COSMIC_MEAN, RHO_CRIT)

print("=" * 70)
print("SESSION 653: C_galactic vs C_cosmic in framework's natural units")
print("=" * 70)
print()
print(f"ρ_crit                = {RHO_CRIT:.2e} g/cm^3  (= ρ_galactic_outer)")
print(f"ρ_galactic            = {RHO_GALACTIC_OUTER:.2e} g/cm^3")
print(f"ρ_cosmic              = {RHO_COSMIC_MEAN:.2e} g/cm^3")
print()
print(f"ρ_cosmic / ρ_crit     = {RHO_COSMIC_MEAN / RHO_CRIT:.2e}")
print(f"ρ_galactic / ρ_crit   = {RHO_GALACTIC_OUTER / RHO_CRIT:.2e}")
print()
print(f"γ                     = {GAMMA}")
print()
print(f"C_galactic = C(ρ_galactic) = tanh(2·ln(2))     = {C_gal:.4e}")
print(f"C_cosmic   = C(ρ_cosmic)   = tanh(2·ln(1+δ))  ≈ {C_cos:.4e}")
print()
print(f"C_cosmic / C_galactic = {C_cos / C_gal:.2e}")
print(f"C_galactic / C_cosmic = {C_gal / C_cos:.2e}")
print()
print("=" * 70)
print("INTERPRETATION")
print("=" * 70)
print()
print("Session 107's mechanism: G_local/G_global = C_cosmic/C_galactic < 1")
print("                         → suppressed growth at low z (predicted)")
print()
print(f"Computed C_cosmic/C_galactic = {C_cos/C_gal:.2e}  << 1   ✓")
print()
print("So Session 107's identification of the ratio direction is")
print("CORRECT given the natural choice ρ_crit = ρ_galactic_outer.")
print("The suppression mechanism is the framework's own prediction.")
print()
print("But DR1 observed enhancement at low z. Therefore:")
print(" - Branch 2 (suppressor class dead) is the framework's own verdict")
print(" - Branch 1 (sign-flip recoverable) requires INVERTING the")
print("   coupling-to-coherence map, which is a framework REINTERPRETATION,")
print("   not a recomputation of ratios.")
print()
print("The numerical answer doesn't escape the failure; the framework's")
print("equations dictate suppression, and DR1 falsifies suppression.")
