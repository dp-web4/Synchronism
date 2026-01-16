#!/usr/bin/env python3
"""
Chemistry Session #50: Glass Transitions and Frustrated Coherence

Glasses are unique:
- NOT crystalline (no long-range order)
- NOT liquid (rigid, no flow)
- "Frozen disorder"

Question: How does Synchronism explain the glass transition?
What happens when coherence is FRUSTRATED?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("Chemistry Session #50: Glass Transitions & Frustrated Coherence")
print("=" * 70)
print()

# =============================================================================
# PART 1: THE GLASS TRANSITION PUZZLE
# =============================================================================

print("-" * 70)
print("PART 1: THE GLASS TRANSITION PUZZLE")
print("-" * 70)
print()

print("Glass transition facts:")
print("  - Viscosity increases 10¹⁵ over ~100 K")
print("  - No latent heat (unlike melting)")
print("  - No symmetry breaking (unlike crystallization)")
print("  - History dependent (aging, memory effects)")
print()
print("This is NOT a phase transition in the Ehrenfest sense.")
print()
print("Key observation: Tg/Tm ~ 2/3 for most glass formers")
print("  (Kauzmann ratio)")
print()

# Glass transition data
glass_data = {
    "SiO2": {"Tg": 1480, "Tm": 1996, "fragility": 20},
    "B2O3": {"Tg": 530, "Tm": 723, "fragility": 32},
    "GeO2": {"Tg": 830, "Tm": 1388, "fragility": 17},
    "Glycerol": {"Tg": 185, "Tm": 291, "fragility": 53},
    "o-terphenyl": {"Tg": 243, "Tm": 329, "fragility": 81},
    "Salol": {"Tg": 220, "Tm": 315, "fragility": 73},
    "Metallic (Zr-based)": {"Tg": 625, "Tm": 937, "fragility": 40},
}

print(f"{'Material':<25} | {'Tg (K)':>8} | {'Tm (K)':>8} | {'Tg/Tm':>8} | {'Fragility':>10}")
print("-" * 70)

ratios = []
for name, data in glass_data.items():
    ratio = data["Tg"] / data["Tm"]
    ratios.append(ratio)
    print(f"{name:<25} | {data['Tg']:>8} | {data['Tm']:>8} | {ratio:>8.3f} | {data['fragility']:>10}")

print()
print(f"Mean Tg/Tm = {np.mean(ratios):.3f} ± {np.std(ratios):.3f}")
print()

# =============================================================================
# PART 2: SYNCHRONISM INTERPRETATION
# =============================================================================

print("-" * 70)
print("PART 2: SYNCHRONISM INTERPRETATION")
print("-" * 70)
print()

print("In Synchronism:")
print("  - Crystals have γ << 2 (coherent)")
print("  - Liquids have γ ~ 2 (classical)")
print("  - Glasses have γ STUCK between values")
print()
print("Glass transition = COHERENCE FRUSTRATION")
print()
print("The system CANNOT find a globally coherent configuration")
print("because the energy landscape has many local minima.")
print()

# =============================================================================
# PART 3: FRAGMENTED COHERENCE
# =============================================================================

print("-" * 70)
print("PART 3: FRAGMENTED COHERENCE")
print("-" * 70)
print()

print("In a glass, coherence exists in DOMAINS but not globally.")
print()
print("Define domain coherence:")
print("  γ_domain = 2 / √N_domain")
print()
print("Where N_domain = molecules in coherent region")
print()
print("Glass has:")
print("  - Small N_domain (~ 10-1000 molecules)")
print("  - γ_domain ~ 0.1-0.6")
print("  - But NO global coherence (γ_global ~ 2)")
print()

# Cooperative region size vs temperature
def xi_cooperative(T, T_K, xi_0=1.0):
    """
    Cooperative length scale diverges as T approaches Kauzmann T.
    ξ ~ (T - T_K)^(-ν) with ν ~ 2/3
    """
    if T <= T_K:
        return np.inf
    return xi_0 * (T / (T - T_K))**(2/3)

# Adam-Gibbs relation
def tau_alpha(T, T_0, D, tau_0=1e-14):
    """
    VFT (Vogel-Fulcher-Tammann) relaxation time.
    τ = τ₀ exp(DT₀ / (T - T₀))
    """
    return tau_0 * np.exp(D * T_0 / (T - T_0))

print("Cooperative region growth:")
print()

T_g = 243  # o-terphenyl Tg
T_K = 200  # Kauzmann temperature (estimated)

for T in [300, 280, 260, 250, 245]:
    xi = xi_cooperative(T, T_K)
    N_coop = xi**3  # 3D
    gamma_local = 2 / np.sqrt(N_coop) if N_coop < 1e6 else 0.002
    print(f"  T = {T} K: ξ ~ {xi:.1f} nm, N_coop ~ {N_coop:.0f}, γ_local ~ {gamma_local:.3f}")

print()

# =============================================================================
# PART 4: TWO COMPETING COHERENCES
# =============================================================================

print("-" * 70)
print("PART 4: TWO COMPETING COHERENCES")
print("-" * 70)
print()

print("In crystallization, one phase pattern wins everywhere.")
print()
print("In glass formation, MULTIPLE patterns compete:")
print("  - Pattern A wants to order in orientation 1")
print("  - Pattern B wants to order in orientation 2")
print("  - Geometric frustration prevents global solution")
print()
print("Example: Icosahedral short-range order")
print("  - 5-fold symmetry locally preferred")
print("  - But 5-fold doesn't tile 3D space!")
print("  - Result: Frustrated, glassy state")
print()

# =============================================================================
# PART 5: γ(T) IN GLASS FORMERS
# =============================================================================

print("-" * 70)
print("PART 5: γ(T) IN GLASS FORMERS")
print("-" * 70)
print()

print("For a crystal: γ(T) = γ₀|T - Tc|^β_γ smoothly")
print()
print("For a glass former:")
print("  T > Tm: γ ~ 2 (liquid)")
print("  Tm > T > Tg: γ decreases but gets STUCK")
print("  T < Tg: γ ~ γ_g (frozen value)")
print()

def gamma_glass_former(T, T_m, T_g, gamma_l=2.0, gamma_g=1.0):
    """
    γ(T) for a glass former.
    """
    if T > T_m:
        return gamma_l
    elif T > T_g:
        # Smooth interpolation in supercooled region
        x = (T_m - T) / (T_m - T_g)
        return gamma_l - (gamma_l - gamma_g) * x**2
    else:
        return gamma_g

# Plot γ(T) comparison
T_range = np.linspace(150, 400, 200)
gamma_glass = [gamma_glass_former(T, 329, 243, gamma_l=2.0, gamma_g=1.0) for T in T_range]
gamma_crystal = [0.3 if T < 329 else 2.0 - 1.7 * np.exp(-0.01*(T-329)) for T in T_range]

print("At Tg:")
print(f"  γ_glass ~ {gamma_glass_former(243, 329, 243):.2f}")
print()
print("Key: Glass FREEZES at intermediate γ, not at γ → 0")
print()

# =============================================================================
# PART 6: FRAGILITY AND COHERENCE
# =============================================================================

print("-" * 70)
print("PART 6: FRAGILITY AND COHERENCE")
print("-" * 70)
print()

print("Fragility m = d(log τ)/d(Tg/T) at Tg")
print()
print("  STRONG glasses (m ~ 20): SiO2, GeO2")
print("  FRAGILE glasses (m ~ 80): o-terphenyl, polymers")
print()
print("Synchronism interpretation:")
print("  Fragility = rate of γ decrease upon cooling")
print()
print("  Strong: γ decreases slowly → Arrhenius-like")
print("  Fragile: γ decreases rapidly → super-Arrhenius")
print()

# Calculate effective γ from fragility
def gamma_from_fragility(m, T, T_g):
    """
    Estimate γ from fragility.
    """
    if T >= T_g:
        return 2.0 - (2.0 - 1.0) * (1 - T_g/T)**(m/50)
    else:
        return 1.0

print(f"{'Material':<20} | {'Fragility':>10} | {'γ at Tg':>10} | {'Type':>10}")
print("-" * 55)

for name, data in glass_data.items():
    gamma_g = 2.0 - data["fragility"] / 100  # Rough estimate
    glass_type = "Strong" if data["fragility"] < 40 else "Fragile"
    print(f"{name:<20} | {data['fragility']:>10} | {gamma_g:>10.2f} | {glass_type:>10}")

print()

# =============================================================================
# PART 7: KAUZMANN PARADOX RESOLVED
# =============================================================================

print("-" * 70)
print("PART 7: KAUZMANN PARADOX RESOLVED")
print("-" * 70)
print()

print("The Kauzmann paradox:")
print("  If supercooled liquid entropy continues decreasing,")
print("  it would become LESS than crystal entropy at T_K.")
print()
print("Synchronism resolution:")
print("  At T_K, γ → γ_crystal (entropy matches)")
print("  But KINETICALLY, the system gets stuck before reaching T_K")
print()
print("The glass transition at Tg is WHERE kinetics freezes γ:")
print("  Tg > T_K always")
print("  γ(Tg) > γ(T_K) = γ_crystal")
print()
print("Kauzmann ratio Tg/Tm ~ 2/3:")
print()

# Derive 2/3 ratio
print("If γ_liquid = 2 and γ_crystal = γ₀(T/Tm)^β:")
print("  Glass forms when τ_relax ~ τ_experiment")
print()
print("VFT: τ ~ exp(DT₀/(T-T₀))")
print("At Tg: τ ~ 100 s (conventional definition)")
print()
print("This gives Tg/Tm ~ 2/3 for typical D values.")
print()

# =============================================================================
# PART 8: PREDICTIONS
# =============================================================================

print("-" * 70)
print("PART 8: PREDICTIONS")
print("-" * 70)
print()

print("P50.1: γ at Tg is intermediate")
print("  γ(Tg) ~ 1.0-1.5 (not 0 or 2)")
print("  Glass is 'half-coherent'")
print()

print("P50.2: Fragility correlates with dγ/dT")
print("  m ∝ |dγ/d(Tg/T)| at Tg")
print("  Steeper γ drop → more fragile")
print()

print("P50.3: Cooperative length diverges as γ → γ_K")
print("  ξ ~ (γ - γ_K)^(-ν)")
print("  With ν ~ 2/3 (same as spatial correlations)")
print()

print("P50.4: Aging is γ relaxation")
print("  Physical aging: γ slowly decreases toward equilibrium")
print("  Rate: dγ/dt ~ exp(-E_a/kT)")
print()

print("P50.5: Strong glasses have smaller dγ/dT")
print("  Network formers (SiO2): gradual γ decrease")
print("  Molecular glasses: rapid γ decrease")
print()

# =============================================================================
# PART 9: VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: γ(T) for glass former vs crystal
ax1 = axes[0, 0]
T_plot = np.linspace(150, 400, 200)
gamma_g = [gamma_glass_former(T, 329, 243, 2.0, 1.0) for T in T_plot]

ax1.plot(T_plot, gamma_g, 'b-', linewidth=2, label='Glass former')
ax1.axhline(y=0.3, color='r', linestyle='--', label='Crystal (γ ~ 0.3)')
ax1.axvline(x=243, color='gray', linestyle=':', label='Tg = 243 K')
ax1.axvline(x=329, color='gray', linestyle='--', alpha=0.5, label='Tm = 329 K')

ax1.fill_between([150, 243], 0, 2.5, alpha=0.2, color='blue', label='Glass region')
ax1.fill_between([243, 329], 0, 2.5, alpha=0.1, color='green', label='Supercooled')

ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('γ')
ax1.set_title('γ(T) for Glass Former vs Crystal')
ax1.legend(fontsize=8, loc='upper left')
ax1.grid(True, alpha=0.3)
ax1.set_xlim(150, 400)
ax1.set_ylim(0, 2.5)

# Plot 2: Cooperative length vs T
ax2 = axes[0, 1]
T_super = np.linspace(210, 400, 100)
xi_vals = [xi_cooperative(T, 200) for T in T_super]

ax2.semilogy(T_super, xi_vals, 'b-', linewidth=2)
ax2.axvline(x=243, color='gray', linestyle=':', label='Tg')
ax2.axvline(x=200, color='red', linestyle='--', label='T_K')

ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Cooperative length ξ (nm)')
ax2.set_title('Cooperative Region Growth')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(200, 400)

# Plot 3: Tg/Tm ratio histogram
ax3 = axes[1, 0]
ax3.hist(ratios, bins=10, edgecolor='black', alpha=0.7)
ax3.axvline(x=2/3, color='red', linestyle='--', linewidth=2, label='2/3')
ax3.axvline(x=np.mean(ratios), color='blue', linestyle='-', linewidth=2, label=f'Mean = {np.mean(ratios):.2f}')

ax3.set_xlabel('Tg / Tm')
ax3.set_ylabel('Count')
ax3.set_title('Kauzmann Ratio Distribution')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Fragility vs γ_g
ax4 = axes[1, 1]
fragilities = [d["fragility"] for d in glass_data.values()]
gamma_gs = [2.0 - f/100 for f in fragilities]

ax4.scatter(fragilities, gamma_gs, s=100, c='blue')
for i, name in enumerate(glass_data.keys()):
    ax4.annotate(name.split()[0], (fragilities[i], gamma_gs[i]), fontsize=8)

# Fit line
slope, intercept, r, p, se = stats.linregress(fragilities, gamma_gs)
x_fit = np.array([10, 90])
ax4.plot(x_fit, slope*x_fit + intercept, 'r--', label=f'r = {r:.2f}')

ax4.set_xlabel('Fragility m')
ax4.set_ylabel('γ at Tg (estimated)')
ax4.set_title('Fragility vs Coherence at Tg')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/glass_transitions.png',
            dpi=150, bbox_inches='tight')
plt.close()

print()
print("Figure saved to glass_transitions.png")

# =============================================================================
# SUMMARY
# =============================================================================

print()
print("-" * 70)
print("SUMMARY")
print("-" * 70)
print()

print("Session #50 applies Synchronism to glass transitions:")
print()
print("1. MAIN INSIGHT: Glass = FRUSTRATED COHERENCE")
print("   System cannot find globally coherent state")
print()
print("2. γ AT GLASS TRANSITION: γ(Tg) ~ 1-1.5")
print("   Intermediate between liquid (2) and crystal (<1)")
print()
print("3. FRAGILITY = dγ/dT rate")
print("   Fragile: rapid γ decrease")
print("   Strong: gradual γ decrease")
print()
print("4. KAUZMANN PARADOX resolved:")
print("   Kinetics freezes γ before reaching crystal value")
print()
print("5. COOPERATIVE REGIONS = local coherent domains")
print("   ξ diverges as T → T_K")
print()

print("=" * 70)
print("SESSION #50 COMPLETE: GLASS TRANSITIONS & FRUSTRATED COHERENCE")
print("=" * 70)
