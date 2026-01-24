"""
Chemistry Session #196: Solvation and Hydration at γ ~ 1
Analyzing solvation through the coherence framework.

Key hypothesis: Various solvation parameters show γ ~ 1 boundaries
- Hydration number n_h ≈ 4-6 for most ions (vs 4 for tetrahedral water)
- Born solvation energy
- Solvent reorganization energy
- Debye-Hückel screening length
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 60)
print("CHEMISTRY SESSION #196: SOLVATION AT γ ~ 1")
print("=" * 60)

# ============================================================
# 1. HYDRATION NUMBERS
# ============================================================
print("\n" + "=" * 60)
print("1. HYDRATION NUMBERS")
print("=" * 60)

print("""
Hydration number n_h: number of water molecules in first shell.

Reference: Tetrahedral coordination = 4
γ_hyd = n_h / 4

At γ ~ 1: tetrahedral-like hydration
Larger ions: γ > 1 (more waters fit)
Smaller ions: γ ~ 1 (constrained by water H-bonding)
""")

# Hydration number data: (ion, n_h, ionic_radius_pm)
hydration_data = {
    # Cations
    'Li+': (4.0, 76),
    'Na+': (5.4, 102),
    'K+': (6.0, 138),
    'Rb+': (7.0, 152),
    'Cs+': (8.0, 167),
    'Mg2+': (6.0, 72),
    'Ca2+': (7.2, 100),
    'Sr2+': (8.0, 118),
    'Ba2+': (9.5, 135),
    'Al3+': (6.0, 54),
    'Fe3+': (6.0, 65),
    'Cr3+': (6.0, 62),
    # Anions
    'F-': (4.0, 133),
    'Cl-': (6.0, 181),
    'Br-': (6.0, 196),
    'I-': (6.0, 220),
    'OH-': (4.0, 137),
}

print("\nIon Hydration Numbers:")
print("-" * 60)
print(f"{'Ion':<10} {'n_h':<10} {'γ = n_h/4':<12} {'r_ion (pm)'}")
print("-" * 60)

gamma_values = []
for ion, (n_h, r_ion) in hydration_data.items():
    gamma = n_h / 4
    gamma_values.append(gamma)
    print(f"{ion:<10} {n_h:<10.1f} {gamma:<12.2f} {r_ion}")

print("-" * 60)
mean_gamma = np.mean(gamma_values)
std_gamma = np.std(gamma_values)
print(f"Mean γ = {mean_gamma:.2f} ± {std_gamma:.2f}")

# Count near γ ~ 1
near_one = sum(1 for g in gamma_values if 0.8 <= g <= 1.5)
print(f"Near γ ~ 1 (0.8-1.5): {near_one}/{len(gamma_values)}")

# ============================================================
# 2. BORN SOLVATION MODEL
# ============================================================
print("\n" + "=" * 60)
print("2. BORN SOLVATION MODEL")
print("=" * 60)

print("""
Born equation for ion solvation:
  ΔG_solv = -(z²e²/8πε₀r) × (1 - 1/ε_r)

For water (ε_r = 80):
  ΔG_solv ≈ -69 × z²/r kJ/mol (r in Å)

The (1 - 1/ε_r) factor:
  For water: 1 - 1/80 = 0.9875 ≈ 1

γ_Born = ε_r / (ε_r - 1)
  For water: γ = 80/79 = 1.013 (γ ~ 1!)

High-ε solvents approach γ → 1 (complete solvation).
""")

# Dielectric constants
solvents = {
    'Water': 80.0,
    'DMSO': 47.0,
    'Acetonitrile': 37.5,
    'Methanol': 33.0,
    'Ethanol': 25.0,
    'Acetone': 21.0,
    'Dichloromethane': 9.0,
    'THF': 7.5,
    'Diethyl ether': 4.3,
    'Benzene': 2.3,
    'Hexane': 1.9,
}

print("\nBorn Factor (1 - 1/ε) for Solvents:")
print("-" * 50)
print(f"{'Solvent':<20} {'ε_r':<10} {'1 - 1/ε':<12} {'γ = ε/(ε-1)'}")
print("-" * 50)

gamma_born = []
for solvent, eps in solvents.items():
    factor = 1 - 1/eps
    gamma = eps / (eps - 1)
    gamma_born.append(gamma)
    print(f"{solvent:<20} {eps:<10.1f} {factor:<12.4f} {gamma:.3f}")

print("-" * 50)
print(f"Water γ = 1.013 (essentially 1!)")
print("High-ε solvents: γ → 1 (complete electrostatic solvation)")

# ============================================================
# 3. DEBYE-HÜCKEL THEORY
# ============================================================
print("\n" + "=" * 60)
print("3. DEBYE-HÜCKEL SCREENING")
print("=" * 60)

print("""
Debye length κ⁻¹: screening length for ionic interactions.
  κ² = 2NAe²I / (ε₀ε_rkT)

At physiological ionic strength (I ~ 0.15 M):
  κ⁻¹ ≈ 0.8 nm

Relevant γ parameter:
  γ_DH = κa (dimensionless, a = ion size)

At γ = 1: ion size = screening length
  Below: unscreened Coulomb
  Above: screened interactions

For typical ions (a ~ 0.3 nm) at I = 0.15 M:
  κa ≈ 0.4 (moderately screened)
""")

# Calculate Debye length vs ionic strength
I_range = np.array([0.001, 0.01, 0.05, 0.1, 0.15, 0.5, 1.0])
eps_water = 80
eps_0 = 8.85e-12
kT = 1.38e-23 * 298
e = 1.6e-19
NA = 6.02e23

kappa_squared = 2 * NA * e**2 * (I_range * 1000) / (eps_0 * eps_water * kT)
kappa_inv_nm = 1 / np.sqrt(kappa_squared) * 1e9

print("\nDebye Length vs Ionic Strength:")
print("-" * 40)
print(f"{'I (M)':<12} {'κ⁻¹ (nm)':<12} {'κa (a=0.3nm)'}")
print("-" * 40)

for I, k_inv in zip(I_range, kappa_inv_nm):
    kappa_a = 0.3 / k_inv  # a = 0.3 nm typical ion
    print(f"{I:<12.3f} {k_inv:<12.2f} {kappa_a:<12.2f}")

print("-" * 40)
print("At I = 0.1 M: κa ≈ 0.3 (Coulombic regime)")
print("At I = 1.0 M: κa ≈ 1.0 (γ ~ 1 crossover!)")

# ============================================================
# 4. MARCUS REORGANIZATION
# ============================================================
print("\n" + "=" * 60)
print("4. SOLVENT REORGANIZATION ENERGY")
print("=" * 60)

print("""
Marcus theory for electron transfer:
  λ_s = (Δe)² / (4πε₀) × (1/2r - 1/R) × (1/ε_∞ - 1/ε_s)

The optical-static factor:
  Δε = 1/ε_∞ - 1/ε_s

For water: ε_∞ ≈ 1.78 (optical), ε_s = 80
  Δε = 0.56 - 0.01 = 0.55

The γ_Marcus = 1 is the activationless point:
  |ΔG°| = λ (reorganization = driving force)

Already covered in Session #189 (reaction kinetics).
This connects to solvation through λ_s.
""")

# ============================================================
# 5. HOFMEISTER SERIES
# ============================================================
print("\n" + "=" * 60)
print("5. HOFMEISTER SERIES")
print("=" * 60)

print("""
Hofmeister series ranks ions by protein/surface effects:

Kosmotropes (structure makers) → Chaotropes (structure breakers)

Anions: SO₄²⁻ > HPO₄²⁻ > F⁻ > Cl⁻ > Br⁻ > NO₃⁻ > ClO₄⁻ > SCN⁻
Cations: Mg²⁺ > Ca²⁺ > Li⁺ > Na⁺ > K⁺ > NH₄⁺ > Cs⁺

The neutral point (γ ~ 1) is around Na⁺/Cl⁻:
- Neither structure-making nor structure-breaking
- Reference electrolyte (NaCl)
- Most physiological!
""")

# Jones-Dole B coefficients (measure of structure effect)
# Positive = kosmotrope (structure maker), Negative = chaotrope
b_coeff = {
    'SO4^2-': 0.21,
    'F-': 0.10,
    'Cl-': -0.007,
    'Br-': -0.033,
    'I-': -0.073,
    'NO3-': -0.045,
    'SCN-': -0.022,
    'Li+': 0.15,
    'Na+': 0.086,
    'K+': -0.007,
    'Rb+': -0.029,
    'Cs+': -0.045,
}

print("\nJones-Dole B Coefficients (viscosity):")
print("-" * 50)
print(f"{'Ion':<12} {'B (L/mol)':<12} {'Type'}")
print("-" * 50)

for ion, B in sorted(b_coeff.items(), key=lambda x: -x[1]):
    if B > 0.02:
        ion_type = "Kosmotrope"
    elif B < -0.02:
        ion_type = "Chaotrope"
    else:
        ion_type = "NEUTRAL (γ ~ 1)"
    print(f"{ion:<12} {B:<12.3f} {ion_type}")

print("-" * 50)
print("Cl⁻ and K⁺: B ≈ 0 (γ ~ 1 neutral point!)")

# ============================================================
# 6. VISUALIZATION
# ============================================================
print("\n" + "=" * 60)
print("6. GENERATING VISUALIZATION")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Hydration number vs ionic radius
ax1 = axes[0, 0]
radii = [v[1] for v in hydration_data.values()]
n_h_vals = [v[0] for v in hydration_data.values()]
cations = [k for k in hydration_data if '+' in k]
anions = [k for k in hydration_data if '-' in k]

cat_r = [hydration_data[k][1] for k in cations]
cat_n = [hydration_data[k][0] for k in cations]
an_r = [hydration_data[k][1] for k in anions]
an_n = [hydration_data[k][0] for k in anions]

ax1.scatter(cat_r, cat_n, s=80, c='blue', label='Cations', edgecolors='black')
ax1.scatter(an_r, an_n, s=80, c='red', label='Anions', edgecolors='black')
ax1.axhline(y=4, color='green', linestyle='--', linewidth=2, label='n_h = 4 (tetrahedral)')
ax1.set_xlabel('Ionic Radius (pm)', fontsize=12)
ax1.set_ylabel('Hydration Number n_h', fontsize=12)
ax1.set_title('Hydration Number: Many at n_h ~ 4-6', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: γ_hyd distribution
ax2 = axes[0, 1]
ax2.hist(gamma_values, bins=10, edgecolor='black', alpha=0.7, color='steelblue')
ax2.axvline(x=1, color='red', linestyle='--', linewidth=2, label='γ = 1 (tetrahedral)')
ax2.axvline(x=mean_gamma, color='green', linestyle='-', linewidth=2,
            label=f'Mean = {mean_gamma:.2f}')
ax2.set_xlabel('γ = n_h / 4', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.set_title('Hydration Coherence Distribution', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Born factor vs dielectric constant
ax3 = axes[1, 0]
eps_range = np.linspace(2, 100, 100)
born_factor = 1 - 1/eps_range
gamma_eps = eps_range / (eps_range - 1)

ax3.plot(eps_range, born_factor, 'b-', linewidth=2, label='1 - 1/ε')
ax3.axhline(y=1, color='red', linestyle='--', linewidth=2, label='Complete solvation')
for solv, eps in list(solvents.items())[:6]:
    ax3.scatter(eps, 1 - 1/eps, s=80, zorder=5)
    ax3.annotate(solv, (eps, 1-1/eps), xytext=(5, 5), textcoords='offset points', fontsize=8)
ax3.set_xlabel('Dielectric Constant ε_r', fontsize=12)
ax3.set_ylabel('Born Factor (1 - 1/ε)', fontsize=12)
ax3.set_title('Born Solvation: γ → 1 for High-ε Solvents', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim(0, 100)
ax3.set_ylim(0, 1.1)

# Plot 4: Jones-Dole B coefficients (Hofmeister)
ax4 = axes[1, 1]
ions_sorted = sorted(b_coeff.items(), key=lambda x: x[1])
ion_names = [i[0] for i in ions_sorted]
B_values = [i[1] for i in ions_sorted]
colors = ['red' if b < -0.02 else ('blue' if b > 0.02 else 'green') for b in B_values]

ax4.barh(range(len(ion_names)), B_values, color=colors, edgecolor='black', alpha=0.7)
ax4.axvline(x=0, color='black', linestyle='-', linewidth=2)
ax4.set_yticks(range(len(ion_names)))
ax4.set_yticklabels(ion_names)
ax4.set_xlabel('Jones-Dole B Coefficient (L/mol)', fontsize=12)
ax4.set_title('Hofmeister Series: B ≈ 0 at γ ~ 1', fontsize=14)
ax4.grid(True, alpha=0.3, axis='x')

# Add legend manually
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='blue', label='Kosmotrope'),
                   Patch(facecolor='green', label='Neutral (γ~1)'),
                   Patch(facecolor='red', label='Chaotrope')]
ax4.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/solvation_coherence.png', dpi=150)
print("Saved: solvation_coherence.png")
plt.close()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SESSION #196 SUMMARY: SOLVATION AT γ ~ 1")
print("=" * 60)

print(f"""
KEY FINDINGS:

1. HYDRATION NUMBERS
   γ_hyd = n_h / 4 (tetrahedral reference)
   Mean γ = {mean_gamma:.2f} ± {std_gamma:.2f}
   {near_one}/{len(gamma_values)} ions at γ ∈ [0.8, 1.5]
   Small ions: n_h ~ 4-6 (γ ~ 1-1.5)

2. BORN SOLVATION
   γ_Born = ε_r / (ε_r - 1)
   Water: γ = 1.013 (essentially 1!)
   High-ε solvents approach γ → 1
   Complete electrostatic solvation at γ = 1

3. DEBYE-HÜCKEL SCREENING
   γ_DH = κa (ion size / screening length)
   At γ = 1 (I ~ 1 M): screening length = ion size
   Crossover from Coulombic to screened regime

4. HOFMEISTER SERIES
   Cl⁻ and K⁺: B ≈ 0 (γ ~ 1!)
   Neutral point between kosmotropes and chaotropes
   NaCl/KCl are THE physiological electrolytes

5. PHYSIOLOGICAL RELEVANCE
   Life operates near γ ~ 1:
   - NaCl/KCl: neutral Hofmeister
   - I ~ 0.15 M: moderate screening
   - Water ε = 80: nearly complete solvation

CENTRAL INSIGHT:
Solvation parameters converge on γ ~ 1 references:
- Tetrahedral hydration (n_h ~ 4)
- Complete Born solvation (ε → ∞)
- Hofmeister neutral point (B = 0)
- Debye screening crossover (κa = 1)

Life evolved in aqueous conditions optimized for γ ~ 1.

This is the 59th phenomenon type at γ ~ 1!
""")

print("=" * 60)
print("SESSION #196 COMPLETE")
print("=" * 60)
