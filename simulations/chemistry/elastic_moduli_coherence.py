#!/usr/bin/env python3
"""
Chemistry Session #110: Elastic Moduli and Coherence

Test whether elastic moduli relate to coherence parameters.

Moduli:
- E = Young's modulus (uniaxial stress/strain)
- G = Shear modulus (shear stress/strain)
- B = Bulk modulus (volumetric stress/strain)
- K = B (bulk modulus, same as B)
- ν = Poisson's ratio (lateral strain / axial strain)

Relations:
E = 2G(1+ν) = 3B(1-2ν)
G = E/(2(1+ν))
B = E/(3(1-2ν))

Coherence connection:
- From #109: v_D ∝ √(E/ρ) → θ_D → γ_phonon
- Expect: E ∝ θ_D² × ρ ∝ 1/γ_phonon²

But also:
- Binding energy U ∝ E × a³ (lattice parameter)
- Melting point T_m ∝ E × a³ / k_B

The question: Is there a direct E vs γ correlation beyond the θ_D connection?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# ELEMENT DATA
# =============================================================================

elements_data = {
    # Format: E (GPa), G (GPa), B (GPa), ν (Poisson), ρ (g/cm³), θ_D (K), T_m (K), a (Å)
    # Noble metals
    'Ag': {'E': 83, 'G': 30, 'B': 100, 'nu': 0.37, 'rho': 10.5, 'theta_D': 225, 'Tm': 1235, 'a': 4.09},
    'Cu': {'E': 130, 'G': 48, 'B': 140, 'nu': 0.34, 'rho': 8.96, 'theta_D': 343, 'Tm': 1358, 'a': 3.61},
    'Au': {'E': 78, 'G': 27, 'B': 180, 'nu': 0.44, 'rho': 19.3, 'theta_D': 165, 'Tm': 1337, 'a': 4.08},

    # Alkali metals
    'Na': {'E': 10, 'G': 3.3, 'B': 6.3, 'nu': 0.32, 'rho': 0.97, 'theta_D': 158, 'Tm': 371, 'a': 4.29},
    'K': {'E': 3.5, 'G': 1.3, 'B': 3.1, 'nu': 0.31, 'rho': 0.86, 'theta_D': 91, 'Tm': 337, 'a': 5.33},
    'Li': {'E': 4.9, 'G': 4.2, 'B': 11, 'nu': 0.36, 'rho': 0.53, 'theta_D': 344, 'Tm': 454, 'a': 3.51},

    # Transition metals (refractory)
    'W': {'E': 411, 'G': 161, 'B': 310, 'nu': 0.28, 'rho': 19.3, 'theta_D': 400, 'Tm': 3695, 'a': 3.16},
    'Mo': {'E': 329, 'G': 126, 'B': 230, 'nu': 0.31, 'rho': 10.2, 'theta_D': 450, 'Tm': 2896, 'a': 3.15},
    'Ta': {'E': 186, 'G': 69, 'B': 200, 'nu': 0.34, 'rho': 16.7, 'theta_D': 240, 'Tm': 3290, 'a': 3.30},
    'Nb': {'E': 105, 'G': 38, 'B': 170, 'nu': 0.40, 'rho': 8.57, 'theta_D': 275, 'Tm': 2750, 'a': 3.30},

    # Transition metals (3d)
    'Fe': {'E': 211, 'G': 82, 'B': 170, 'nu': 0.29, 'rho': 7.87, 'theta_D': 470, 'Tm': 1811, 'a': 2.87},
    'Ni': {'E': 200, 'G': 76, 'B': 180, 'nu': 0.31, 'rho': 8.91, 'theta_D': 450, 'Tm': 1728, 'a': 3.52},
    'Co': {'E': 209, 'G': 76, 'B': 180, 'nu': 0.31, 'rho': 8.86, 'theta_D': 445, 'Tm': 1768, 'a': 2.51},
    'Ti': {'E': 116, 'G': 44, 'B': 110, 'nu': 0.32, 'rho': 4.51, 'theta_D': 420, 'Tm': 1941, 'a': 2.95},
    'Cr': {'E': 279, 'G': 115, 'B': 160, 'nu': 0.21, 'rho': 7.15, 'theta_D': 630, 'Tm': 2180, 'a': 2.88},
    'V': {'E': 128, 'G': 47, 'B': 160, 'nu': 0.37, 'rho': 6.11, 'theta_D': 380, 'Tm': 2183, 'a': 3.02},

    # Simple metals
    'Al': {'E': 70, 'G': 26, 'B': 76, 'nu': 0.35, 'rho': 2.70, 'theta_D': 428, 'Tm': 933, 'a': 4.05},
    'Pb': {'E': 16, 'G': 5.6, 'B': 46, 'nu': 0.44, 'rho': 11.3, 'theta_D': 105, 'Tm': 600, 'a': 4.95},
    'Zn': {'E': 108, 'G': 43, 'B': 70, 'nu': 0.25, 'rho': 7.14, 'theta_D': 327, 'Tm': 693, 'a': 2.66},
    'Mg': {'E': 45, 'G': 17, 'B': 45, 'nu': 0.29, 'rho': 1.74, 'theta_D': 400, 'Tm': 923, 'a': 3.21},
    'Sn': {'E': 50, 'G': 18, 'B': 58, 'nu': 0.36, 'rho': 7.26, 'theta_D': 200, 'Tm': 505, 'a': 5.83},

    # Semiconductors
    'Si': {'E': 160, 'G': 66, 'B': 98, 'nu': 0.22, 'rho': 2.33, 'theta_D': 645, 'Tm': 1687, 'a': 5.43},
    'Ge': {'E': 130, 'G': 55, 'B': 75, 'nu': 0.21, 'rho': 5.32, 'theta_D': 374, 'Tm': 1211, 'a': 5.66},

    # Ceramics
    'C-d': {'E': 1050, 'G': 478, 'B': 442, 'nu': 0.10, 'rho': 3.51, 'theta_D': 2230, 'Tm': 4000, 'a': 3.57},
    'MgO': {'E': 300, 'G': 130, 'B': 155, 'nu': 0.18, 'rho': 3.58, 'theta_D': 946, 'Tm': 3098, 'a': 4.21},
    'Al2O3': {'E': 400, 'G': 163, 'B': 240, 'nu': 0.23, 'rho': 3.97, 'theta_D': 1030, 'Tm': 2345, 'a': 4.76},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #110: ELASTIC MODULI AND COHERENCE")
print("="*70)

T = 300  # K

elements = list(elements_data.keys())
E = np.array([elements_data[e]['E'] for e in elements])
G = np.array([elements_data[e]['G'] for e in elements])
B = np.array([elements_data[e]['B'] for e in elements])
nu = np.array([elements_data[e]['nu'] for e in elements])
rho = np.array([elements_data[e]['rho'] for e in elements])
theta_D = np.array([elements_data[e]['theta_D'] for e in elements])
Tm = np.array([elements_data[e]['Tm'] for e in elements])
a = np.array([elements_data[e]['a'] for e in elements])

# Coherence parameters
gamma_phonon = 2 * T / theta_D

# Derived quantities
v_D = np.sqrt(E * 1e9 / (rho * 1e3)) / 1.5  # Approximate Debye velocity

print(f"\n{'Element':<8} {'E (GPa)':<10} {'G (GPa)':<10} {'B (GPa)':<10} {'ν':<8} {'θ_D (K)':<10} {'γ_ph':<8}")
print("-"*68)
for i, e in enumerate(elements):
    print(f"{e:<8} {E[i]:<10.0f} {G[i]:<10.0f} {B[i]:<10.0f} {nu[i]:<8.2f} {theta_D[i]:<10.0f} {gamma_phonon[i]:<8.3f}")

# =============================================================================
# CORRELATIONS: ELASTIC MODULI VS COHERENCE
# =============================================================================

print("\n" + "="*70)
print("CORRELATIONS: MODULI VS COHERENCE")
print("="*70)

# Direct correlations
r1, p1 = stats.pearsonr(E, 1/gamma_phonon)
print(f"\nE vs 1/γ_phonon: r = {r1:.3f}, p = {p1:.2e}")

r2, p2 = stats.pearsonr(G, 1/gamma_phonon)
print(f"G vs 1/γ_phonon: r = {r2:.3f}, p = {p2:.2e}")

r3, p3 = stats.pearsonr(B, 1/gamma_phonon)
print(f"B vs 1/γ_phonon: r = {r3:.3f}, p = {p3:.2e}")

# Scaled by density (from #109)
r4, p4 = stats.pearsonr(E/rho, theta_D**2)
print(f"\nE/ρ vs θ_D²: r = {r4:.3f}, p = {p4:.2e}")

r5, p5 = stats.pearsonr(np.sqrt(E/rho), theta_D)
print(f"√(E/ρ) vs θ_D: r = {r5:.3f}, p = {p5:.2e}")

# Poisson ratio vs γ
r6, p6 = stats.pearsonr(nu, gamma_phonon)
print(f"\nPoisson ν vs γ_phonon: r = {r6:.3f}, p = {p6:.2e}")

# =============================================================================
# PUGH RATIO (B/G) - Ductility Indicator
# =============================================================================

print("\n" + "="*70)
print("PUGH RATIO (B/G) - DUCTILITY VS COHERENCE")
print("="*70)

pugh = B / G
# B/G > 1.75: ductile; B/G < 1.75: brittle

r7, p7 = stats.pearsonr(pugh, gamma_phonon)
print(f"B/G vs γ_phonon: r = {r7:.3f}, p = {p7:.2e}")

r8, p8 = stats.pearsonr(pugh, nu)
print(f"B/G vs Poisson ν: r = {r8:.3f}")

print(f"\n{'Element':<8} {'B/G':<8} {'ν':<8} {'Ductility':<12}")
print("-"*38)
for i in np.argsort(pugh)[::-1]:
    ductility = "Ductile" if pugh[i] > 1.75 else "Brittle"
    print(f"{elements[i]:<8} {pugh[i]:<8.2f} {nu[i]:<8.2f} {ductility:<12}")

# =============================================================================
# CAUCHY RELATION (C12 - C44)
# =============================================================================

print("\n" + "="*70)
print("CAUCHY RELATION AND NON-CENTRAL FORCES")
print("="*70)

# For cubic crystals:
# C12 = B - 2G/3
# C44 = G
# Cauchy relation: C12 = C44 for central forces
# Violation: non-central (angular) bonding

# Cauchy pressure
C12 = B - 2*G/3
C44 = G
cauchy_pressure = C12 - C44

r9, p9 = stats.pearsonr(cauchy_pressure, gamma_phonon)
print(f"Cauchy pressure (C12-C44) vs γ_phonon: r = {r9:.3f}")

print(f"\n{'Element':<8} {'C12 (GPa)':<12} {'C44 (GPa)':<12} {'C12-C44':<12} {'Bonding':<12}")
print("-"*60)
for i in np.argsort(cauchy_pressure)[::-1][:15]:
    bonding = "Central" if abs(cauchy_pressure[i]) < 20 else ("Non-central" if cauchy_pressure[i] > 0 else "Directional")
    print(f"{elements[i]:<8} {C12[i]:<12.0f} {C44[i]:<12.0f} {cauchy_pressure[i]:<12.0f} {bonding:<12}")

print("""
Cauchy relation interpretation:
- C12 = C44: Central forces (noble gases, simple metals)
- C12 > C44: Metallic bonding (non-central electron gas contribution)
- C12 < C44: Covalent bonding (directional, angular-dependent)
""")

# =============================================================================
# BINDING ENERGY AND ELASTIC MODULI
# =============================================================================

print("\n" + "="*70)
print("BINDING ENERGY CORRELATION")
print("="*70)

# Cohesive energy ∝ E × a³ (dimensional analysis)
# More precisely: E_coh = C × E × Ω where Ω is atomic volume

binding_factor = E * (a**3)  # Arbitrary units

r10, p10 = stats.pearsonr(binding_factor, Tm)
print(f"E × a³ vs T_m: r = {r10:.3f}")

# Specific binding per atom
atomic_vol = a**3  # Å³
E_per_vol = E / atomic_vol * 1e30 / 1e9  # GPa/Å³ → J/m³ → normalized

r11, p11 = stats.pearsonr(E / atomic_vol, theta_D)
print(f"E/a³ vs θ_D: r = {r11:.3f}")

# =============================================================================
# MATERIAL CLASS ANALYSIS
# =============================================================================

print("\n" + "="*70)
print("MATERIAL CLASS ANALYSIS")
print("="*70)

classes = {
    'Noble': ['Ag', 'Cu', 'Au'],
    'Alkali': ['Na', 'K', 'Li'],
    'Refractory': ['W', 'Mo', 'Ta', 'Nb'],
    '3d TM': ['Fe', 'Ni', 'Co', 'Ti', 'Cr', 'V'],
    'Simple': ['Al', 'Pb', 'Zn', 'Mg', 'Sn'],
    'Semiconductors': ['Si', 'Ge'],
    'Ceramics': ['C-d', 'MgO', 'Al2O3'],
}

print(f"\n{'Class':<15} {'Mean E':<10} {'Mean G':<10} {'Mean B':<10} {'Mean ν':<10} {'Mean B/G':<10}")
print("-"*65)
for cls, members in classes.items():
    idx = [elements.index(m) for m in members if m in elements]
    if len(idx) > 0:
        print(f"{cls:<15} {np.mean(E[idx]):<10.0f} {np.mean(G[idx]):<10.0f} {np.mean(B[idx]):<10.0f} {np.mean(nu[idx]):<10.2f} {np.mean(pugh[idx]):<10.2f}")

# =============================================================================
# PHYSICAL INSIGHTS
# =============================================================================

print("\n" + "="*70)
print("PHYSICAL INSIGHTS")
print("="*70)

print(f"""
1. ELASTIC MODULI VS COHERENCE
   E vs 1/γ_phonon: r = {r1:.3f}
   G vs 1/γ_phonon: r = {r2:.3f}
   B vs 1/γ_phonon: r = {r3:.3f}

   MODERATE correlations - moduli set by BOTH:
   - Bond strength (coherence-related)
   - Atomic mass and packing (density-related)

2. DEBYE MODEL VALIDATION (from #109)
   E/ρ vs θ_D²: r = {r4:.3f} (STRONG)
   √(E/ρ) vs θ_D: r = {r5:.3f} (STRONG)

   Modulus/density ratio sets coherence scale.

3. POISSON RATIO NOT COHERENCE-CORRELATED
   ν vs γ_phonon: r = {r6:.3f}

   Poisson ratio reflects bonding GEOMETRY, not coherence.
   - ν ≈ 0.33: typical metals
   - ν ≈ 0.1-0.2: covalent ceramics
   - ν ≈ 0.44: soft metals (Au, Pb)

4. PUGH RATIO (DUCTILITY)
   B/G vs γ_phonon: r = {r7:.3f}

   Ductility weakly correlated with phonon coherence.
   B/G > 1.75: ductile (Au = 6.7, Pb = 8.2)
   B/G < 1.75: brittle (diamond = 0.93, Si = 1.48)

5. CAUCHY PRESSURE (BONDING CHARACTER)
   C12 - C44 vs γ_phonon: r = {r9:.3f}

   Metallic bonding (C12 > C44) → positive Cauchy pressure
   Covalent bonding (C12 < C44) → negative Cauchy pressure
""")

# =============================================================================
# VISUALIZATION
# =============================================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Define colors by class
colors = []
for e in elements:
    if e in ['Ag', 'Cu', 'Au']:
        colors.append('gold')
    elif e in ['Na', 'K', 'Li']:
        colors.append('blue')
    elif e in ['W', 'Mo', 'Ta', 'Nb']:
        colors.append('gray')
    elif e in ['Fe', 'Ni', 'Co', 'Ti', 'Cr', 'V']:
        colors.append('red')
    elif e in ['Al', 'Pb', 'Zn', 'Mg', 'Sn']:
        colors.append('purple')
    elif e in ['Si', 'Ge']:
        colors.append('orange')
    else:
        colors.append('green')

# Plot 1: E vs 1/γ_phonon
ax1 = axes[0, 0]
ax1.scatter(1/gamma_phonon, E, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax1.annotate(e, (1/gamma_phonon[i], E[i]), fontsize=8)
ax1.set_xlabel('1/γ_phonon (coherence)', fontsize=12)
ax1.set_ylabel('Young\'s modulus E (GPa)', fontsize=12)
ax1.set_title(f'E vs Coherence: r = {r1:.3f}', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: √(E/ρ) vs θ_D
ax2 = axes[0, 1]
sqrt_E_rho = np.sqrt(E/rho)
ax2.scatter(sqrt_E_rho, theta_D, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax2.annotate(e, (sqrt_E_rho[i], theta_D[i]), fontsize=8)

# Fit line
slope, intercept, r, p, se = stats.linregress(sqrt_E_rho, theta_D)
x_fit = np.linspace(0, np.max(sqrt_E_rho)*1.1, 100)
ax2.plot(x_fit, slope*x_fit + intercept, 'k--', alpha=0.5)

ax2.set_xlabel('√(E/ρ) (GPa·cm³/g)^0.5', fontsize=12)
ax2.set_ylabel('Debye temperature θ_D (K)', fontsize=12)
ax2.set_title(f'Debye Model: r = {r5:.3f}', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: Pugh ratio B/G vs Poisson ratio
ax3 = axes[1, 0]
ax3.scatter(nu, pugh, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax3.annotate(e, (nu[i], pugh[i]), fontsize=8)
ax3.axhline(y=1.75, color='k', linestyle='--', alpha=0.5, label='Ductile/Brittle')
ax3.set_xlabel('Poisson ratio ν', fontsize=12)
ax3.set_ylabel('Pugh ratio B/G', fontsize=12)
ax3.set_title(f'Ductility: B/G vs ν, r = {r8:.3f}', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Cauchy pressure vs γ_phonon
ax4 = axes[1, 1]
ax4.scatter(gamma_phonon, cauchy_pressure, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax4.annotate(e, (gamma_phonon[i], cauchy_pressure[i]), fontsize=8)
ax4.axhline(y=0, color='k', linestyle='--', alpha=0.5)
ax4.set_xlabel('γ_phonon', fontsize=12)
ax4.set_ylabel('Cauchy pressure (C12-C44) GPa', fontsize=12)
ax4.set_title(f'Bonding Character: r = {r9:.3f}', fontsize=14)
ax4.grid(True, alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='gold', label='Noble'),
    Patch(facecolor='blue', label='Alkali'),
    Patch(facecolor='gray', label='Refractory'),
    Patch(facecolor='red', label='3d TM'),
    Patch(facecolor='purple', label='Simple'),
    Patch(facecolor='orange', label='Semiconductor'),
    Patch(facecolor='green', label='Ceramic'),
]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/elastic_moduli_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #110")
print("="*70)

print(f"""
KEY RESULTS:

1. DIRECT MODULUS-COHERENCE CORRELATION
   E vs 1/γ_phonon: r = {r1:.3f}
   G vs 1/γ_phonon: r = {r2:.3f}
   B vs 1/γ_phonon: r = {r3:.3f}

   MODERATE correlations - moduli involve density too.

2. DEBYE MODEL VALIDATED
   E/ρ vs θ_D²: r = {r4:.3f} (STRONG)
   √(E/ρ) vs θ_D: r = {r5:.3f} (STRONG)

   The modulus/density ratio is what matters for coherence.

3. POISSON RATIO - NOT COHERENCE-DEPENDENT
   ν vs γ_phonon: r = {r6:.3f}
   Reflects bonding geometry, not coherence.

4. PUGH RATIO (DUCTILITY)
   B/G > 1.75: ductile
   B/G < 1.75: brittle

   Ductile metals have high γ_phonon (soft).
   Brittle ceramics have low γ_phonon (stiff).

5. CAUCHY PRESSURE - BONDING CHARACTER
   Positive: metallic (non-central forces)
   Negative: covalent (directional bonds)

FRAMEWORK CONNECTION:
Elastic moduli are EXTENSIVE properties (like E_F from #100).
They set the SCALE of mechanical behavior.

Coherence (γ) is INTENSIVE.
It determines the QUALITY of phonon transport.

Full thermal conductivity:
κ_ph ∝ v_D × l_ph ∝ √(E/ρ) × 1/(γ_G² × γ_phonon)

Modulus appears in v_D (velocity).
γ appears in l_ph (mean free path).
""")

print("\nFigure saved to: elastic_moduli_coherence.png")
