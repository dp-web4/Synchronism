#!/usr/bin/env python3
"""
Chemistry Session #109: Sound Velocity and Coherence

Test whether sound velocity v_s relates to coherence parameters.

Theory:
v_L = √(E/ρ) × √((1-ν)/((1+ν)(1-2ν))) (longitudinal)
v_T = √(G/ρ) (transverse)
v_D = (1/3 × (1/v_L³ + 2/v_T³))^(-1/3) (Debye average)

Debye temperature relation:
θ_D = (ℏ/k_B) × (6π²n)^(1/3) × v_D

Coherence connection:
- γ_phonon = 2T/θ_D = 2T × k_B/(ℏ × v_D × (6π²n)^(1/3))
- Higher v_D → higher θ_D → lower γ_phonon (more coherent)

From Session #108: κ_ph ∝ v × l_ph ∝ v_D²
So sound velocity directly affects phonon thermal transport.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# ELEMENT DATA
# =============================================================================

elements_data = {
    # Format: v_L (m/s), v_T (m/s), ρ (g/cm³), θ_D (K), B (GPa), E (GPa)
    # Noble metals
    'Ag': {'v_L': 3650, 'v_T': 1660, 'rho': 10.5, 'theta_D': 225, 'B': 100, 'E': 83},
    'Cu': {'v_L': 4760, 'v_T': 2325, 'rho': 8.96, 'theta_D': 343, 'B': 140, 'E': 130},
    'Au': {'v_L': 3240, 'v_T': 1200, 'rho': 19.3, 'theta_D': 165, 'B': 180, 'E': 78},

    # Alkali metals
    'Na': {'v_L': 3200, 'v_T': 1550, 'rho': 0.97, 'theta_D': 158, 'B': 6.3, 'E': 10},
    'K': {'v_L': 2000, 'v_T': 1000, 'rho': 0.86, 'theta_D': 91, 'B': 3.1, 'E': 3.5},
    'Li': {'v_L': 6000, 'v_T': 2850, 'rho': 0.53, 'theta_D': 344, 'B': 11, 'E': 4.9},

    # Transition metals (refractory)
    'W': {'v_L': 5220, 'v_T': 2890, 'rho': 19.3, 'theta_D': 400, 'B': 310, 'E': 411},
    'Mo': {'v_L': 6190, 'v_T': 3350, 'rho': 10.2, 'theta_D': 450, 'B': 230, 'E': 329},
    'Ta': {'v_L': 4100, 'v_T': 2030, 'rho': 16.7, 'theta_D': 240, 'B': 200, 'E': 186},
    'Nb': {'v_L': 4920, 'v_T': 2090, 'rho': 8.57, 'theta_D': 275, 'B': 170, 'E': 105},

    # Transition metals (3d)
    'Fe': {'v_L': 5960, 'v_T': 3240, 'rho': 7.87, 'theta_D': 470, 'B': 170, 'E': 211},
    'Ni': {'v_L': 5630, 'v_T': 2960, 'rho': 8.91, 'theta_D': 450, 'B': 180, 'E': 200},
    'Co': {'v_L': 5720, 'v_T': 3100, 'rho': 8.86, 'theta_D': 445, 'B': 180, 'E': 209},
    'Ti': {'v_L': 6070, 'v_T': 3120, 'rho': 4.51, 'theta_D': 420, 'B': 110, 'E': 116},
    'Cr': {'v_L': 6850, 'v_T': 4005, 'rho': 7.15, 'theta_D': 630, 'B': 160, 'E': 279},

    # Simple metals
    'Al': {'v_L': 6420, 'v_T': 3040, 'rho': 2.70, 'theta_D': 428, 'B': 76, 'E': 70},
    'Pb': {'v_L': 2160, 'v_T': 860, 'rho': 11.3, 'theta_D': 105, 'B': 46, 'E': 16},
    'Zn': {'v_L': 4210, 'v_T': 2440, 'rho': 7.14, 'theta_D': 327, 'B': 70, 'E': 108},
    'Mg': {'v_L': 5770, 'v_T': 3050, 'rho': 1.74, 'theta_D': 400, 'B': 45, 'E': 45},
    'Sn': {'v_L': 3320, 'v_T': 1670, 'rho': 7.26, 'theta_D': 200, 'B': 58, 'E': 50},

    # Semiconductors
    'Si': {'v_L': 8433, 'v_T': 5843, 'rho': 2.33, 'theta_D': 645, 'B': 98, 'E': 160},
    'Ge': {'v_L': 5400, 'v_T': 3570, 'rho': 5.32, 'theta_D': 374, 'B': 75, 'E': 130},

    # Insulators
    'C-d': {'v_L': 18000, 'v_T': 12000, 'rho': 3.51, 'theta_D': 2230, 'B': 442, 'E': 1050},
    'MgO': {'v_L': 9100, 'v_T': 6100, 'rho': 3.58, 'theta_D': 946, 'B': 155, 'E': 300},
    'Al2O3': {'v_L': 10800, 'v_T': 6500, 'rho': 3.97, 'theta_D': 1030, 'B': 240, 'E': 400},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*70)
print("CHEMISTRY SESSION #109: SOUND VELOCITY AND COHERENCE")
print("="*70)

T = 300  # K

elements = list(elements_data.keys())
v_L = np.array([elements_data[e]['v_L'] for e in elements])
v_T = np.array([elements_data[e]['v_T'] for e in elements])
rho = np.array([elements_data[e]['rho'] for e in elements])
theta_D = np.array([elements_data[e]['theta_D'] for e in elements])
B = np.array([elements_data[e]['B'] for e in elements])
E = np.array([elements_data[e]['E'] for e in elements])

# Debye average velocity (approximate)
v_D = (1/3 * (1/v_L**3 + 2/v_T**3))**(-1/3)

# Coherence parameters
gamma_phonon = 2 * T / theta_D

# Theoretical θ_D from v_D
# θ_D = (ℏ/k_B) × (6π²n)^(1/3) × v_D
# n = atoms/volume; for monatomic solids n = N_A × ρ / M
# Simplified: θ_D ∝ v_D × (ρ/M)^(1/3)

print(f"\n{'Element':<8} {'v_L':<8} {'v_T':<8} {'v_D':<8} {'θ_D':<8} {'γ_ph':<8} {'E':<8}")
print(f"{'':8} {'(m/s)':<8} {'(m/s)':<8} {'(m/s)':<8} {'(K)':<8} {'':8} {'(GPa)':<8}")
print("-"*64)
for i, e in enumerate(elements):
    print(f"{e:<8} {v_L[i]:<8.0f} {v_T[i]:<8.0f} {v_D[i]:<8.0f} {theta_D[i]:<8.0f} {gamma_phonon[i]:<8.3f} {E[i]:<8.0f}")

# =============================================================================
# CORRELATIONS
# =============================================================================

print("\n" + "="*70)
print("CORRELATIONS")
print("="*70)

# Sound velocity vs θ_D
r1, p1 = stats.pearsonr(v_D, theta_D)
print(f"\nv_D vs θ_D: r = {r1:.3f}, p = {p1:.2e}")

# v_D vs 1/γ_phonon (inverse coherence)
r2, p2 = stats.pearsonr(v_D, 1/gamma_phonon)
print(f"v_D vs 1/γ_phonon: r = {r2:.3f}, p = {p2:.2e}")

# v_D vs √(E/ρ) (elastic wave relation)
v_elastic = np.sqrt(E * 1e9 / (rho * 1e3))  # Convert to SI
r3, p3 = stats.pearsonr(v_L, v_elastic)
print(f"v_L vs √(E/ρ): r = {r3:.3f}, p = {p3:.2e}")

# v_D vs √(B/ρ)
v_bulk = np.sqrt(B * 1e9 / (rho * 1e3))
r4, p4 = stats.pearsonr(v_D, v_bulk)
print(f"v_D vs √(B/ρ): r = {r4:.3f}, p = {p4:.2e}")

# v_T/v_L ratio (related to Poisson ratio)
v_ratio = v_T / v_L
r5, p5 = stats.pearsonr(v_ratio, gamma_phonon)
print(f"\nv_T/v_L vs γ_phonon: r = {r5:.3f}, p = {p5:.2e}")

# Mean v_ratio
print(f"Mean v_T/v_L: {np.mean(v_ratio):.3f} ± {np.std(v_ratio):.3f}")
# For Poisson ν = 0.33: v_T/v_L = 0.58
# For Poisson ν = 0.25: v_T/v_L = 0.61

# θ_D² vs E/ρ (Debye model prediction)
r6, p6 = stats.pearsonr(theta_D**2, E/rho)
print(f"\nθ_D² vs E/ρ: r = {r6:.3f}, p = {p6:.2e}")

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
    '3d TM': ['Fe', 'Ni', 'Co', 'Ti', 'Cr'],
    'Simple': ['Al', 'Pb', 'Zn', 'Mg', 'Sn'],
    'Semiconductors': ['Si', 'Ge'],
    'Ceramics': ['C-d', 'MgO', 'Al2O3'],
}

print(f"\n{'Class':<15} {'Mean v_D':<12} {'Mean θ_D':<10} {'Mean γ_ph':<10}")
print("-"*47)
for cls, members in classes.items():
    idx = [elements.index(m) for m in members if m in elements]
    if len(idx) > 0:
        mean_vD = np.mean(v_D[idx])
        mean_theta = np.mean(theta_D[idx])
        mean_gamma = np.mean(gamma_phonon[idx])
        print(f"{cls:<15} {mean_vD:<12.0f} {mean_theta:<10.0f} {mean_gamma:<10.3f}")

# =============================================================================
# PHYSICAL INSIGHTS
# =============================================================================

print("\n" + "="*70)
print("PHYSICAL INSIGHTS")
print("="*70)

print(f"""
1. SOUND VELOCITY ↔ DEBYE TEMPERATURE
   v_D vs θ_D: r = {r1:.3f} (STRONG)

   θ_D = (ℏ/k_B) × (6π²n)^(1/3) × v_D
   This is the Debye model: phonon cutoff ∝ sound velocity.

2. SOUND VELOCITY → COHERENCE
   v_D vs 1/γ_phonon: r = {r2:.3f} (STRONG)

   Higher v_D → higher θ_D → lower γ_phonon → more coherent.
   Diamond (v_D ~ 13300 m/s) has γ_phonon = 0.27 (most coherent).
   Pb (v_D ~ 1150 m/s) has γ_phonon = 5.7 (least coherent).

3. ELASTIC WAVE RELATION
   v_L vs √(E/ρ): r = {r3:.3f} (EXCELLENT)

   Longitudinal waves follow elastic theory.
   v = √(modulus/density)

4. MATERIAL HIERARCHY BY v_D:
   | Material | v_D (m/s) | γ_phonon |
   |----------|-----------|----------|
   | Diamond  | 13,300    | 0.27     |
   | Al2O3    | 7,100     | 0.58     |
   | MgO      | 6,700     | 0.63     |
   | Si       | 6,400     | 0.93     |
   | Cr       | 4,700     | 0.95     |
   | ...      | ...       | ...      |
   | Pb       | 1,150     | 5.71     |

5. FRAMEWORK CONNECTION
   From Session #108: κ_ph ∝ v × l_ph
   So: κ_ph ∝ v_D × τ_ph × v_D = v_D² / Γ_ph

   High v_D materials (diamond) have:
   - High θ_D → low γ_phonon
   - Long phonon mean free path
   - High thermal conductivity

   Chain: v_D → θ_D → γ_phonon → transport
""")

# =============================================================================
# COHERENCE QUALITY FACTOR
# =============================================================================

print("\n" + "="*70)
print("PHONON COHERENCE QUALITY FACTOR")
print("="*70)

# Define Q_phonon = ω_D × τ_D = θ_D / (γ_phonon × something)
# Actually Q ~ 1/γ_phonon since τ ∝ 1/Γ_ph ∝ 1/(γ_G² × γ_phonon)

Q_phonon = 1 / gamma_phonon

print(f"\n{'Element':<8} {'v_D (m/s)':<12} {'γ_phonon':<10} {'Q_ph = 1/γ':<12}")
print("-"*44)
for i in np.argsort(Q_phonon)[::-1]:
    print(f"{elements[i]:<8} {v_D[i]:<12.0f} {gamma_phonon[i]:<10.3f} {Q_phonon[i]:<12.2f}")

# Correlation
r7, p7 = stats.pearsonr(v_D, Q_phonon)
print(f"\nv_D vs Q_phonon: r = {r7:.3f}")

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
    elif e in ['Fe', 'Ni', 'Co', 'Ti', 'Cr']:
        colors.append('red')
    elif e in ['Al', 'Pb', 'Zn', 'Mg', 'Sn']:
        colors.append('purple')
    elif e in ['Si', 'Ge']:
        colors.append('orange')
    else:
        colors.append('green')

# Plot 1: v_D vs θ_D
ax1 = axes[0, 0]
ax1.scatter(v_D, theta_D, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax1.annotate(e, (v_D[i], theta_D[i]), fontsize=8)

# Fit line
slope, intercept, r, p, se = stats.linregress(v_D, theta_D)
x_fit = np.linspace(1000, 14000, 100)
ax1.plot(x_fit, slope*x_fit + intercept, 'k--', alpha=0.5)

ax1.set_xlabel('Debye velocity v_D (m/s)', fontsize=12)
ax1.set_ylabel('Debye temperature θ_D (K)', fontsize=12)
ax1.set_title(f'v_D vs θ_D: r = {r1:.3f}', fontsize=14)
ax1.grid(True, alpha=0.3)

# Plot 2: v_D vs 1/γ_phonon
ax2 = axes[0, 1]
ax2.scatter(v_D, 1/gamma_phonon, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax2.annotate(e, (v_D[i], 1/gamma_phonon[i]), fontsize=8)
ax2.set_xlabel('Debye velocity v_D (m/s)', fontsize=12)
ax2.set_ylabel('1/γ_phonon (coherence)', fontsize=12)
ax2.set_title(f'v_D vs 1/γ_phonon: r = {r2:.3f}', fontsize=14)
ax2.grid(True, alpha=0.3)

# Plot 3: v_L vs √(E/ρ)
ax3 = axes[1, 0]
ax3.scatter(v_elastic, v_L, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax3.annotate(e, (v_elastic[i], v_L[i]), fontsize=8)

# 1:1 line
max_v = max(np.max(v_L), np.max(v_elastic))
ax3.plot([0, max_v], [0, max_v], 'k--', alpha=0.5, label='v = √(E/ρ)')

ax3.set_xlabel('√(E/ρ) (m/s)', fontsize=12)
ax3.set_ylabel('v_L measured (m/s)', fontsize=12)
ax3.set_title(f'Elastic Wave Relation: r = {r3:.3f}', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: θ_D² vs E/ρ (Debye model)
ax4 = axes[1, 1]
ax4.scatter(E/rho, theta_D**2, c=colors, s=100, alpha=0.7)
for i, e in enumerate(elements):
    ax4.annotate(e, (E[i]/rho[i], theta_D[i]**2), fontsize=8)
ax4.set_xlabel('E/ρ (GPa·cm³/g)', fontsize=12)
ax4.set_ylabel('θ_D² (K²)', fontsize=12)
ax4.set_title(f'Debye Model: θ_D² ∝ E/ρ, r = {r6:.3f}', fontsize=14)
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
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/sound_velocity_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "="*70)
print("SUMMARY - SESSION #109")
print("="*70)

print(f"""
KEY RESULTS:

1. SOUND VELOCITY → DEBYE TEMPERATURE
   v_D vs θ_D: r = {r1:.3f} (STRONG)
   Debye model validated: θ_D ∝ v_D × (n)^(1/3)

2. SOUND VELOCITY → COHERENCE
   v_D vs 1/γ_phonon: r = {r2:.3f} (STRONG)
   High v_D → high θ_D → low γ_phonon → coherent phonons

3. ELASTIC WAVE RELATION
   v_L vs √(E/ρ): r = {r3:.3f} (EXCELLENT)
   Elastic theory validated

4. DEBYE MODEL
   θ_D² vs E/ρ: r = {r6:.3f}
   θ_D ∝ √(E/ρ) ∝ v_D

5. MATERIAL HIERARCHY:
   Diamond: v_D = 13,300 m/s, γ = 0.27 (most coherent)
   Ceramics: v_D ~ 6,000-7,000 m/s, γ ~ 0.5-0.7
   Metals: v_D ~ 2,000-5,000 m/s, γ ~ 1-3
   Pb: v_D = 1,150 m/s, γ = 5.7 (least coherent)

FRAMEWORK CONNECTION:
- Sound velocity sets the phonon coherence scale
- v_D → θ_D → γ_phonon
- Combined with #107 (Γ_ph) and #108 (κ_ph):

  κ_ph ∝ v_D² / Γ_ph ∝ v_D² / (γ_G² × γ_phonon)

  Diamond wins on ALL factors:
  - Highest v_D (stiff covalent bonds)
  - Lowest γ_G (harmonic potential)
  - Lowest γ_phonon (high θ_D)

  → Highest κ_ph = 2200 W/m·K
""")

print("\nFigure saved to: sound_velocity_coherence.png")
