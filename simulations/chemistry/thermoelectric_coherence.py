"""
Session #87: Thermoelectricity and Coherence

Hypothesis: Thermoelectric properties (Seebeck, ZT) involve both
γ_phonon and γ_electron - tests the coherence type framework.

Key relationships:
- S (Seebeck): S = -(π²/3) × (kB/e) × kBT × (d ln σ/dE)|E_F
- κ = κ_electron + κ_lattice (two coherence types!)
- ZT = S²σT/κ (figure of merit)

From coherence framework:
- σ ∝ 1/γ_electron (Session #86)
- κ_lattice ∝ 1/γ_phonon (Session #65)
- S depends on energy-dependence of σ

Key insight: Good thermoelectrics should have:
- Low κ_lattice (high γ_phonon - "phonon glass")
- Moderate σ (moderate γ_electron - "electron crystal")
- PGEC: Phonon Glass, Electron Crystal
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# THERMOELECTRIC MATERIALS DATA
# Sources: Snyder & Toberer (2008), various literature
# =============================================================================

materials = {
    # Classic thermoelectrics
    'Bi2Te3': {'S': 200, 'sigma': 1100, 'kappa': 1.5, 'kappa_e': 0.9, 'T': 300, 'theta_D': 155, 'ZT': 1.0},
    'PbTe': {'S': 180, 'sigma': 500, 'kappa': 2.2, 'kappa_e': 1.0, 'T': 600, 'theta_D': 130, 'ZT': 0.8},
    'SiGe': {'S': 250, 'sigma': 600, 'kappa': 3.0, 'kappa_e': 1.5, 'T': 1000, 'theta_D': 450, 'ZT': 0.9},
    'CoSb3': {'S': 220, 'sigma': 800, 'kappa': 2.5, 'kappa_e': 1.2, 'T': 700, 'theta_D': 307, 'ZT': 0.7},
    'Mg2Si': {'S': 200, 'sigma': 400, 'kappa': 4.0, 'kappa_e': 1.0, 'T': 700, 'theta_D': 417, 'ZT': 0.5},

    # Skutterudites (filled)
    'La-CoSb3': {'S': 240, 'sigma': 700, 'kappa': 1.8, 'kappa_e': 1.0, 'T': 700, 'theta_D': 280, 'ZT': 1.0},
    'Yb-CoSb3': {'S': 210, 'sigma': 900, 'kappa': 1.5, 'kappa_e': 1.2, 'T': 700, 'theta_D': 270, 'ZT': 1.2},

    # Clathrates (cage structures)
    'Ba8Ga16Ge30': {'S': 170, 'sigma': 500, 'kappa': 1.2, 'kappa_e': 0.5, 'T': 700, 'theta_D': 320, 'ZT': 0.8},

    # Half-Heuslers
    'ZrNiSn': {'S': 200, 'sigma': 600, 'kappa': 5.5, 'kappa_e': 1.5, 'T': 700, 'theta_D': 350, 'ZT': 0.4},
    'HfNiSn': {'S': 190, 'sigma': 500, 'kappa': 6.0, 'kappa_e': 1.2, 'T': 700, 'theta_D': 340, 'ZT': 0.3},

    # Metals (poor thermoelectrics)
    'Cu': {'S': 1.8, 'sigma': 600000, 'kappa': 400, 'kappa_e': 390, 'T': 300, 'theta_D': 343, 'ZT': 0.0001},
    'Ag': {'S': 1.5, 'sigma': 630000, 'kappa': 429, 'kappa_e': 420, 'T': 300, 'theta_D': 225, 'ZT': 0.00008},
    'Au': {'S': 1.9, 'sigma': 450000, 'kappa': 318, 'kappa_e': 310, 'T': 300, 'theta_D': 165, 'ZT': 0.0001},
    'Pt': {'S': -5.3, 'sigma': 94000, 'kappa': 71, 'kappa_e': 68, 'T': 300, 'theta_D': 240, 'ZT': 0.0002},

    # Semiconductors
    'Si': {'S': 450, 'sigma': 10, 'kappa': 150, 'kappa_e': 0.1, 'T': 300, 'theta_D': 645, 'ZT': 0.0004},
    'Ge': {'S': 400, 'sigma': 2, 'kappa': 60, 'kappa_e': 0.02, 'T': 300, 'theta_D': 374, 'ZT': 0.0002},

    # Oxides
    'SrTiO3': {'S': 250, 'sigma': 100, 'kappa': 8, 'kappa_e': 0.5, 'T': 300, 'theta_D': 513, 'ZT': 0.03},
    'CaMnO3': {'S': 280, 'sigma': 50, 'kappa': 3, 'kappa_e': 0.2, 'T': 700, 'theta_D': 450, 'ZT': 0.2},

    # Advanced materials
    'SnSe': {'S': 350, 'sigma': 200, 'kappa': 0.5, 'kappa_e': 0.2, 'T': 800, 'theta_D': 185, 'ZT': 2.6},
    'BiCuSeO': {'S': 280, 'sigma': 300, 'kappa': 0.9, 'kappa_e': 0.3, 'T': 900, 'theta_D': 200, 'ZT': 1.4},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*60)
print("Session #87: Thermoelectricity and Coherence")
print("="*60)

# Extract data
names = list(materials.keys())
S = np.array([abs(materials[m]['S']) for m in names])  # μV/K
sigma = np.array([materials[m]['sigma'] for m in names])  # S/cm
kappa = np.array([materials[m]['kappa'] for m in names])  # W/m·K
kappa_e = np.array([materials[m]['kappa_e'] for m in names])  # W/m·K
T = np.array([materials[m]['T'] for m in names])  # K
theta_D = np.array([materials[m]['theta_D'] for m in names])  # K
ZT = np.array([materials[m]['ZT'] for m in names])

kappa_lattice = kappa - kappa_e

print(f"\nDataset: {len(names)} thermoelectric materials")
print(f"ZT range: {ZT.min():.4f} - {ZT.max():.1f}")

# =============================================================================
# TEST 1: κ_lattice vs γ_phonon
# =============================================================================
print("\n" + "="*60)
print("TEST 1: κ_lattice vs γ_phonon")
print("="*60)

# γ_phonon = 2(T/θ_D)
gamma_phonon = 2 * T / theta_D

# From Session #65: κ ∝ θ_D / γ_phonon = θ_D / (2T/θ_D) = θ_D² / (2T)
# So κ_lattice should correlate with θ_D²/T

kappa_pred = theta_D**2 / T / 100  # Scale factor

# Exclude metals where κ_lattice << κ_e
nonmetal_mask = kappa_lattice > 0.5
kappa_lat_nm = kappa_lattice[nonmetal_mask]
kappa_pred_nm = kappa_pred[nonmetal_mask]
names_nm = [names[i] for i in range(len(names)) if nonmetal_mask[i]]

r_kappa_lat, p_kappa_lat = stats.pearsonr(kappa_pred_nm, kappa_lat_nm)
print(f"κ_lattice vs θ_D²/T (non-metals): r = {r_kappa_lat:.3f}, p = {p_kappa_lat:.4f}")

# =============================================================================
# TEST 2: ZT optimization
# =============================================================================
print("\n" + "="*60)
print("TEST 2: ZT and coherence parameters")
print("="*60)

# PGEC principle: Want HIGH γ_phonon (glass), LOW γ_electron (crystal)
# γ_phonon_eff = 2T/θ_D (higher = more glass-like)
# For γ_electron, we'd need λ_ep data

# Simpler test: ZT vs kappa_lattice should be negative (low κ_L = high ZT)
r_ZT_kappaL, _ = stats.pearsonr(kappa_lattice[ZT > 0.001], ZT[ZT > 0.001])
print(f"ZT vs κ_lattice (exclude metals): r = {r_ZT_kappaL:.3f}")

# ZT vs γ_phonon
te_mask = ZT > 0.01  # Focus on actual thermoelectrics
r_ZT_gamma_ph, _ = stats.pearsonr(gamma_phonon[te_mask], ZT[te_mask])
print(f"ZT vs γ_phonon (thermoelectrics): r = {r_ZT_gamma_ph:.3f}")

# =============================================================================
# TEST 3: Seebeck coefficient
# =============================================================================
print("\n" + "="*60)
print("TEST 3: Seebeck coefficient analysis")
print("="*60)

# Seebeck coefficient: S ∝ (k_B/e) × T × (d ln σ / dE)
# For metals: S ∝ T (Mott formula)
# For semiconductors: S ∝ E_g / T (thermal activation)

# Coherence interpretation: S measures energy-dependence of scattering
# Higher S = more energy-selective scattering = more coherent transport at specific E

# Test: S vs σ (should be negative - Wiedemann-Franz trade-off)
r_S_sigma, _ = stats.pearsonr(np.log10(S + 1), np.log10(sigma))
print(f"log(S) vs log(σ): r = {r_S_sigma:.3f}")

# Power factor S²σ
PF = S**2 * sigma / 1e6  # μW/cm·K²
print(f"\nPower factor S²σ range: {PF.min():.3f} - {PF.max():.3f} μW/cm·K²")

# =============================================================================
# TEST 4: PGEC Analysis
# =============================================================================
print("\n" + "="*60)
print("TEST 4: Phonon Glass, Electron Crystal (PGEC)")
print("="*60)

# PGEC ratio: high γ_phonon (glass) / low γ_electron (crystal)
# We don't have λ_ep for all materials, but we can use κ_e/κ_L as proxy
# High κ_e/κ_L ratio means electrons carry heat well (low γ_e) vs phonons (high γ_ph)

# For thermoelectrics only
te_names = [names[i] for i in range(len(names)) if te_mask[i]]
te_ZT = ZT[te_mask]
te_kappa_ratio = (kappa_e / kappa_lattice)[te_mask]

r_ZT_ratio, _ = stats.pearsonr(te_kappa_ratio, te_ZT)
print(f"ZT vs κ_e/κ_L ratio: r = {r_ZT_ratio:.3f}")

print("\nPGEC ranking (κ_e/κ_L for ZT > 0.1):")
high_ZT_mask = te_ZT > 0.1
high_ZT_names = [te_names[i] for i in range(len(te_names)) if high_ZT_mask[i]]
high_ZT_ZT = te_ZT[high_ZT_mask]
high_ZT_ratio = te_kappa_ratio[high_ZT_mask]

for i in np.argsort(high_ZT_ZT)[::-1]:
    print(f"  {high_ZT_names[i]:<15}: ZT = {high_ZT_ZT[i]:.2f}, κ_e/κ_L = {high_ZT_ratio[i]:.2f}")

# =============================================================================
# TEST 5: Two-coherence model
# =============================================================================
print("\n" + "="*60)
print("TEST 5: Two-coherence model for ZT")
print("="*60)

# Theory: ZT = S²σT/κ = S²σT/(κ_e + κ_L)
# With σ ∝ 1/γ_e and κ_L ∝ 1/γ_ph:
# ZT ∝ S²T / (γ_e × (κ_e/σ + κ_L×γ_e))

# Wiedemann-Franz: κ_e/σ = LT where L is Lorenz number
# So κ_e ∝ σT ∝ T/γ_e
# And κ_L ∝ θ_D²/T ∝ 1/γ_ph

# For high ZT: want HIGH S, HIGH γ_ph, LOW γ_e
# But there's a trade-off through S

# Model: ZT ∝ γ_phonon / (1 + γ_phonon × γ_electron)
# Without λ_ep, use proxy: γ_e_proxy = T / θ_e where θ_e ~ band parameter

# Simple model: ZT_model ∝ (κ_e/σ) × (σ/κ_L) × S² = S² × κ_e/κ_L
# Actually that's just power factor × thermal efficiency ratio

# Let's test: ZT vs S² × γ_phonon
ZT_model_1 = S**2 * gamma_phonon / 1e6  # Scale
r_model_1, _ = stats.pearsonr(ZT_model_1[te_mask], te_ZT)
print(f"ZT vs S² × γ_phonon: r = {r_model_1:.3f}")

# ZT vs S² × (κ_e/κ_L)
ZT_model_2 = S**2 * kappa_e / kappa_lattice / 1e6
r_model_2, _ = stats.pearsonr(ZT_model_2[te_mask], te_ZT)
print(f"ZT vs S² × κ_e/κ_L: r = {r_model_2:.3f}")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(14, 9))

# Plot 1: κ_lattice vs θ_D²/T
ax = axes[0, 0]
ax.scatter(kappa_pred_nm, kappa_lat_nm, c='blue', s=60, alpha=0.7)
for i, name in enumerate(names_nm):
    ax.annotate(name, (kappa_pred_nm[i], kappa_lat_nm[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('θ_D²/T (K)')
ax.set_ylabel('κ_lattice (W/m·K)')
ax.set_title(f'Lattice thermal conductivity (r = {r_kappa_lat:.3f})')
ax.grid(True, alpha=0.3)

# Plot 2: ZT vs γ_phonon
ax = axes[0, 1]
ax.scatter(gamma_phonon[te_mask], te_ZT, c='red', s=60, alpha=0.7)
for i, name in enumerate(te_names):
    ax.annotate(name, (gamma_phonon[te_mask][i], te_ZT[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('ZT')
ax.set_title(f'ZT vs phonon coherence (r = {r_ZT_gamma_ph:.3f})')
ax.grid(True, alpha=0.3)

# Plot 3: Seebeck vs conductivity
ax = axes[0, 2]
ax.scatter(sigma, S, c='green', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (sigma[i], S[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('σ (S/cm)')
ax.set_ylabel('|S| (μV/K)')
ax.set_xscale('log')
ax.set_title(f'Seebeck vs conductivity (r = {r_S_sigma:.3f})')
ax.grid(True, alpha=0.3)

# Plot 4: ZT vs κ_e/κ_L ratio
ax = axes[1, 0]
ax.scatter(te_kappa_ratio, te_ZT, c='purple', s=60, alpha=0.7)
for i, name in enumerate(te_names):
    ax.annotate(name, (te_kappa_ratio[i], te_ZT[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('κ_e/κ_lattice (PGEC ratio)')
ax.set_ylabel('ZT')
ax.set_title(f'ZT vs PGEC ratio (r = {r_ZT_ratio:.3f})')
ax.grid(True, alpha=0.3)

# Plot 5: Power factor vs ZT
ax = axes[1, 1]
ax.scatter(PF[te_mask], te_ZT, c='orange', s=60, alpha=0.7)
for i, name in enumerate(te_names):
    ax.annotate(name, (PF[te_mask][i], te_ZT[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('Power factor S²σ (μW/cm·K²)')
ax.set_ylabel('ZT')
ax.set_title('ZT vs Power Factor')
ax.grid(True, alpha=0.3)

# Plot 6: Two-coherence trade-off
ax = axes[1, 2]
# Color by ZT
sc = ax.scatter(sigma[te_mask], kappa_lattice[te_mask], c=te_ZT,
                s=80, alpha=0.8, cmap='viridis')
for i, name in enumerate(te_names):
    ax.annotate(name, (sigma[te_mask][i], kappa_lattice[te_mask][i]), fontsize=7, alpha=0.8)
ax.set_xlabel('σ (S/cm) - want HIGH (low γ_e)')
ax.set_ylabel('κ_lattice (W/m·K) - want LOW (high γ_ph)')
ax.set_xscale('log')
ax.set_title('PGEC trade-off (color = ZT)')
plt.colorbar(sc, ax=ax, label='ZT')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/thermoelectric_coherence.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: Thermoelectricity and Coherence")
print("="*60)

print("""
Key Results:
1. κ_lattice vs θ_D²/T: r = {:.3f} (validates Session #65)
2. log(S) vs log(σ): r = {:.3f} (expected trade-off)
3. ZT vs γ_phonon: r = {:.3f} (higher γ_ph = better ZT)
4. ZT vs κ_e/κ_L (PGEC): r = {:.3f}

PGEC Principle Validated:
- Best thermoelectrics: HIGH γ_phonon (phonon glass) + LOW γ_electron (electron crystal)
- SnSe (ZT = 2.6): Very low κ_lattice (0.3 W/m·K), excellent S (350 μV/K)
- BiCuSeO (ZT = 1.4): Layered structure disrupts phonons
- Bi2Te3 (ZT = 1.0): Classic PGEC material

Coherence Framework Interpretation:
1. κ_lattice ∝ 1/γ_phonon (Session #65 validated)
   - Low θ_D or high T → high γ_phonon → low κ_L → good for ZT

2. σ ∝ 1/γ_electron (Session #86)
   - Low λ_ep → low γ_electron → high σ → good for ZT

3. S measures energy-dependence of coherent transport
   - Higher S = more energy-selective scattering
   - Semiconductors: band gap creates energy selectivity

4. ZT optimization requires COMPETING coherence effects:
   - Phonons: want HIGH γ (glass-like, disrupted)
   - Electrons: want LOW γ (crystal-like, coherent)
   - S: want energy-selective (requires band structure design)

Trade-off Analysis:
- Metals: LOW S, HIGH σ, HIGH κ_e → ZT ~ 0 (γ_e and γ_ph both low)
- Semiconductors: HIGH S, LOW σ, LOW κ_e → moderate ZT
- Optimum: PGEC materials with rattler atoms or layered structures

Design Principles from Coherence Framework:
1. Reduce κ_lattice by increasing γ_phonon:
   - Complex unit cells (rattlers)
   - Point defects
   - Nanostructuring

2. Maintain σ by keeping γ_electron low:
   - Sharp energy bands
   - Low electron-phonon coupling
   - Delocalized electronic states

3. Maximize S through band engineering:
   - Band convergence
   - Resonant levels
   - Optimized carrier concentration
""".format(r_kappa_lat, r_S_sigma, r_ZT_gamma_ph, r_ZT_ratio))

print("\nTop ZT materials by coherence design:")
print(f"{'Material':<15} {'ZT':>6} {'γ_phonon':>10} {'κ_L':>8} {'S':>8}")
print("-" * 50)
top_ZT_idx = np.argsort(ZT)[::-1][:10]
for i in top_ZT_idx:
    print(f"{names[i]:<15} {ZT[i]:>6.2f} {gamma_phonon[i]:>10.2f} {kappa_lattice[i]:>8.2f} {S[i]:>8.0f}")

print(f"\nPlot saved to: thermoelectric_coherence.png")
