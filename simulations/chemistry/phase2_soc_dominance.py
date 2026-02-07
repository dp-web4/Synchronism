#!/usr/bin/env python3
"""
Phase 2 Session #4: Spin-Orbit Coupling Dominance

The coherence framework fails for magnetic properties of heavy elements
because spin-orbit coupling (SOC) dominates. SOC is an ATOMIC property
(∝ Z⁴) that overwhelms the COLLECTIVE coherence effect.

Questions:
1. Can we quantify WHEN SOC dominates over coherence?
2. Is there a dominance parameter D that predicts failure?
3. For light elements (3d), does coherence still contribute?
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

print("=" * 70)
print("PHASE 2 SESSION #4: SPIN-ORBIT COUPLING DOMINANCE")
print("=" * 70)

# ==============================================================================
# Dataset: Magnetic materials with atomic number information
# ==============================================================================

# From Sessions #94, #99: magnetostriction and anisotropy data
# Extended with SOC estimates from atomic spectroscopy

magnetic_data = {
    # 3d transition metals
    'Fe':    {'Z': 26, 'K1_MJ': 0.048,  'lambda_s': -7e-6,   'theta_D': 470, 'Tc': 1043,
              'SOC_meV': 50,   'class': '3d'},
    'Co':    {'Z': 27, 'K1_MJ': 0.45,   'lambda_s': -60e-6,  'theta_D': 445, 'Tc': 1388,
              'SOC_meV': 65,   'class': '3d'},
    'Ni':    {'Z': 28, 'K1_MJ': 0.005,  'lambda_s': -34e-6,  'theta_D': 450, 'Tc': 627,
              'SOC_meV': 80,   'class': '3d'},
    'Cr':    {'Z': 24, 'K1_MJ': 0.0001, 'lambda_s': 0,       'theta_D': 630, 'Tc': 311,
              'SOC_meV': 35,   'class': '3d'},

    # Rare earth metals — SOC is massive
    'Gd':   {'Z': 64, 'K1_MJ': 0.012,  'lambda_s': 0,       'theta_D': 182, 'Tc': 292,
              'SOC_meV': 200,  'class': 'RE'},  # Gd: L=0, so actually LOW anisotropy
    'Tb':   {'Z': 65, 'K1_MJ': 8.3,    'lambda_s': 8800e-6, 'theta_D': 177, 'Tc': 219,
              'SOC_meV': 300,  'class': 'RE'},
    'Dy':   {'Z': 66, 'K1_MJ': 4.5,    'lambda_s': 6000e-6, 'theta_D': 186, 'Tc': 85,
              'SOC_meV': 350,  'class': 'RE'},
    'Ho':   {'Z': 67, 'K1_MJ': 3.1,    'lambda_s': 4000e-6, 'theta_D': 190, 'Tc': 20,
              'SOC_meV': 400,  'class': 'RE'},
    'Er':   {'Z': 68, 'K1_MJ': 1.3,    'lambda_s': 2000e-6, 'theta_D': 195, 'Tc': 19,
              'SOC_meV': 450,  'class': 'RE'},

    # RE-TM compounds
    'SmCo5':    {'Z': 62, 'K1_MJ': 17.0,  'lambda_s': 1500e-6, 'theta_D': 400, 'Tc': 1020,
                  'SOC_meV': 280, 'class': 'RE-TM'},
    'Nd2Fe14B': {'Z': 60, 'K1_MJ': 4.9,   'lambda_s': 900e-6,  'theta_D': 350, 'Tc': 585,
                  'SOC_meV': 220, 'class': 'RE-TM'},

    # L10 ordered alloys (5d elements)
    'FePt':  {'Z': 78, 'K1_MJ': 6.6,   'lambda_s': 500e-6,  'theta_D': 300, 'Tc': 750,
              'SOC_meV': 600,  'class': '5d'},
    'CoPt':  {'Z': 78, 'K1_MJ': 4.9,   'lambda_s': 400e-6,  'theta_D': 290, 'Tc': 840,
              'SOC_meV': 600,  'class': '5d'},

    # Ferrites (low SOC, oxide)
    'Fe3O4':     {'Z': 26, 'K1_MJ': 0.011,  'lambda_s': 20e-6,   'theta_D': 600, 'Tc': 858,
                   'SOC_meV': 50,   'class': 'ferrite'},
    'CoFe2O4':   {'Z': 27, 'K1_MJ': 0.27,   'lambda_s': 250e-6,  'theta_D': 550, 'Tc': 793,
                   'SOC_meV': 65,   'class': 'ferrite'},
    'YIG':       {'Z': 26, 'K1_MJ': 0.0006, 'lambda_s': 2e-6,    'theta_D': 500, 'Tc': 560,
                   'SOC_meV': 50,   'class': 'ferrite'},
}

T = 300
names = list(magnetic_data.keys())
Z = np.array([magnetic_data[m]['Z'] for m in names])
K1 = np.array([magnetic_data[m]['K1_MJ'] for m in names])
lambda_s = np.array([abs(magnetic_data[m]['lambda_s']) for m in names])
theta_D = np.array([magnetic_data[m]['theta_D'] for m in names])
Tc = np.array([magnetic_data[m]['Tc'] for m in names])
SOC = np.array([magnetic_data[m]['SOC_meV'] for m in names])
classes = [magnetic_data[m]['class'] for m in names]

gamma_phonon = 2 * T / theta_D
gamma_spin = 2 * T / Tc

print(f"\n{'Material':<12} {'Z':<5} {'|K1|(MJ)':<10} {'SOC(meV)':<10} {'γ_ph':<8} {'γ_sp':<8} {'Class'}")
print("-" * 70)
for i, m in enumerate(names):
    print(f"{m:<12} {Z[i]:<5} {K1[i]:<10.3f} {SOC[i]:<10.0f} {gamma_phonon[i]:<8.2f} {gamma_spin[i]:<8.2f} {classes[i]}")

# ==============================================================================
# Analysis 1: K1 vs γ_phonon — does coherence predict anisotropy?
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 1: K1 vs γ_phonon (coherence prediction)")
print("=" * 70)

log_K1 = np.log10(K1 + 1e-6)  # avoid log(0)

r_K_gamma, p_K_gamma = stats.pearsonr(gamma_phonon, log_K1)
print(f"\nAll materials (N={len(names)}):")
print(f"  log|K1| vs γ_phonon: r = {r_K_gamma:.3f}  (p = {p_K_gamma:.3e})")

# By class
for cls in ['3d', 'RE', 'RE-TM', '5d', 'ferrite']:
    mask = np.array([c == cls for c in classes])
    if np.sum(mask) >= 3:
        r, p = stats.pearsonr(gamma_phonon[mask], log_K1[mask])
        print(f"  {cls} only (N={np.sum(mask)}): r = {r:.3f}  (p = {p:.3e})")

# ==============================================================================
# Analysis 2: K1 vs SOC — does atomic structure predict anisotropy?
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 2: K1 vs SOC (atomic prediction)")
print("=" * 70)

log_SOC = np.log10(SOC)

r_K_SOC, p_K_SOC = stats.pearsonr(log_SOC, log_K1)
print(f"\nAll materials (N={len(names)}):")
print(f"  log|K1| vs log(SOC): r = {r_K_SOC:.3f}  (p = {p_K_SOC:.3e})")

# K1 vs Z^4 (SOC ∝ Z^4 scaling)
r_K_Z4, p_K_Z4 = stats.pearsonr(np.log10(Z**4), log_K1)
print(f"  log|K1| vs log(Z⁴):  r = {r_K_Z4:.3f}  (p = {p_K_Z4:.3e})")

# K1 vs SOC^2 (second-order perturbation theory: K ∝ SOC²/W)
r_K_SOC2, p_K_SOC2 = stats.pearsonr(np.log10(SOC**2), log_K1)
print(f"  log|K1| vs log(SOC²): r = {r_K_SOC2:.3f}  (p = {p_K_SOC2:.3e})")

# ==============================================================================
# Analysis 3: Dominance parameter D
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 3: DOMINANCE PARAMETER D = SOC / k_B × θ_D")
print("=" * 70)

# D = SOC energy / phonon energy = ξ_SOC / (k_B × θ_D)
# When D >> 1: SOC dominates → γ_phonon prediction fails
# When D << 1: coherence dominates → γ_phonon prediction works
k_B_meV = 0.08617  # meV/K
D = SOC / (k_B_meV * theta_D)

print(f"\n{'Material':<12} {'SOC(meV)':<10} {'θ_D(K)':<8} {'k_Bθ_D(meV)':<12} {'D=SOC/k_Bθ_D':<14} {'Class'}")
print("-" * 70)
for i, m in enumerate(names):
    phonon_energy = k_B_meV * theta_D[i]
    print(f"{m:<12} {SOC[i]:<10.0f} {theta_D[i]:<8.0f} {phonon_energy:<12.1f} {D[i]:<14.1f} {classes[i]}")

# Classify by dominance
D_threshold = 5.0
soc_dominated = D > D_threshold
coh_relevant = D <= D_threshold

print(f"\nDominance classification (threshold D = {D_threshold}):")
print(f"  SOC-dominated (D > {D_threshold}): {sum(soc_dominated)} materials")
print(f"  Coherence-relevant (D ≤ {D_threshold}): {sum(coh_relevant)} materials")

print(f"\nSOC-dominated materials:")
for i, m in enumerate(names):
    if soc_dominated[i]:
        print(f"  {m}: D = {D[i]:.1f}")

print(f"\nCoherence-relevant materials:")
for i, m in enumerate(names):
    if coh_relevant[i]:
        print(f"  {m}: D = {D[i]:.1f}")

# ==============================================================================
# Analysis 4: Can γ_phonon predict K1 within the coherence-relevant subset?
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 4: γ_phonon PREDICTION WITHIN D ≤ {0} MATERIALS".format(D_threshold))
print("=" * 70)

if np.sum(coh_relevant) >= 3:
    r_coh, p_coh = stats.pearsonr(gamma_phonon[coh_relevant], log_K1[coh_relevant])
    print(f"\nCoherence-relevant (N={np.sum(coh_relevant)}):")
    print(f"  log|K1| vs γ_phonon: r = {r_coh:.3f}  (p = {p_coh:.3e})")
else:
    r_coh = 0
    print("\nToo few coherence-relevant materials for correlation")

if np.sum(soc_dominated) >= 3:
    r_soc, p_soc = stats.pearsonr(gamma_phonon[soc_dominated], log_K1[soc_dominated])
    print(f"\nSOC-dominated (N={np.sum(soc_dominated)}):")
    print(f"  log|K1| vs γ_phonon: r = {r_soc:.3f}  (p = {p_soc:.3e})")
else:
    r_soc = 0

# ==============================================================================
# Analysis 5: The Gd anomaly — RE metal with LOW anisotropy
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 5: THE GADOLINIUM ANOMALY")
print("=" * 70)

print("""
Gd is a RARE EARTH metal (Z=64, 4f⁷ configuration) but has K1 = 0.012 MJ/m³,
lower than many 3d transition metals!

Why? Because Gd has L = 0 (half-filled 4f shell).
  4f⁷: S = 7/2, L = 0, J = 7/2

Spin-orbit coupling requires BOTH spin AND orbital angular momentum.
SOC = ξ × L · S

With L = 0, there is NO first-order SOC contribution to anisotropy.
The anisotropy comes entirely from HIGHER-ORDER effects (crystal field).

This is the STRONGEST evidence that K1 is SOC-dominated:
Even within the RE series, K1 tracks ORBITAL angular momentum, not
atomic number or lattice coherence.

RE anisotropy hierarchy:
  Gd (L=0): K1 = 0.012 MJ/m³  (no orbital contribution)
  Tb (L=3): K1 = 8.3 MJ/m³    (large orbital moment)
  Dy (L=5): K1 = 4.5 MJ/m³    (maximum orbital moment)
  Ho (L=6): K1 = 3.1 MJ/m³
  Er (L=6): K1 = 1.3 MJ/m³    (reduced by crystal field)

The ordering is Tb > Dy > Ho > Er >> Gd, which tracks |α_J × (J-1/2)|,
the Stevens factor for single-ion anisotropy, NOT θ_D or γ_phonon.
""")

# Stevens factors (α_J × J(J+1) related quantity)
# These are the single-ion anisotropy parameters
RE_data = {
    'Gd': (0, 0.012),       # L=0, alpha_J = 0
    'Tb': (3, 8.3),         # L=3, alpha_J = -1/99
    'Dy': (5, 4.5),         # L=5, alpha_J = -2/315
    'Ho': (6, 3.1),         # L=6, alpha_J = -1/450
    'Er': (6, 1.3),         # L=6, alpha_J = 4/1575
}

L_vals = np.array([RE_data[m][0] for m in RE_data])
K1_RE = np.array([RE_data[m][1] for m in RE_data])

# Correlation of K1 with L (excluding Gd, L=0)
L_nonzero = L_vals[L_vals > 0]
K1_nonzero = K1_RE[L_vals > 0]
if len(L_nonzero) >= 3:
    r_K_L, p_K_L = stats.pearsonr(L_nonzero, K1_nonzero)
    print(f"K1 vs L (RE metals, L>0): r = {r_K_L:.3f}")

# ==============================================================================
# Analysis 6: Magnetostriction — same SOC story
# ==============================================================================
print("\n" + "=" * 70)
print("ANALYSIS 6: MAGNETOSTRICTION λ_s (Same SOC Pattern)")
print("=" * 70)

# Only materials with nonzero λ_s
nonzero_ls = lambda_s > 0
if np.sum(nonzero_ls) >= 3:
    log_ls = np.log10(lambda_s[nonzero_ls] + 1e-10)

    r_ls_gamma, _ = stats.pearsonr(gamma_phonon[nonzero_ls], log_ls)
    r_ls_SOC, _ = stats.pearsonr(log_SOC[nonzero_ls], log_ls)
    r_ls_Z4, _ = stats.pearsonr(np.log10(Z[nonzero_ls]**4), log_ls)
    r_ls_D, _ = stats.pearsonr(np.log10(D[nonzero_ls]), log_ls)

    print(f"\nMagnetostriction |λ_s| correlations (N={np.sum(nonzero_ls)}):")
    print(f"  log|λ_s| vs γ_phonon:  r = {r_ls_gamma:.3f}  (coherence)")
    print(f"  log|λ_s| vs log(SOC):  r = {r_ls_SOC:.3f}  (atomic)")
    print(f"  log|λ_s| vs log(Z⁴):   r = {r_ls_Z4:.3f}  (atomic number)")
    print(f"  log|λ_s| vs log(D):    r = {r_ls_D:.3f}  (dominance parameter)")

    print(f"\n  {'SOC DOMINATES' if abs(r_ls_SOC) > abs(r_ls_gamma) else 'COHERENCE CONTRIBUTES'} magnetostriction")
    print(f"  SOC predictor: r = {r_ls_SOC:.3f}")
    print(f"  γ predictor:   r = {r_ls_gamma:.3f}")
    print(f"  Ratio: {abs(r_ls_SOC)/max(abs(r_ls_gamma), 0.001):.1f}×")

# ==============================================================================
# The SOC ∝ Z^4 relationship
# ==============================================================================
print("\n" + "=" * 70)
print("THE Z⁴ SCALING LAW")
print("=" * 70)

# SOC energy ξ ∝ Z^4 for hydrogen-like atoms
# Real atoms: ξ ∝ Z_eff^4 where Z_eff accounts for screening
# Test: does our data follow this?

r_SOC_Z, _ = stats.pearsonr(np.log10(Z), np.log10(SOC))
slope_SOC, intercept_SOC = np.polyfit(np.log10(Z), np.log10(SOC), 1)
print(f"\nlog(SOC) vs log(Z): r = {r_SOC_Z:.3f}")
print(f"Power law: SOC ∝ Z^{slope_SOC:.2f}")
print(f"Expected: SOC ∝ Z^4 for hydrogen-like; actual Z^{slope_SOC:.2f}")
print(f"(Screening reduces effective power law from 4 to {slope_SOC:.2f})")

# ==============================================================================
# Summary: When does coherence work for magnetic properties?
# ==============================================================================
print("\n" + "=" * 70)
print("SUMMARY: WHEN DOES γ WORK FOR MAGNETIC PROPERTIES?")
print("=" * 70)

print(f"""
DOMINANCE PARAMETER: D = ξ_SOC / (k_B × θ_D)

  D < {D_threshold:.0f}:  Coherence regime — γ_phonon may contribute
    Examples: {', '.join([names[i] for i in range(len(names)) if not soc_dominated[i]])}
    Within this subset: r(K1 vs γ) = {r_coh:.3f}

  D > {D_threshold:.0f}:  SOC regime — γ_phonon is irrelevant
    Examples: {', '.join([names[i] for i in range(len(names)) if soc_dominated[i]])}
    Within this subset: r(K1 vs γ) = {r_soc:.3f} (spurious)

OVERALL COMPARISON:
  K1 vs γ_phonon (coherence): r = {r_K_gamma:.3f}
  K1 vs SOC (atomic):         r = {r_K_SOC:.3f}
  K1 vs Z⁴ (atomic number):   r = {r_K_Z4:.3f}

  SOC outperforms γ by {abs(r_K_SOC)/max(abs(r_K_gamma), 0.001):.1f}× for predicting K1.

PHYSICAL EXPLANATION:
  Magnetocrystalline anisotropy arises from spin-orbit coupling:
    K ∝ ξ_SOC² / W  (second-order perturbation theory)

  where ξ_SOC ∝ Z⁴ is the SOC constant and W is the bandwidth.
  The coherence parameter γ has NO direct relationship to ξ_SOC.

  For 3d metals (Z ~ 26): ξ ~ 50 meV, k_Bθ_D ~ 40 meV → D ~ 1
    SOC and lattice energies are comparable → some coherence contribution

  For RE metals (Z ~ 65): ξ ~ 300 meV, k_Bθ_D ~ 15 meV → D ~ 20
    SOC completely dominates → coherence is irrelevant

  For 5d compounds (Z ~ 78): ξ ~ 600 meV, k_Bθ_D ~ 25 meV → D ~ 25
    SOC completely dominates → coherence is irrelevant

THE GADOLINIUM TEST:
  Gd (Z=64, L=0) has K1 = 0.012 — comparable to 3d metals!
  This proves K1 is NOT determined by atomic number alone,
  but by L×S (orbital angular momentum × spin).
  Even Z⁴ fails when L=0. SOC requires orbital angular momentum.

FRAMEWORK RECOMMENDATION:
  Before applying γ to any magnetic property:
  1. Calculate D = ξ_SOC / (k_B × θ_D) for the relevant element
  2. If D > 5: Do NOT use γ_phonon — use SOC theory instead
  3. If D < 5: γ_phonon MAY contribute, but check for SOC corrections
  4. Special case: If L=0 (half-filled subshell), SOC effect vanishes
""")

# ==============================================================================
# Visualization
# ==============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Phase 2 Session #4: Spin-Orbit Coupling Dominance\nWhen Atomic Physics Overwhelms Collective Coherence',
             fontsize=14, fontweight='bold')

class_colors = {'3d': 'blue', 'RE': 'red', 'RE-TM': 'orange', '5d': 'purple', 'ferrite': 'green'}

# Plot 1: K1 vs γ_phonon (FAILS)
ax = axes[0, 0]
for cls in class_colors:
    mask = np.array([c == cls for c in classes])
    if any(mask):
        ax.scatter(gamma_phonon[mask], K1[mask], c=class_colors[cls], s=80, alpha=0.7, label=cls)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_phonon[i], K1[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('|K₁| (MJ/m³)')
ax.set_yscale('log')
ax.set_title(f'K₁ vs γ_phonon: r = {r_K_gamma:.3f}\n(coherence FAILS)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2: K1 vs SOC (WORKS)
ax = axes[0, 1]
for cls in class_colors:
    mask = np.array([c == cls for c in classes])
    if any(mask):
        ax.scatter(SOC[mask], K1[mask], c=class_colors[cls], s=80, alpha=0.7, label=cls)
for i, name in enumerate(names):
    ax.annotate(name, (SOC[i], K1[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('SOC energy ξ (meV)')
ax.set_ylabel('|K₁| (MJ/m³)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f'K₁ vs SOC: r = {r_K_SOC:.3f}\n(atomic WORKS)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 3: Dominance parameter D distribution
ax = axes[0, 2]
bar_colors = [class_colors[c] for c in classes]
sorted_idx = np.argsort(D)
ax.barh(range(len(names)), D[sorted_idx], color=[bar_colors[i] for i in sorted_idx], alpha=0.7)
ax.set_yticks(range(len(names)))
ax.set_yticklabels([names[i] for i in sorted_idx], fontsize=8)
ax.axvline(x=D_threshold, color='black', linestyle='--', linewidth=2, label=f'D = {D_threshold}')
ax.set_xlabel('D = ξ_SOC / k_B θ_D')
ax.set_title('Dominance Parameter D\n(left: coherence-relevant, right: SOC-dominated)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='x')

# Plot 4: K1 vs γ — colored by D
ax = axes[1, 0]
scatter = ax.scatter(gamma_phonon, K1, c=np.log10(D), cmap='RdYlBu_r', s=100, alpha=0.8, edgecolors='k', linewidth=0.5)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_phonon[i], K1[i]), fontsize=6, alpha=0.6)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('|K₁| (MJ/m³)')
ax.set_yscale('log')
ax.set_title('K₁ vs γ, colored by log(D)\n(blue=coherence, red=SOC)')
plt.colorbar(scatter, ax=ax, label='log₁₀(D)')
ax.grid(True, alpha=0.3)

# Plot 5: SOC vs Z — the Z^n scaling
ax = axes[1, 1]
for cls in class_colors:
    mask = np.array([c == cls for c in classes])
    if any(mask):
        ax.scatter(Z[mask], SOC[mask], c=class_colors[cls], s=80, alpha=0.7, label=cls)
Z_fit = np.linspace(Z.min(), Z.max(), 100)
SOC_fit = 10**(slope_SOC * np.log10(Z_fit) + intercept_SOC)
ax.plot(Z_fit, SOC_fit, 'k--', alpha=0.5, label=f'SOC ∝ Z^{slope_SOC:.1f}')
ax.set_xlabel('Atomic Number Z')
ax.set_ylabel('SOC energy ξ (meV)')
ax.set_title(f'SOC vs Z: ξ ∝ Z^{slope_SOC:.1f}\n(r = {r_SOC_Z:.3f})')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 6: Summary
ax = axes[1, 2]
ax.text(0.5, 0.92, 'SOC Dominance Summary', fontsize=14, ha='center', fontweight='bold', transform=ax.transAxes)

ax.text(0.5, 0.78, 'Predicting K₁ (magnetic anisotropy):', fontsize=10, ha='center', transform=ax.transAxes)
ax.text(0.5, 0.68, f'γ_phonon: r = {r_K_gamma:.3f}', fontsize=11, ha='center', transform=ax.transAxes, color='blue')
ax.text(0.5, 0.58, f'SOC (ξ):  r = {r_K_SOC:.3f}', fontsize=11, ha='center', transform=ax.transAxes, color='red')
ax.text(0.5, 0.48, f'Z⁴:       r = {r_K_Z4:.3f}', fontsize=11, ha='center', transform=ax.transAxes, color='darkred')

ax.text(0.5, 0.33, 'Dominance threshold: D = 5', fontsize=10, ha='center', transform=ax.transAxes)
ax.text(0.5, 0.22, f'D < 5: γ may work (3d, ferrites)', fontsize=10, ha='center', color='blue', transform=ax.transAxes)
ax.text(0.5, 0.12, f'D > 5: SOC dominates (RE, 5d)', fontsize=10, ha='center', color='red', transform=ax.transAxes)

ax.text(0.5, 0.02, 'Gd anomaly: Z=64 but K₁≈3d metals (L=0!)', fontsize=9, ha='center',
        color='darkgreen', fontweight='bold', transform=ax.transAxes)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phase2_soc_dominance.png',
            dpi=150, bbox_inches='tight')
plt.close()
print("\nFigure saved: phase2_soc_dominance.png")
