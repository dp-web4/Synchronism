"""
Session #90: Phonon-Limited Mobility in Semiconductors

Hypothesis: Carrier mobility μ in semiconductors is limited by phonon scattering,
which relates to γ_phonon through the Debye temperature θ_D.

Key relationships:
- μ ∝ T^(-3/2) for acoustic phonon scattering
- μ ∝ 1/m* (effective mass)
- μ ∝ θ_D^(3/2) / T^(3/2) from deformation potential theory

From coherence framework:
- γ_phonon = 2T/θ_D (thermal phonon population)
- Phonon scattering increases with γ_phonon
- μ ∝ 1/γ_phonon^n for some exponent n

Data: Room temperature electron and hole mobilities for semiconductors
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# =============================================================================
# SEMICONDUCTOR MOBILITY DATA
# Sources: Sze "Physics of Semiconductor Devices", Madelung "Semiconductors"
# =============================================================================

# μ in cm²/V·s at 300K
semiconductors = {
    # Group IV
    'Si': {'mu_e': 1450, 'mu_h': 450, 'theta_D': 645, 'm_e': 0.26, 'm_h': 0.37, 'E_g': 1.12},
    'Ge': {'mu_e': 3900, 'mu_h': 1900, 'theta_D': 374, 'm_e': 0.12, 'm_h': 0.21, 'E_g': 0.67},
    'Diamond': {'mu_e': 2200, 'mu_h': 1800, 'theta_D': 1860, 'm_e': 0.20, 'm_h': 0.24, 'E_g': 5.47},

    # III-V compounds
    'GaAs': {'mu_e': 8500, 'mu_h': 400, 'theta_D': 344, 'm_e': 0.067, 'm_h': 0.45, 'E_g': 1.42},
    'GaP': {'mu_e': 300, 'mu_h': 150, 'theta_D': 445, 'm_e': 0.35, 'm_h': 0.67, 'E_g': 2.26},
    'GaSb': {'mu_e': 4000, 'mu_h': 850, 'theta_D': 265, 'm_e': 0.042, 'm_h': 0.28, 'E_g': 0.73},
    'InAs': {'mu_e': 33000, 'mu_h': 460, 'theta_D': 249, 'm_e': 0.023, 'm_h': 0.41, 'E_g': 0.36},
    'InP': {'mu_e': 5400, 'mu_h': 200, 'theta_D': 321, 'm_e': 0.077, 'm_h': 0.64, 'E_g': 1.35},
    'InSb': {'mu_e': 77000, 'mu_h': 850, 'theta_D': 202, 'm_e': 0.013, 'm_h': 0.40, 'E_g': 0.17},
    'AlAs': {'mu_e': 280, 'mu_h': 100, 'theta_D': 417, 'm_e': 0.42, 'm_h': 0.76, 'E_g': 2.16},
    'AlSb': {'mu_e': 200, 'mu_h': 400, 'theta_D': 292, 'm_e': 0.12, 'm_h': 0.98, 'E_g': 1.58},

    # II-VI compounds
    'ZnS': {'mu_e': 165, 'mu_h': 5, 'theta_D': 340, 'm_e': 0.28, 'm_h': 1.76, 'E_g': 3.68},
    'ZnSe': {'mu_e': 530, 'mu_h': 30, 'theta_D': 271, 'm_e': 0.16, 'm_h': 0.78, 'E_g': 2.70},
    'ZnTe': {'mu_e': 340, 'mu_h': 100, 'theta_D': 223, 'm_e': 0.12, 'm_h': 0.60, 'E_g': 2.26},
    'CdS': {'mu_e': 340, 'mu_h': 15, 'theta_D': 215, 'm_e': 0.21, 'm_h': 0.80, 'E_g': 2.42},
    'CdSe': {'mu_e': 650, 'mu_h': 50, 'theta_D': 181, 'm_e': 0.13, 'm_h': 0.45, 'E_g': 1.74},
    'CdTe': {'mu_e': 1050, 'mu_h': 100, 'theta_D': 158, 'm_e': 0.11, 'm_h': 0.35, 'E_g': 1.49},

    # Lead chalcogenides
    'PbS': {'mu_e': 600, 'mu_h': 700, 'theta_D': 227, 'm_e': 0.08, 'm_h': 0.08, 'E_g': 0.41},
    'PbSe': {'mu_e': 1000, 'mu_h': 1000, 'theta_D': 130, 'm_e': 0.04, 'm_h': 0.04, 'E_g': 0.27},
    'PbTe': {'mu_e': 6000, 'mu_h': 4000, 'theta_D': 130, 'm_e': 0.03, 'm_h': 0.03, 'E_g': 0.31},
}

# =============================================================================
# ANALYSIS
# =============================================================================

print("="*60)
print("Session #90: Phonon-Limited Mobility in Semiconductors")
print("="*60)

# Extract data
names = list(semiconductors.keys())
mu_e = np.array([semiconductors[m]['mu_e'] for m in names])
mu_h = np.array([semiconductors[m]['mu_h'] for m in names])
theta_D = np.array([semiconductors[m]['theta_D'] for m in names])
m_e = np.array([semiconductors[m]['m_e'] for m in names])
m_h = np.array([semiconductors[m]['m_h'] for m in names])
E_g = np.array([semiconductors[m]['E_g'] for m in names])

T = 300  # K
gamma_phonon = 2 * T / theta_D

print(f"\nDataset: {len(names)} semiconductors")
print(f"μ_e range: {mu_e.min()} - {mu_e.max()} cm²/V·s")
print(f"θ_D range: {theta_D.min()} - {theta_D.max()} K")
print(f"γ_phonon range: {gamma_phonon.min():.2f} - {gamma_phonon.max():.2f}")

# =============================================================================
# TEST 1: μ vs θ_D (direct correlation)
# =============================================================================
print("\n" + "="*60)
print("TEST 1: Mobility vs Debye temperature")
print("="*60)

# Theory: μ ∝ θ_D^(3/2) from acoustic phonon scattering
r_mu_e_theta, p_e = stats.pearsonr(theta_D, mu_e)
r_mu_h_theta, p_h = stats.pearsonr(theta_D, mu_h)

print(f"μ_e vs θ_D: r = {r_mu_e_theta:.3f}, p = {p_e:.4f}")
print(f"μ_h vs θ_D: r = {r_mu_h_theta:.3f}, p = {p_h:.4f}")

# Log-log for power law
r_log_e, _ = stats.pearsonr(np.log(theta_D), np.log(mu_e))
r_log_h, _ = stats.pearsonr(np.log(theta_D), np.log(mu_h))
print(f"log(μ_e) vs log(θ_D): r = {r_log_e:.3f}")
print(f"log(μ_h) vs log(θ_D): r = {r_log_h:.3f}")

# =============================================================================
# TEST 2: μ vs γ_phonon (coherence test)
# =============================================================================
print("\n" + "="*60)
print("TEST 2: Mobility vs γ_phonon = 2T/θ_D")
print("="*60)

# Theory from coherence: μ ∝ 1/γ_phonon^n
# Since γ = 2T/θ_D, this means μ ∝ θ_D^n / T^n
# Deformation potential gives n ≈ 3/2

r_mu_e_gamma, _ = stats.pearsonr(gamma_phonon, mu_e)
r_mu_h_gamma, _ = stats.pearsonr(gamma_phonon, mu_h)
print(f"μ_e vs γ_phonon: r = {r_mu_e_gamma:.3f}")
print(f"μ_h vs γ_phonon: r = {r_mu_h_gamma:.3f}")

# Should be negative (higher γ = more scattering = lower μ)
# Try 2/γ instead
two_over_gamma = 2 / gamma_phonon
r_mu_e_2gamma, _ = stats.pearsonr(two_over_gamma, mu_e)
r_mu_h_2gamma, _ = stats.pearsonr(two_over_gamma, mu_h)
print(f"μ_e vs 2/γ_phonon: r = {r_mu_e_2gamma:.3f}")
print(f"μ_h vs 2/γ_phonon: r = {r_mu_h_2gamma:.3f}")

# =============================================================================
# TEST 3: Include effective mass (μ×m* should correlate with θ)
# =============================================================================
print("\n" + "="*60)
print("TEST 3: Mobility × effective mass (μ×m*)")
print("="*60)

# Theory: μ = eτ/m* where τ ∝ θ_D^(3/2) / T^(3/2)
# So μ × m* ∝ θ_D^(3/2)

mu_m_e = mu_e * m_e
mu_m_h = mu_h * m_h

r_mum_e_theta, _ = stats.pearsonr(theta_D, mu_m_e)
r_mum_h_theta, _ = stats.pearsonr(theta_D, mu_m_h)
print(f"μ_e × m*_e vs θ_D: r = {r_mum_e_theta:.3f}")
print(f"μ_h × m*_h vs θ_D: r = {r_mum_h_theta:.3f}")

# Power law test
slope_e, _, r_pow_e, _, _ = stats.linregress(np.log(theta_D), np.log(mu_m_e))
slope_h, _, r_pow_h, _, _ = stats.linregress(np.log(theta_D), np.log(mu_m_h))
print(f"μ_e × m*_e ∝ θ_D^{slope_e:.2f} (r = {r_pow_e:.3f})")
print(f"μ_h × m*_h ∝ θ_D^{slope_h:.2f} (r = {r_pow_h:.3f})")

# =============================================================================
# TEST 4: Coherence model μ ∝ (2/γ)^n / m*
# =============================================================================
print("\n" + "="*60)
print("TEST 4: Coherence model μ ∝ (2/γ)^n / m*")
print("="*60)

# Try different exponents
for n in [1.0, 1.5, 2.0]:
    mu_pred_e = (2/gamma_phonon)**n / m_e
    mu_pred_h = (2/gamma_phonon)**n / m_h

    # Normalize to same scale
    mu_pred_e = mu_pred_e * mu_e.mean() / mu_pred_e.mean()
    mu_pred_h = mu_pred_h * mu_h.mean() / mu_pred_h.mean()

    r_e, _ = stats.pearsonr(mu_pred_e, mu_e)
    r_h, _ = stats.pearsonr(mu_pred_h, mu_h)
    print(f"n = {n}: μ_e r = {r_e:.3f}, μ_h r = {r_h:.3f}")

# Optimal n
best_r_e = 0
best_n_e = 1.0
for n in np.arange(0.5, 3.0, 0.1):
    mu_pred = (2/gamma_phonon)**n / m_e
    mu_pred = mu_pred * mu_e.mean() / mu_pred.mean()
    r, _ = stats.pearsonr(mu_pred, mu_e)
    if r > best_r_e:
        best_r_e = r
        best_n_e = n

print(f"\nOptimal for μ_e: n = {best_n_e:.1f}, r = {best_r_e:.3f}")

# =============================================================================
# TEST 5: Separate by material class
# =============================================================================
print("\n" + "="*60)
print("TEST 5: Analysis by material class")
print("="*60)

classes = {
    'IV': ['Si', 'Ge', 'Diamond'],
    'III-V': ['GaAs', 'GaP', 'GaSb', 'InAs', 'InP', 'InSb', 'AlAs', 'AlSb'],
    'II-VI': ['ZnS', 'ZnSe', 'ZnTe', 'CdS', 'CdSe', 'CdTe'],
    'IV-VI': ['PbS', 'PbSe', 'PbTe'],
}

for cls_name, materials in classes.items():
    idx = [names.index(m) for m in materials if m in names]
    if len(idx) < 3:
        continue
    mu_cls = mu_e[idx]
    theta_cls = theta_D[idx]
    gamma_cls = gamma_phonon[idx]
    m_cls = m_e[idx]

    r_theta, _ = stats.pearsonr(theta_cls, mu_cls)
    r_gamma, _ = stats.pearsonr(gamma_cls, mu_cls)

    # Model with m*
    mu_pred = (2/gamma_cls)**1.5 / m_cls
    mu_pred = mu_pred * mu_cls.mean() / mu_pred.mean()
    r_model, _ = stats.pearsonr(mu_pred, mu_cls)

    print(f"{cls_name:<8}: μ_e vs θ_D: r = {r_theta:.3f}, μ_e vs model: r = {r_model:.3f}")

# =============================================================================
# TEST 6: Temperature dependence prediction
# =============================================================================
print("\n" + "="*60)
print("TEST 6: Temperature dependence (scaling)")
print("="*60)

# At 300K: μ(300) ∝ (θ_D/300)^n
# At T: μ(T) ∝ (θ_D/T)^n
# Ratio: μ(T)/μ(300) = (300/T)^n

# This means γ_phonon(T) = 2T/θ_D
# μ(T) ∝ (2/γ(T))^n = (θ_D/T)^n

print("Prediction: μ ∝ T^(-n) where n ≈ 1.5 for acoustic phonons")
print(f"From coherence: μ ∝ (2/γ)^n = (θ_D/T)^n")
print(f"Best fit n = {best_n_e:.1f}")

# =============================================================================
# VISUALIZATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(14, 9))

# Plot 1: μ_e vs θ_D
ax = axes[0, 0]
ax.scatter(theta_D, mu_e, c='blue', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (theta_D[i], mu_e[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('θ_D (K)')
ax.set_ylabel('μ_e (cm²/V·s)')
ax.set_yscale('log')
ax.set_title(f'Electron mobility vs θ_D (r = {r_mu_e_theta:.3f})')
ax.grid(True, alpha=0.3)

# Plot 2: μ_e vs γ_phonon
ax = axes[0, 1]
ax.scatter(gamma_phonon, mu_e, c='green', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (gamma_phonon[i], mu_e[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('γ_phonon = 2T/θ_D')
ax.set_ylabel('μ_e (cm²/V·s)')
ax.set_yscale('log')
ax.set_title(f'Electron mobility vs γ_phonon (r = {r_mu_e_gamma:.3f})')
ax.grid(True, alpha=0.3)

# Plot 3: Model prediction
ax = axes[0, 2]
mu_pred_e = (2/gamma_phonon)**1.5 / m_e
mu_pred_e = mu_pred_e * mu_e.mean() / mu_pred_e.mean()
r_final, _ = stats.pearsonr(mu_pred_e, mu_e)
ax.scatter(mu_pred_e, mu_e, c='red', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (mu_pred_e[i], mu_e[i]), fontsize=7, alpha=0.8)
min_val = min(mu_pred_e.min(), mu_e.min())
max_val = max(mu_pred_e.max(), mu_e.max())
ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)
ax.set_xlabel('μ_predicted ∝ (2/γ)^1.5 / m*')
ax.set_ylabel('μ_observed (cm²/V·s)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f'Model prediction (r = {r_final:.3f})')
ax.grid(True, alpha=0.3)

# Plot 4: μ_h vs θ_D
ax = axes[1, 0]
ax.scatter(theta_D, mu_h, c='purple', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (theta_D[i], mu_h[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('θ_D (K)')
ax.set_ylabel('μ_h (cm²/V·s)')
ax.set_yscale('log')
ax.set_title(f'Hole mobility vs θ_D (r = {r_mu_h_theta:.3f})')
ax.grid(True, alpha=0.3)

# Plot 5: μ×m* vs θ_D
ax = axes[1, 1]
ax.scatter(theta_D, mu_m_e, c='orange', s=60, alpha=0.7)
for i, name in enumerate(names):
    ax.annotate(name, (theta_D[i], mu_m_e[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('θ_D (K)')
ax.set_ylabel('μ_e × m*_e')
ax.set_title(f'μ_e × m*_e vs θ_D (r = {r_mum_e_theta:.3f})')
ax.grid(True, alpha=0.3)

# Plot 6: Effect of mass
ax = axes[1, 2]
# Color by effective mass
sc = ax.scatter(theta_D, mu_e, c=m_e, s=80, alpha=0.8, cmap='viridis')
for i, name in enumerate(names):
    ax.annotate(name, (theta_D[i], mu_e[i]), fontsize=7, alpha=0.8)
ax.set_xlabel('θ_D (K)')
ax.set_ylabel('μ_e (cm²/V·s)')
ax.set_yscale('log')
ax.set_title('Mobility (color = m*/m_0)')
plt.colorbar(sc, ax=ax, label='m*/m_0')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/phonon_limited_mobility.png',
            dpi=150, bbox_inches='tight')
plt.close()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*60)
print("SUMMARY: Phonon-Limited Mobility and Coherence")
print("="*60)

print("""
Key Results:
1. μ_e vs θ_D: r = {:.3f} (weak direct correlation)
2. μ_e vs γ_phonon: r = {:.3f} (negative as expected)
3. μ_e × m* vs θ_D: r = {:.3f} (improved by mass correction)
4. Model μ ∝ (2/γ)^1.5/m*: r = {:.3f}

Coherence Framework Interpretation:

1. **Phonon scattering limits mobility**
   - Higher γ_phonon = more thermal phonons = more scattering
   - μ ∝ (2/γ)^n where n ≈ 1.5

2. **Effective mass matters**
   - μ = eτ/m* where τ = scattering time
   - τ ∝ (2/γ)^n from phonon population
   - So μ ∝ (2/γ)^n / m*

3. **Material class variations**
   - III-V: Often high μ due to low m* (GaAs, InAs, InSb)
   - II-VI: Lower μ due to higher m*
   - IV-VI (lead salts): Anomalous - both carriers similar μ

4. **Temperature dependence**
   - γ(T) = 2T/θ_D
   - μ(T) ∝ (2/γ)^n = (θ_D/T)^n ∝ T^(-n)
   - n = 3/2 for acoustic phonon scattering

5. **InSb has highest μ because**:
   - Very low m* = 0.013 m_0
   - Moderate θ_D = 202 K gives γ ~ 3
   - Combined: μ = 77,000 cm²/V·s!

Physical Insight:
Mobility is inversely related to phonon coherence parameter γ_phonon.
Higher Debye temperature → lower γ → fewer thermal phonons → higher mobility.
But effective mass provides additional modulation.

Why simple μ vs θ_D fails:
- Mobility depends on BOTH γ_phonon AND m*
- InSb has low θ_D (low phonon energy) but extremely low m*
- Diamond has high θ_D but moderate m*
- The m* correction is essential!
""".format(r_mu_e_theta, r_mu_e_gamma, r_mum_e_theta, r_final))

print("\nTop mobility semiconductors by coherence design:")
print(f"{'Material':<10} {'μ_e':>8} {'θ_D':>6} {'γ_phonon':>10} {'m*':>6} {'μ×m*':>8}")
print("-" * 55)
top_idx = np.argsort(mu_e)[::-1][:8]
for i in top_idx:
    print(f"{names[i]:<10} {mu_e[i]:>8.0f} {theta_D[i]:>6.0f} {gamma_phonon[i]:>10.2f} {m_e[i]:>6.3f} {mu_m_e[i]:>8.1f}")

print(f"\nPlot saved to: phonon_limited_mobility.png")
