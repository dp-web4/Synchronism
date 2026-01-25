"""
Chemistry Session #202: Dielectric Relaxation at γ ~ 1
Analyzing frequency-dependent polarization through the coherence framework.

Key hypothesis: ωτ = 1 IS the γ ~ 1 crossover
- Below ωτ = 1: dipoles follow field (coherent)
- Above ωτ = 1: dipoles can't follow (incoherent)
- At ωτ = 1: maximum energy dissipation (ε'' peak)

γ parameters:
- ωτ (frequency × relaxation time)
- Cole-Cole α parameter
- Arrhenius Ea/RT for τ
"""

import numpy as np
import matplotlib.pyplot as plt

print("=" * 60)
print("CHEMISTRY SESSION #202: DIELECTRIC RELAXATION AT γ ~ 1")
print("=" * 60)

# ============================================================
# 1. DEBYE RELAXATION MODEL
# ============================================================
print("\n" + "=" * 60)
print("1. DEBYE RELAXATION MODEL")
print("=" * 60)

print("""
Debye equation:
  ε*(ω) = ε_∞ + (ε_s - ε_∞) / (1 + iωτ)

where:
  ε_s = static (low-frequency) permittivity
  ε_∞ = high-frequency permittivity
  τ = relaxation time

Real and imaginary parts:
  ε'(ω) = ε_∞ + (ε_s - ε_∞) / (1 + ω²τ²)
  ε''(ω) = (ε_s - ε_∞) × ωτ / (1 + ω²τ²)

γ_Debye = ωτ (dimensionless frequency)

At γ = 1 (ωτ = 1):
  - ε'' is MAXIMUM (peak loss)
  - ε' is halfway between ε_s and ε_∞
  - Maximum energy dissipation
  - This IS the coherence transition!
""")

# Calculate Debye response
omega_tau = np.logspace(-3, 3, 200)  # ωτ range
eps_s, eps_inf = 80, 5  # Water-like values
delta_eps = eps_s - eps_inf

eps_prime = eps_inf + delta_eps / (1 + omega_tau**2)
eps_double = delta_eps * omega_tau / (1 + omega_tau**2)

# Find maximum ε''
idx_max = np.argmax(eps_double)
print(f"\nMaximum ε'' at ωτ = {omega_tau[idx_max]:.3f}")
print(f"ε''_max = {eps_double[idx_max]:.1f}")
print(f"ε' at maximum = {eps_prime[idx_max]:.1f} (= (ε_s + ε_∞)/2 = {(eps_s+eps_inf)/2:.1f})")

# ============================================================
# 2. RELAXATION TIMES OF LIQUIDS
# ============================================================
print("\n" + "=" * 60)
print("2. RELAXATION TIMES OF LIQUIDS")
print("=" * 60)

print("""
Relaxation time τ depends on:
  - Viscosity η
  - Molecular size (Debye-Stokes-Einstein)
  - Temperature

τ = 4πηr³ / (kT)  (Debye equation)

For water at 25°C: τ ≈ 8 ps
This gives peak frequency f = 1/(2πτ) ≈ 20 GHz
""")

# Relaxation time data (ps)
tau_data = {
    # (τ at 25°C in ps, activation energy kJ/mol)
    'Water': (8.3, 17),
    'Methanol': (51, 15),
    'Ethanol': (170, 20),
    'Propanol': (400, 24),
    'Glycerol': (1e6, 75),
    'Acetone': (3.2, 8),
    'Acetonitrile': (3.4, 9),
    'DMSO': (21, 18),
    'Formamide': (37, 23),
    'NMF': (120, 28),
}

print("\nDielectric Relaxation Times:")
print("-" * 60)
print(f"{'Liquid':<15} {'τ (ps)':<12} {'f_peak (GHz)':<15} {'E_a (kJ/mol)'}")
print("-" * 60)

for liquid, (tau, E_a) in tau_data.items():
    f_peak = 1e3 / (2 * np.pi * tau)  # GHz
    print(f"{liquid:<15} {tau:<12.1f} {f_peak:<15.2f} {E_a}")

print("-" * 60)

# ============================================================
# 3. TEMPERATURE DEPENDENCE
# ============================================================
print("\n" + "=" * 60)
print("3. TEMPERATURE DEPENDENCE (ARRHENIUS)")
print("=" * 60)

print("""
Arrhenius behavior:
  τ = τ₀ × exp(E_a/RT)

γ_Arr = E_a / RT

At γ = 1: thermal energy equals activation barrier
For water: E_a ≈ 17 kJ/mol
At T = 298K: RT = 2.5 kJ/mol
So γ ≈ 6.8 (moderately activated)

The relaxation process is thermally activated.
""")

R = 8.314  # J/(mol·K)
T = 298  # K
RT = R * T / 1000  # kJ/mol

print("\nArrhenius γ = E_a/RT:")
print("-" * 45)

gamma_Arr = []
for liquid, (tau, E_a) in tau_data.items():
    gamma = E_a / RT
    gamma_Arr.append(gamma)
    print(f"{liquid:<15} γ = {gamma:.1f}")

print("-" * 45)
print(f"Mean γ_Arr = {np.mean(gamma_Arr):.1f} ± {np.std(gamma_Arr):.1f}")

# ============================================================
# 4. COLE-COLE DISTRIBUTION
# ============================================================
print("\n" + "=" * 60)
print("4. COLE-COLE BROADENING")
print("=" * 60)

print("""
Cole-Cole equation (non-Debye):
  ε* = ε_∞ + (ε_s - ε_∞) / [1 + (iωτ)^(1-α)]

where α is the distribution parameter:
  α = 0: ideal Debye (single τ)
  α > 0: distribution of relaxation times

γ_CC = 1 - α (Cole-Cole parameter)

At γ = 1 (α = 0): ideal Debye relaxation
α > 0: stretched exponential, broader peak
""")

# Cole-Cole α values
cole_cole_data = {
    'Water': 0.02,
    'Methanol': 0.05,
    'Ethanol': 0.08,
    'Propanol': 0.12,
    'Glycerol': 0.25,
    'Acetone': 0.01,
    'DMSO': 0.03,
}

print("\nCole-Cole α Parameters:")
print("-" * 45)
print(f"{'Liquid':<15} {'α':<10} {'γ = 1-α':<12} {'Type'}")
print("-" * 45)

gamma_CC = []
for liquid, alpha in cole_cole_data.items():
    gamma = 1 - alpha
    gamma_CC.append(gamma)
    if alpha < 0.05:
        liq_type = "Debye-like (γ~1)"
    elif alpha < 0.15:
        liq_type = "Slightly broadened"
    else:
        liq_type = "Distributed"
    print(f"{liquid:<15} {alpha:<10.2f} {gamma:<12.2f} {liq_type}")

print("-" * 45)
print(f"Mean γ = 1-α = {np.mean(gamma_CC):.2f} ± {np.std(gamma_CC):.2f}")

# ============================================================
# 5. HAVRILIAK-NEGAMI
# ============================================================
print("\n" + "=" * 60)
print("5. HAVRILIAK-NEGAMI EQUATION")
print("=" * 60)

print("""
Havriliak-Negami (general form):
  ε* = ε_∞ + (ε_s - ε_∞) / [1 + (iωτ)^α]^β

Parameters:
  α: broadening (0-1, 1 = no broadening)
  β: asymmetry (0-1, 1 = symmetric)

Special cases:
  α=1, β=1: Debye
  α<1, β=1: Cole-Cole
  α=1, β<1: Davidson-Cole
  α<1, β<1: Havriliak-Negami

γ_HN = α × β (combined shape factor)
At γ = 1: ideal Debye behavior
""")

# HN parameters for polymers (more complex)
HN_data = {
    # (α, β)
    'PVAc': (0.82, 0.52),
    'PMMA': (0.70, 0.45),
    'Polystyrene': (0.75, 0.50),
    'PVC': (0.78, 0.48),
    'Glycerol': (0.90, 0.85),
    'Propylene carbonate': (0.95, 0.92),
}

print("\nHavriliak-Negami Parameters:")
print("-" * 55)
print(f"{'Material':<20} {'α':<8} {'β':<8} {'γ = αβ'}")
print("-" * 55)

gamma_HN = []
for material, (alpha, beta) in HN_data.items():
    gamma = alpha * beta
    gamma_HN.append(gamma)
    print(f"{material:<20} {alpha:<8.2f} {beta:<8.2f} {gamma:.2f}")

print("-" * 55)
print(f"Mean γ_HN = αβ = {np.mean(gamma_HN):.2f} ± {np.std(gamma_HN):.2f}")

# ============================================================
# 6. DIELECTRIC STRENGTH
# ============================================================
print("\n" + "=" * 60)
print("6. DIELECTRIC STRENGTH")
print("=" * 60)

print("""
Dielectric strength Δε = ε_s - ε_∞

Kirkwood correlation factor g:
  g = ε_s × (2ε_s + ε_∞) × M / [ρ × N_A × μ² × (ε_∞ + 2)²]

For uncorrelated dipoles: g = 1 (γ ~ 1!)
g > 1: parallel correlation (H-bonding)
g < 1: antiparallel correlation

Water: g ≈ 2.7 (strong H-bond correlation)
""")

# Kirkwood g factors
g_factor_data = {
    'Water': 2.7,
    'Methanol': 2.9,
    'Ethanol': 2.8,
    'Acetone': 1.0,
    'Acetonitrile': 0.95,
    'Chloroform': 0.85,
    'Benzene': 0.0,
    'DMSO': 1.8,
}

print("\nKirkwood g Factors:")
print("-" * 50)
print(f"{'Liquid':<15} {'g':<10} {'Type'}")
print("-" * 50)

for liquid, g in g_factor_data.items():
    if 0.9 <= g <= 1.1:
        liq_type = "Uncorrelated (γ~1)"
    elif g > 1.1:
        liq_type = "Parallel (H-bond)"
    else:
        liq_type = "Antiparallel"
    print(f"{liquid:<15} {g:<10.2f} {liq_type}")

# ============================================================
# 7. VISUALIZATION
# ============================================================
print("\n" + "=" * 60)
print("7. GENERATING VISUALIZATION")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Debye relaxation curves
ax1 = axes[0, 0]
ax1.semilogx(omega_tau, eps_prime, 'b-', linewidth=2, label="ε' (real)")
ax1.semilogx(omega_tau, eps_double, 'r-', linewidth=2, label="ε'' (loss)")
ax1.axvline(x=1, color='green', linestyle='--', linewidth=2, label='ωτ = 1 (γ ~ 1)')
ax1.scatter([1], [eps_double[np.argmin(np.abs(omega_tau-1))]], s=100, c='green', zorder=5)
ax1.set_xlabel('ωτ (dimensionless)', fontsize=12)
ax1.set_ylabel('Permittivity', fontsize=12)
ax1.set_title('Debye Relaxation: ε\'\' Peak at ωτ = 1', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Cole-Cole plot
ax2 = axes[0, 1]
ax2.plot(eps_prime, eps_double, 'b-', linewidth=2, label='Debye (α=0)')
# Add Cole-Cole with α = 0.2
alpha_cc = 0.2
theta = np.linspace(0, np.pi, 100)
r = (eps_s - eps_inf) / 2
center_x = (eps_s + eps_inf) / 2
center_y = 0
# Cole-Cole is a depressed semicircle
x_cc = center_x - r * np.cos(theta * (1-alpha_cc)) / np.cos(np.pi * alpha_cc / 2)
y_cc = r * np.sin(theta * (1-alpha_cc)) / np.cos(np.pi * alpha_cc / 2)
ax2.plot(eps_prime, eps_double, 'b-', linewidth=2, label='Debye (α=0, γ=1)')
ax2.scatter([(eps_s+eps_inf)/2], [delta_eps/2], s=100, c='red', zorder=5, label='ωτ = 1')
ax2.set_xlabel("ε' (real)", fontsize=12)
ax2.set_ylabel("ε'' (imaginary)", fontsize=12)
ax2.set_title('Cole-Cole Plot: Semicircle at γ ~ 1', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 90)
ax2.set_ylim(0, 50)

# Plot 3: Arrhenius γ distribution
ax3 = axes[1, 0]
liquids = list(tau_data.keys())
ax3.barh(range(len(liquids)), gamma_Arr, color='steelblue', edgecolor='black', alpha=0.7)
ax3.axvline(x=1, color='red', linestyle='--', linewidth=2, label='γ = 1 (kT = E_a)')
ax3.set_yticks(range(len(liquids)))
ax3.set_yticklabels(liquids, fontsize=9)
ax3.set_xlabel('γ = E_a / RT', fontsize=12)
ax3.set_title('Relaxation Activation: γ > 1 (Barrier)', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3, axis='x')

# Plot 4: Cole-Cole γ = 1-α distribution
ax4 = axes[1, 1]
ax4.hist(gamma_CC, bins=8, edgecolor='black', alpha=0.7, color='coral')
ax4.axvline(x=1, color='red', linestyle='--', linewidth=2, label='γ = 1 (Debye)')
ax4.axvline(x=np.mean(gamma_CC), color='green', linestyle='-', linewidth=2,
            label=f'Mean = {np.mean(gamma_CC):.2f}')
ax4.set_xlabel('γ = 1 - α (Cole-Cole)', fontsize=12)
ax4.set_ylabel('Count', fontsize=12)
ax4.set_title('Cole-Cole: α = 0 is γ ~ 1', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/dielectric_relaxation_coherence.png', dpi=150)
print("Saved: dielectric_relaxation_coherence.png")
plt.close()

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("SESSION #202 SUMMARY: DIELECTRIC RELAXATION AT γ ~ 1")
print("=" * 60)

print(f"""
KEY FINDINGS:

1. DEBYE CROSSOVER AT ωτ = 1
   γ_Debye = ωτ (dimensionless frequency)
   At γ = 1: ε'' is MAXIMUM (peak loss)
   Below γ: dipoles follow field (coherent)
   Above γ: dipoles can't follow (incoherent)

2. ARRHENIUS ACTIVATION
   γ_Arr = E_a / RT
   Mean γ = {np.mean(gamma_Arr):.1f} ± {np.std(gamma_Arr):.1f}
   Relaxation is thermally activated

3. COLE-COLE BROADENING
   γ_CC = 1 - α
   Mean γ = {np.mean(gamma_CC):.2f} ± {np.std(gamma_CC):.2f}
   α = 0 (γ = 1): ideal Debye relaxation
   Water, acetone: α ~ 0.02 (nearly ideal!)

4. HAVRILIAK-NEGAMI
   γ_HN = α × β (shape factor)
   Mean γ = {np.mean(gamma_HN):.2f}
   Small molecules: γ ~ 0.8-0.9 (near Debye)
   Polymers: γ ~ 0.4-0.5 (broader distribution)

5. KIRKWOOD FACTOR
   g = 1: uncorrelated dipoles (γ ~ 1!)
   Acetone, acetonitrile: g ~ 1.0
   H-bonding (water): g ~ 2.7 (correlated)

CENTRAL INSIGHT:
The Debye relaxation crossover ωτ = 1 IS γ ~ 1:
- Maximum energy dissipation at the transition
- ε' changes from ε_s to ε_∞ through midpoint
- Dipole dynamics: coherent → incoherent

The Cole-Cole α = 0 (γ = 1) represents ideal
single-relaxation-time behavior - perfect coherence.

This is the 65th phenomenon type at γ ~ 1!
""")

print("=" * 60)
print("SESSION #202 COMPLETE")
print("=" * 60)
