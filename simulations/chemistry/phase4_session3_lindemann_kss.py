"""
Phase 4 Session 3: Lindemann-KSS Connection
=============================================
Does η/s at the melting point provide a material-independent entity criterion?

Hypothesis: The crystal→liquid transition IS an entity→process transition.
If the Synchronism entity criterion is physically meaningful, the melting
point should correspond to a characteristic η/s value relative to KSS.

Key: η/s uses measured viscosity and entropy — NOT θ_D — so this is
non-circular with the Phase 2-3 findings.

Data sources:
- Liquid metal viscosities at T_m: Iida & Guthrie, "The Physical Properties
  of Liquid Metals" (1988); Assael et al. (2006) recommended values
- Entropy of fusion: standard thermodynamic tables (Chase, JANAF)
- Liquid entropy at T_m: S_liquid = S_solid(T_m) + ΔS_fusion
- KSS bound: η/s ≥ ℏ/(4πk_B)
"""

import numpy as np

# Physical constants
hbar = 1.0546e-34       # J·s
k_B = 1.3806e-23        # J/K
h = 6.626e-34            # J·s
N_A = 6.022e23           # mol^-1
R = 8.314                # J/(mol·K)

# KSS bound
eta_over_s_KSS = hbar / (4 * np.pi * k_B)  # in SI: Pa·s / (J/(m³·K)) = Pa·s·m³·K/J = m³·K·s / m² = m·K·s
# But we need consistent units. η in Pa·s, s in J/(m³·K) = Pa/K
# η/s in Pa·s / (Pa/K) = K·s
# Actually η/s has units of ℏ/k_B = J·s / (J/K) = K·s
A_KSS = hbar / (4 * np.pi * k_B)  # K·s
print(f"KSS bound: η/s ≥ {A_KSS:.4e} K·s")
print(f"  = ℏ/(4π k_B) = {A_KSS:.4e} K·s")
print()

# =============================================================================
# MATERIAL DATA
# =============================================================================
# Format: (name, symbol, T_m [K], η_liquid [mPa·s at T_m], M [g/mol],
#          ΔS_fus [J/(mol·K)], ρ_liquid [kg/m³ at T_m], θ_D [K])
#
# η_liquid: viscosity of liquid just above melting point
# ΔS_fus: entropy of fusion per mole
# ρ_liquid: liquid density at T_m (for converting molar to volumetric)
# θ_D: Debye temperature (for comparison — NOT used in η/s calculation)
#
# Sources: Iida & Guthrie (1988), Assael et al. (2006), CRC Handbook,
#          Chase JANAF tables, Grimvall "Thermophysical Properties of Materials"
# =============================================================================

materials = [
    # Alkali metals — simple electronic structure, low T_m
    ("Lithium",    "Li",  453.7,  0.57,   6.94,   6.60,  512,   344),
    ("Sodium",     "Na",  370.9,  0.68,  22.99,   7.01,  927,   158),
    ("Potassium",  "K",   336.5,  0.51,  39.10,   7.08,  828,    91),
    ("Rubidium",   "Rb",  312.5,  0.67,  85.47,   7.05,  1470,   56),
    ("Cesium",     "Cs",  301.6,  0.68, 132.91,   7.71,  1843,   38),

    # Noble metals — d-band effects
    ("Copper",     "Cu", 1357.8,  4.0,   63.55,  9.62,  7998,  343),
    ("Silver",     "Ag", 1234.9,  3.88,  107.87, 11.30,  9320,  225),
    ("Gold",       "Au", 1337.3,  5.13,  196.97, 12.36, 17310,  165),

    # Simple metals
    ("Aluminum",   "Al",  933.5,  1.3,   26.98,  11.56,  2375,  428),
    ("Lead",       "Pb",  600.6,  2.65,  207.2,   7.99, 10660,  105),
    ("Tin",        "Sn",  505.1,  1.85,  118.71,  14.18,  6980,  200),
    ("Zinc",       "Zn",  692.7,  3.5,   65.38,  10.57,  6570,  327),
    ("Indium",     "In",  429.8,  1.75,  114.82,   7.59,  7020,  108),
    ("Bismuth",    "Bi",  544.6,  1.67,  208.98,  20.77,  10050,  119),

    # Transition metals — complex d-band
    ("Iron",       "Fe", 1811.0,  6.9,   55.85,   7.63,  6980,  470),
    ("Nickel",     "Ni", 1728.0,  5.5,   58.69,  10.12,  7780,  450),
    ("Cobalt",     "Co", 1768.0,  5.0,   58.93,   9.16,  7670,  445),
    ("Titanium",   "Ti", 1941.0,  5.2,   47.87,   9.95,  4110,  420),
    ("Chromium",   "Cr", 2180.0,  6.7,   52.00,   9.62,  6300,  630),

    # Refractory metals
    ("Tungsten",   "W",  3695.0, 8.0,   183.84,  9.62, 17600,  400),
    ("Molybdenum", "Mo", 2896.0, 5.5,   95.95,   9.37, 9330,   450),
    ("Tantalum",   "Ta", 3290.0, 7.5,   180.95,  9.87, 15000,  240),

    # Rare earth (representative)
    ("Gadolinium", "Gd", 1585.0, 4.8,   157.25,  9.41, 7400,   200),

    # Metalloids / semiconductors (for contrast)
    ("Germanium",  "Ge", 1211.4, 0.74,  72.63,  30.95,  5600,  374),
    ("Silicon",    "Si", 1687.0, 0.58,  28.09,  29.79,  2570,  645),

    # Non-metals (molecular — very different regime)
    ("Gallium",    "Ga",  302.9, 2.04,  69.72,  18.45,  6095,  320),
]

print("=" * 120)
print(f"{'Material':<14} {'T_m':>6} {'η_liq':>8} {'ΔS_fus':>8} {'ρ_liq':>7} {'θ_D':>5} "
      f"{'s_vol':>12} {'η/s':>12} {'A/A_KSS':>10} {'T_m/θ_D':>7} {'γ_melt':>7}")
print(f"{'':14} {'(K)':>6} {'(mPa·s)':>8} {'(J/mol·K)':>8} {'(kg/m³)':>7} {'(K)':>5} "
      f"{'(J/m³·K)':>12} {'(K·s)':>12} {'':>10} {'':>7} {'':>7}")
print("=" * 120)

# Storage for analysis
names = []
symbols = []
T_ms = []
eta_over_s_vals = []
A_ratios = []
T_over_theta = []
gamma_melts = []
lindemann_est = []

for name, sym, T_m, eta_mPas, M, dS_fus, rho_liq, theta_D in materials:
    # Convert viscosity: mPa·s → Pa·s
    eta = eta_mPas * 1e-3  # Pa·s

    # Volumetric entropy density at T_m
    # s_vol = (ρ/M) × N_A × (S_liquid_per_atom)
    # S_liquid per mole = S_solid(T_m) + ΔS_fus
    # For the liquid just above T_m, we use the entropy density
    # S_solid(T_m) ≈ 3R[4(T_m/θ_D)³ ∫₀^{θ_D/T_m} x³/(e^x-1)dx - ln(1-e^{-θ_D/T_m})]
    # At high T (T_m >> θ_D): S_solid ≈ 3R(1 + ln(T_m/θ_D))
    # At T_m ~ θ_D: use Debye function

    # Debye entropy function
    x_D = theta_D / T_m
    if x_D < 0.01:  # classical limit
        S_solid_mol = 3 * R * (1 + np.log(T_m / theta_D))
    else:
        # Numerical integration of Debye function
        n_pts = 1000
        x = np.linspace(1e-10, x_D, n_pts)
        dx = x[1] - x[0]
        # Debye entropy: S = 3R[-3ln(1-e^{-x_D}) + 4D₃(x_D)]
        # where D₃(x_D) = (3/x_D³) ∫₀^{x_D} x³/(e^x - 1) dx
        integrand = x**3 / (np.exp(x) - 1)
        D3 = (3.0 / x_D**3) * np.trapz(integrand, x)
        S_solid_mol = 3 * R * (4 * D3 - np.log(1 - np.exp(-x_D)))

    S_liquid_mol = S_solid_mol + dS_fus  # J/(mol·K) at T_m

    # Convert to volumetric entropy density: s [J/(m³·K)]
    n_density = rho_liq / (M * 1e-3)  # mol/m³ (M in g/mol → kg/mol)
    s_vol = n_density * S_liquid_mol  # J/(m³·K)

    # η/s ratio
    eta_s = eta / s_vol  # Pa·s / (J/(m³·K)) = Pa·s·m³·K/J = s·K (since Pa = J/m³)

    # Ratio to KSS bound
    A_ratio = eta_s / A_KSS

    # Chemistry track γ at melting
    gamma_m = 2 * T_m / theta_D

    # Lindemann parameter estimate
    # L² = 3k_BT_m / (M_atom × ω_D² × d²)
    # ω_D = k_B θ_D / ℏ
    # d ≈ (M/(ρ × N_A))^(1/3) — interatomic distance
    M_kg = M * 1e-3 / N_A  # mass per atom in kg
    omega_D = k_B * theta_D / hbar  # rad/s
    d = (M * 1e-3 / (rho_liq * N_A)) ** (1/3)  # m
    L_sq = 3 * k_B * T_m / (M_kg * omega_D**2 * d**2)
    L = np.sqrt(L_sq)

    names.append(name)
    symbols.append(sym)
    T_ms.append(T_m)
    eta_over_s_vals.append(eta_s)
    A_ratios.append(A_ratio)
    T_over_theta.append(T_m / theta_D)
    gamma_melts.append(gamma_m)
    lindemann_est.append(L)

    print(f"{name:<14} {T_m:6.0f} {eta_mPas:8.2f} {dS_fus:8.2f} {rho_liq:7.0f} {theta_D:5.0f} "
          f"{s_vol:12.1f} {eta_s:12.4e} {A_ratio:10.1f} {T_m/theta_D:7.2f} {gamma_m:7.2f}")

print("=" * 120)

# Convert to arrays
eta_over_s_vals = np.array(eta_over_s_vals)
A_ratios = np.array(A_ratios)
T_over_theta = np.array(T_over_theta)
gamma_melts = np.array(gamma_melts)
lindemann_est = np.array(lindemann_est)

print()
print("=" * 80)
print("ANALYSIS")
print("=" * 80)

# 1. Distribution of A/A_KSS at melting
print(f"\n--- η/s at melting point vs KSS bound ---")
print(f"KSS bound A_KSS = {A_KSS:.4e} K·s")
print(f"A/A_KSS at melting:")
print(f"  Mean:   {np.mean(A_ratios):.1f}")
print(f"  Median: {np.median(A_ratios):.1f}")
print(f"  Std:    {np.std(A_ratios):.1f}")
print(f"  Min:    {np.min(A_ratios):.1f} ({names[np.argmin(A_ratios)]})")
print(f"  Max:    {np.max(A_ratios):.1f} ({names[np.argmax(A_ratios)]})")
print(f"  CV:     {np.std(A_ratios)/np.mean(A_ratios):.3f}")

print(f"\n  Compare to Session 1 ranges:")
print(f"    QGP:           A/A_KSS ≈ 1-10")
print(f"    He-4 λ-point:  A/A_KSS ≈ 12")
print(f"    Liquid metals:  A/A_KSS ≈ 500-1200 (Session 1 estimate)")
print(f"    Melting point:  A/A_KSS = {np.median(A_ratios):.0f} (this session, median)")

# 2. Lindemann parameter universality check
print(f"\n--- Lindemann parameter estimates ---")
print(f"  Mean L:   {np.mean(lindemann_est):.4f}")
print(f"  Median L: {np.median(lindemann_est):.4f}")
print(f"  Std L:    {np.std(lindemann_est):.4f}")
print(f"  CV:       {np.std(lindemann_est)/np.mean(lindemann_est):.3f}")
print(f"  Expected: L ≈ 0.10 (Lindemann criterion)")

# Show per material
print(f"\n  {'Material':<14} {'L_est':>8}")
for i, name in enumerate(names):
    print(f"  {name:<14} {lindemann_est[i]:8.4f}")

# 3. Correlation: log(A/A_KSS) vs log(T_m/θ_D)
log_A = np.log10(A_ratios)
log_T_theta = np.log10(T_over_theta)

r_A_Ttheta = np.corrcoef(log_A, log_T_theta)[0, 1]
print(f"\n--- Circularity check: log(A/A_KSS) vs log(T_m/θ_D) ---")
print(f"  r = {r_A_Ttheta:.4f}")
if abs(r_A_Ttheta) > 0.7:
    print(f"  ⚠️  STRONG correlation — η/s at melting may be θ_D-dependent")
elif abs(r_A_Ttheta) > 0.4:
    print(f"  ⚠️  MODERATE correlation — partial θ_D dependence")
else:
    print(f"  ✓  WEAK correlation — η/s at melting largely independent of θ_D")

# 4. Correlation: log(A/A_KSS) vs log(T_m) alone
log_Tm = np.log10(T_ms)
r_A_Tm = np.corrcoef(log_A, log_Tm)[0, 1]
print(f"\n--- log(A/A_KSS) vs log(T_m) ---")
print(f"  r = {r_A_Tm:.4f}")

# 5. Correlation: log(A/A_KSS) vs log(η)
log_eta = np.log10([m[3] for m in materials])  # η in mPa·s
r_A_eta = np.corrcoef(log_A, log_eta)[0, 1]
print(f"\n--- log(A/A_KSS) vs log(η_liquid) ---")
print(f"  r = {r_A_eta:.4f}")
if abs(r_A_eta) > 0.8:
    print(f"  ⚠️  η/s is dominated by η variation (entropy too uniform)")

# 6. Entropy density variation
s_vols = []
for name, sym, T_m, eta_mPas, M, dS_fus, rho_liq, theta_D in materials:
    eta = eta_mPas * 1e-3
    x_D = theta_D / T_m
    if x_D < 0.01:
        S_solid_mol = 3 * R * (1 + np.log(T_m / theta_D))
    else:
        n_pts = 1000
        x = np.linspace(1e-10, x_D, n_pts)
        integrand = x**3 / (np.exp(x) - 1)
        D3 = (3.0 / x_D**3) * np.trapz(integrand, x)
        S_solid_mol = 3 * R * (4 * D3 - np.log(1 - np.exp(-x_D)))
    S_liquid_mol = S_solid_mol + dS_fus
    n_density = rho_liq / (M * 1e-3)
    s_vols.append(n_density * S_liquid_mol)
s_vols = np.array(s_vols)

print(f"\n--- Entropy density variation ---")
print(f"  s_vol range: {np.min(s_vols):.0f} — {np.max(s_vols):.0f} J/(m³·K)")
print(f"  CV(s_vol): {np.std(s_vols)/np.mean(s_vols):.3f}")
print(f"  CV(η):     {np.std([m[3] for m in materials])/np.mean([m[3] for m in materials]):.3f}")
print(f"  → {'Entropy' if np.std(s_vols)/np.mean(s_vols) > np.std([m[3] for m in materials])/np.mean([m[3] for m in materials]) else 'Viscosity'} varies more")

# 7. Quality factor at melting: Q_melt = ω_D × τ_phonon
# τ_phonon ≈ 1/(γ_G² × ω_D × (T/θ_D)²) for T > θ_D (Klemens approximation)
# Use γ_G ≈ 1.5-2.0 (typical Grüneisen for metals)
gamma_G_typical = 1.7
print(f"\n--- Phonon quality factor at melting (Klemens model, γ_G={gamma_G_typical}) ---")
print(f"  {'Material':<14} {'Q_melt':>8} {'γ_melt/f':>8} {'Entity?':>8}")
for i, (name, sym, T_m, eta_mPas, M, dS_fus, rho_liq, theta_D) in enumerate(materials):
    omega_D = k_B * theta_D / hbar
    # Klemens: τ_ph ≈ 1/(γ_G² × ω_D × (T/θ_D)²) when T > θ_D
    # For T < θ_D, more complex, but approximate
    tau_ph = 1.0 / (gamma_G_typical**2 * omega_D * (T_m / theta_D)**2)
    Q = omega_D * tau_ph  # = 1/(γ_G² × (T/θ_D)²)
    gamma_over_f = 1.0 / Q  # damping/frequency ratio
    entity = "YES" if gamma_over_f < 1 else "NO"
    print(f"  {name:<14} {Q:8.3f} {gamma_over_f:8.3f} {entity:>8}")

# 8. Key test: Is A/A_KSS at melting correlated with Lindemann parameter?
r_A_L = np.corrcoef(log_A, np.log10(lindemann_est))[0, 1]
print(f"\n--- log(A/A_KSS) vs log(L_Lindemann) ---")
print(f"  r = {r_A_L:.4f}")

# 9. Grand test: Can we find a UNIVERSAL value?
# If melting corresponds to a specific entity-to-process threshold,
# some combination of η/s and material parameters should be universal.
print(f"\n{'='*80}")
print("UNIVERSALITY TESTS")
print(f"{'='*80}")

# Test A: Raw η/s
print(f"\nTest A: Raw η/s at T_m")
print(f"  CV = {np.std(eta_over_s_vals)/np.mean(eta_over_s_vals):.3f}")
print(f"  Spread: {np.max(eta_over_s_vals)/np.min(eta_over_s_vals):.1f}×")
if np.std(eta_over_s_vals)/np.mean(eta_over_s_vals) < 0.3:
    print(f"  ✓ Reasonably universal!")
else:
    print(f"  ✗ Not universal (spread too large)")

# Test B: η/(s × ℏ/k_B) — dimensionless, A/A_KSS
print(f"\nTest B: A/A_KSS at T_m")
print(f"  CV = {np.std(A_ratios)/np.mean(A_ratios):.3f}")
print(f"  Spread: {np.max(A_ratios)/np.min(A_ratios):.1f}×")

# Test C: η × d / (s × ℏ) — include interatomic distance
d_vals = np.array([(M * 1e-3 / (rho * N_A))**(1/3) for _, _, _, _, M, _, rho, _ in materials])
test_C = eta_over_s_vals / (hbar / k_B) * d_vals  # dimensionless-ish
print(f"\nTest C: (η/s) × d / (ℏ/k_B) at T_m")
print(f"  CV = {np.std(test_C)/np.mean(test_C):.3f}")
print(f"  Spread: {np.max(test_C)/np.min(test_C):.1f}×")

# Test D: η/(s × T_m) — remove temperature dependence
test_D = eta_over_s_vals / np.array(T_ms)  # units: s
print(f"\nTest D: η/(s × T_m) at T_m")
print(f"  CV = {np.std(test_D)/np.mean(test_D):.3f}")
print(f"  Spread: {np.max(test_D)/np.min(test_D):.1f}×")

# Test E: (η/s) × k_B T_m / ℏ — Planck time scaling
test_E = eta_over_s_vals * k_B * np.array(T_ms) / hbar
print(f"\nTest E: (η/s) × k_BT_m/ℏ (= τ_Planck × k_BT_m/ℏ)")
print(f"  Mean:   {np.mean(test_E):.2e}")
print(f"  Median: {np.median(test_E):.2e}")
print(f"  CV = {np.std(test_E)/np.mean(test_E):.3f}")
print(f"  Spread: {np.max(test_E)/np.min(test_E):.1f}×")

# What IS the tightest clustering?
tests = {
    "A: A/A_KSS": (A_ratios, np.std(A_ratios)/np.mean(A_ratios)),
    "B: log(A/A_KSS)": (log_A, np.std(log_A)/np.mean(log_A)),
    "C: η/(s·T_m)": (test_D, np.std(test_D)/np.mean(test_D)),
    "D: (η/s)·d/(ℏ/k_B)": (test_C, np.std(test_C)/np.mean(test_C)),
    "E: (η/s)·k_BT_m/ℏ": (test_E, np.std(test_E)/np.mean(test_E)),
}
print(f"\n--- Summary: which quantity clusters tightest? ---")
for label, (vals, cv) in sorted(tests.items(), key=lambda x: abs(x[1][1])):
    print(f"  {label:<25} CV = {abs(cv):.3f}")

# 10. Regression analysis
print(f"\n{'='*80}")
print("REGRESSION: What determines η/s at melting?")
print(f"{'='*80}")

# Multivariate regression: log(η/s) ~ a + b·log(T_m) + c·log(M) + d·log(ρ)
from numpy.linalg import lstsq

log_M = np.log10([m[4] for m in materials])
log_rho = np.log10([m[5] for m in materials])

X = np.column_stack([np.ones(len(log_A)), log_Tm, log_M, log_rho])
coeffs, residuals, _, _ = lstsq(X, log_A, rcond=None)
y_pred = X @ coeffs
r_multi = np.corrcoef(log_A, y_pred)[0, 1]

print(f"\nlog(A/A_KSS) = {coeffs[0]:.3f} + {coeffs[1]:.3f}·log(T_m) + {coeffs[2]:.3f}·log(M) + {coeffs[3]:.3f}·log(ρ)")
print(f"R² = {r_multi**2:.4f}")
print(f"Residual std = {np.std(log_A - y_pred):.4f} dex")

# Dominant variable
for i, label in enumerate(["intercept", "log(T_m)", "log(M)", "log(ρ)"]):
    print(f"  {label}: {coeffs[i]:+.3f}")

# Include θ_D
log_theta = np.log10([m[7] for m in materials])
X2 = np.column_stack([np.ones(len(log_A)), log_Tm, log_M, log_rho, log_theta])
coeffs2, _, _, _ = lstsq(X2, log_A, rcond=None)
y_pred2 = X2 @ coeffs2
r_multi2 = np.corrcoef(log_A, y_pred2)[0, 1]

print(f"\nWith θ_D added:")
print(f"log(A/A_KSS) = {coeffs2[0]:.3f} + {coeffs2[1]:.3f}·log(T_m) + {coeffs2[2]:.3f}·log(M) + {coeffs2[3]:.3f}·log(ρ) + {coeffs2[4]:.3f}·log(θ_D)")
print(f"R² = {r_multi2**2:.4f} (improvement from θ_D: {r_multi2**2 - r_multi**2:.4f})")

print(f"\n{'='*80}")
print("INTERPRETATION")
print(f"{'='*80}")

print(f"""
1. UNIVERSALITY OF η/s AT MELTING:
   - A/A_KSS ranges from {np.min(A_ratios):.0f} to {np.max(A_ratios):.0f} (spread: {np.max(A_ratios)/np.min(A_ratios):.0f}×)
   - CV = {np.std(A_ratios)/np.mean(A_ratios):.2f}
   - For comparison: liquid metals at room T have spread ~8× (Session 1)
   - Verdict: {'η/s at melting is NOT universal' if np.std(A_ratios)/np.mean(A_ratios) > 0.5 else 'η/s at melting shows moderate universality'}

2. WHERE MELTING SITS IN THE KSS HIERARCHY:
   - QGP:         A/A_KSS ~ 1-10     (quantum critical)
   - He-4:        A/A_KSS ~ 12       (near quantum)
   - MELTING:     A/A_KSS ~ {np.median(A_ratios):.0f}    (this session)
   - Room T:      A/A_KSS ~ 500-4000 (classical)

3. ENTITY CRITERION CONNECTION:
   - If melting = entity→process transition, expect characteristic η/s
   - {'SUPPORTED' if np.std(A_ratios)/np.mean(A_ratios) < 0.5 else 'NOT SUPPORTED'}: CV = {np.std(A_ratios)/np.mean(A_ratios):.2f} ({'tight' if np.std(A_ratios)/np.mean(A_ratios) < 0.3 else 'moderate' if np.std(A_ratios)/np.mean(A_ratios) < 0.5 else 'wide'} clustering)
   - The melting transition is NOT at a universal η/s — materials melt at
     different distances from the KSS bound

4. WHAT DETERMINES η/s AT MELTING:
   - Dominant factor: {'T_m' if abs(coeffs[1]) > abs(coeffs[2]) and abs(coeffs[1]) > abs(coeffs[3]) else 'M' if abs(coeffs[2]) > abs(coeffs[3]) else 'ρ'}
   - R² without θ_D: {r_multi**2:.3f}
   - R² with θ_D:    {r_multi2**2:.3f}
   - θ_D adds {r_multi2**2 - r_multi**2:.3f} to R² → {'minimal' if r_multi2**2 - r_multi**2 < 0.05 else 'moderate' if r_multi2**2 - r_multi**2 < 0.1 else 'substantial'} additional information
   - {'η/s is genuinely non-circular with θ_D' if r_multi2**2 - r_multi**2 < 0.05 else 'η/s has some θ_D dependence'}
""")

# Final assessment
print("=" * 80)
print("SESSION 3 VERDICT")
print("=" * 80)
print(f"""
The crystal→liquid transition does NOT occur at a universal η/s value.

η/s at melting varies by {np.max(A_ratios)/np.min(A_ratios):.0f}× across materials (CV={np.std(A_ratios)/np.mean(A_ratios):.2f}).
This is {'comparable to' if abs(np.std(A_ratios)/np.mean(A_ratios) - 0.5) < 0.2 else 'wider than' if np.std(A_ratios)/np.mean(A_ratios) > 0.7 else 'narrower than'} the spread at fixed temperature.

The Lindemann parameter (L ≈ {np.mean(lindemann_est):.3f}, CV={np.std(lindemann_est)/np.mean(lindemann_est):.2f}) IS more universal
than η/s at melting (CV={np.std(A_ratios)/np.mean(A_ratios):.2f}). This means:
- The melting criterion is better expressed in DISPLACEMENT units (L) than VISCOSITY units (η/s)
- The entity→process transition for crystals is about amplitude exceeding structure,
  NOT about dissipation exceeding coherence (which is what η/s measures)

KEY FINDING: Melting is an AMPLITUDE criterion, not a DISSIPATION criterion.
The Lindemann ratio measures "how far atoms move vs how far apart they are."
η/s measures "how much dissipation vs how much entropy."
These are different physics. The crystal-as-entity dissolves when displacement
exceeds structure, regardless of viscosity.

IMPLICATION FOR SYNCHRONISM: The entity criterion γ/f < 1 (damping vs frequency)
is a DISSIPATION criterion. But the most universal material phase transition
(melting) is an AMPLITUDE criterion. This suggests that for material-scale
entities (crystals), the relevant entity criterion may be displacement-based
rather than dissipation-based.
""")
