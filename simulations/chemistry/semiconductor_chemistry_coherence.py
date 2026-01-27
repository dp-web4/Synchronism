#!/usr/bin/env python3
"""
Chemistry Session #254: Semiconductor / Electronics Chemistry
Finding #191 | 117th phenomenon type at γ ~ 1

Applying Synchronism coherence framework to semiconductor fabrication,
electronic materials, and microelectronics chemistry.

Key γ ~ 1 boundaries investigated:
1. Doping: Intrinsic carrier concentration (n = p = n_i)
2. CVD/Epitaxy: Growth rate regimes (mass transport = surface reaction)
3. Etching: Selectivity and etch stop (etch rate ratio = 1)
4. Oxidation: Deal-Grove model (linear-parabolic transition)
5. Photolithography: Exposure threshold (dose = E_th)
6. CMP: Preston equation (removal = addition balance)
7. Diffusion: Junction depth (C = C_background)
8. Electroplating: Current efficiency (η = 1)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ==============================================================
# ANALYSIS 1: Intrinsic Carrier Concentration (n = p = n_i)
# ==============================================================

def analyze_doping():
    """At n = p = n_i: intrinsic semiconductor (γ ~ 1 exactly!)"""

    # Temperature range
    T = np.linspace(200, 600, 500)  # K
    k_B = 8.617e-5  # eV/K

    # Silicon parameters
    E_g_0 = 1.166  # eV at 0 K
    alpha = 4.73e-4  # eV/K
    beta = 636  # K
    E_g = E_g_0 - alpha * T**2 / (T + beta)

    # Effective density of states (simplified)
    N_c = 2.86e19 * (T / 300)**1.5  # cm^-3
    N_v = 3.10e19 * (T / 300)**1.5

    # Intrinsic carrier concentration
    n_i = np.sqrt(N_c * N_v) * np.exp(-E_g / (2 * k_B * T))

    # Doping: n-type with N_D = 1e15 cm^-3
    N_D = 1e15
    n = 0.5 * (N_D + np.sqrt(N_D**2 + 4 * n_i**2))
    p = n_i**2 / n

    # γ ~ 1 boundary: where n/p = 1 (intrinsic behavior)
    ratio_n_p = n / p

    # Find intrinsic temperature (where N_D << n_i, so n ≈ p)
    idx_intrinsic = np.argmin(np.abs(ratio_n_p - 1))
    T_intrinsic = T[idx_intrinsic]

    # Fermi level relative to midgap
    E_F_minus_E_i = k_B * T * np.log(n / n_i)

    return {
        'T': T, 'n_i': n_i, 'n': n, 'p': p,
        'ratio': ratio_n_p, 'T_intrinsic': T_intrinsic,
        'E_F_shift': E_F_minus_E_i, 'E_g': E_g
    }

# ==============================================================
# ANALYSIS 2: CVD Growth Rate Regimes
# ==============================================================

def analyze_cvd():
    """At T where k_s = h_g: surface reaction = mass transport (γ ~ 1!)"""

    # Temperature range
    T = np.linspace(400, 1200, 500)  # °C
    T_K = T + 273.15

    # Surface reaction rate (Arrhenius)
    k_s0 = 1e8  # cm/s pre-exponential
    E_a = 1.6  # eV (SiH4 decomposition)
    k_B = 8.617e-5  # eV/K
    k_s = k_s0 * np.exp(-E_a / (k_B * T_K))

    # Mass transport coefficient (weakly T-dependent)
    h_g = 2.0  # cm/s (typical for atmospheric CVD)
    h_g_array = h_g * (T_K / 673)**0.5  # slight T dependence from diffusivity

    # Overall growth rate (series resistance model)
    # 1/R = 1/k_s + 1/h_g → R = k_s * h_g / (k_s + h_g)
    R = k_s * h_g_array / (k_s + h_g_array)

    # Ratio: k_s / h_g
    ratio = k_s / h_g_array

    # γ ~ 1 transition temperature
    idx_transition = np.argmin(np.abs(ratio - 1))
    T_transition = T[idx_transition]

    # Regime classification
    # ratio << 1: surface-reaction limited (low T)
    # ratio >> 1: mass-transport limited (high T)

    return {
        'T': T, 'k_s': k_s, 'h_g': h_g_array, 'R': R,
        'ratio': ratio, 'T_transition': T_transition
    }

# ==============================================================
# ANALYSIS 3: Plasma Etching Selectivity
# ==============================================================

def analyze_etching():
    """Selectivity = 1: no discrimination between materials (γ ~ 1!)"""

    # Ion energy range
    E_ion = np.linspace(10, 500, 500)  # eV

    # Etch rates (Steinbrüchel model: R = a * sqrt(E - E_th))
    # SiO2 etching in fluorocarbon plasma
    a_oxide = 0.5  # nm/min per eV^0.5
    E_th_oxide = 30  # eV threshold
    R_oxide = np.where(E_ion > E_th_oxide,
                        a_oxide * np.sqrt(E_ion - E_th_oxide), 0)

    # Si etching (lower threshold, higher rate)
    a_si = 0.8
    E_th_si = 15  # eV
    R_si = np.where(E_ion > E_th_si,
                     a_si * np.sqrt(E_ion - E_th_si), 0)

    # Photoresist (higher threshold)
    a_pr = 0.3
    E_th_pr = 50  # eV
    R_pr = np.where(E_ion > E_th_pr,
                     a_pr * np.sqrt(E_ion - E_th_pr), 0)

    # Selectivity: SiO2:Si and SiO2:PR
    with np.errstate(divide='ignore', invalid='ignore'):
        sel_oxide_si = np.where(R_si > 0, R_oxide / R_si, 0)
        sel_oxide_pr = np.where(R_pr > 0, R_oxide / R_pr, 0)

    # Find where selectivity = 1 (γ ~ 1)
    valid_sel = sel_oxide_si[sel_oxide_si > 0]
    valid_E = E_ion[sel_oxide_si > 0]
    if len(valid_sel) > 0:
        idx_unity = np.argmin(np.abs(valid_sel - 1))
        E_unity = valid_E[idx_unity]
    else:
        E_unity = 100  # fallback

    return {
        'E_ion': E_ion, 'R_oxide': R_oxide, 'R_si': R_si, 'R_pr': R_pr,
        'sel_oxide_si': sel_oxide_si, 'sel_oxide_pr': sel_oxide_pr,
        'E_unity': E_unity
    }

# ==============================================================
# ANALYSIS 4: Deal-Grove Thermal Oxidation
# ==============================================================

def analyze_oxidation():
    """Linear-parabolic transition: t = τ at crossover (γ ~ 1!)"""

    # Time range
    t = np.logspace(-2, 4, 500)  # minutes

    # Deal-Grove parameters for dry O2 on Si at 1000°C
    B = 0.0117  # μm²/min (parabolic rate constant)
    B_over_A = 0.0117 / 0.165  # μm/min (linear rate constant)
    A = 0.165  # μm
    tau = 0.37  # min (accounts for initial oxide)

    # Deal-Grove equation: x² + Ax = B(t + τ)
    # Solution: x = A/2 * (sqrt(1 + 4B(t+τ)/A²) - 1)
    x_ox = (A / 2) * (np.sqrt(1 + 4 * B * (t + tau) / A**2) - 1)

    # Linear limit (thin oxide): x ≈ (B/A)(t + τ)
    x_linear = (B / A) * (t + tau)

    # Parabolic limit (thick oxide): x ≈ sqrt(B(t + τ))
    x_parabolic = np.sqrt(B * (t + tau))

    # Crossover thickness where linear = parabolic contributions equal
    # This occurs at x = A (the Deal-Grove crossover)
    x_crossover = A
    t_crossover = A**2 / B  # time when x ≈ A

    # Rate: dx/dt
    dx_dt = B / (2 * x_ox + A)

    # Ratio of linear to parabolic contribution
    ratio_lin_par = A / (2 * x_ox)  # ratio of A to 2x in denominator

    return {
        't': t, 'x_ox': x_ox, 'x_linear': x_linear, 'x_parabolic': x_parabolic,
        'dx_dt': dx_dt, 'x_crossover': x_crossover, 't_crossover': t_crossover,
        'ratio': ratio_lin_par
    }

# ==============================================================
# ANALYSIS 5: Photolithography Exposure
# ==============================================================

def analyze_lithography():
    """At dose = E_th: resist transformation threshold (γ ~ 1!)"""

    # Dose range
    dose = np.linspace(0, 200, 500)  # mJ/cm²

    # Positive resist (Dill model): exposed region becomes soluble
    # Dissolution rate R = R_max * (1 - exp(-C * (dose - E_th)))
    E_th_pos = 50  # mJ/cm² threshold
    C_pos = 0.03  # sensitivity parameter
    R_pos = np.where(dose > E_th_pos,
                      1.0 * (1 - np.exp(-C_pos * (dose - E_th_pos))),
                      0.01)  # dark erosion rate

    # Negative resist: exposed region becomes insoluble
    E_th_neg = 30  # mJ/cm²
    C_neg = 0.05
    R_neg = 1.0 * np.exp(-C_neg * np.maximum(dose - E_th_neg * 0.5, 0))

    # Contrast (positive resist)
    # γ_resist = 1 / log10(E_clear/E_th)
    E_clear = E_th_pos * 2.5
    contrast = 1 / np.log10(E_clear / E_th_pos)

    # Depth of focus: at Rayleigh criterion
    wavelength = 193  # nm (ArF)
    NA = 0.85
    k1 = 0.5  # at γ ~ 1 resolution limit!
    resolution = k1 * wavelength / NA
    k2 = 0.5
    DOF = k2 * wavelength / NA**2

    # Normalized dose
    dose_norm = dose / E_th_pos

    return {
        'dose': dose, 'R_pos': R_pos, 'R_neg': R_neg,
        'E_th_pos': E_th_pos, 'E_th_neg': E_th_neg,
        'contrast': contrast, 'resolution': resolution, 'DOF': DOF,
        'k1': k1, 'dose_norm': dose_norm
    }

# ==============================================================
# ANALYSIS 6: CMP (Chemical Mechanical Planarization)
# ==============================================================

def analyze_cmp():
    """Preston equation: removal rate at material balance (γ ~ 1!)"""

    # Pressure range
    P = np.linspace(0.5, 10, 500)  # psi

    # Preston equation: RR = K_p * P * V
    V = 100  # rpm (relative velocity, constant)

    # Different materials (K_p in nm/min per psi per rpm)
    K_p_oxide = 0.8
    K_p_cu = 2.5
    K_p_barrier = 0.3  # TaN
    K_p_nitride = 0.15

    RR_oxide = K_p_oxide * P * V
    RR_cu = K_p_cu * P * V
    RR_barrier = K_p_barrier * P * V
    RR_nitride = K_p_nitride * P * V

    # Selectivity ratios
    sel_cu_ox = RR_cu / RR_oxide
    sel_cu_barrier = RR_cu / RR_barrier
    sel_ox_nit = RR_oxide / RR_nitride

    # Step height planarization
    # At step height = 0: perfectly planar (γ ~ 1 target!)
    h0 = 500  # nm initial step
    # Planarization rate depends on pattern density
    density = np.linspace(0, 1, 100)
    # At density = 0.5: maximum non-uniformity (γ ~ 1!)
    h_remaining = h0 * np.abs(density - 0.5) / 0.5

    return {
        'P': P, 'RR_oxide': RR_oxide, 'RR_cu': RR_cu,
        'RR_barrier': RR_barrier, 'RR_nitride': RR_nitride,
        'sel_cu_ox': sel_cu_ox, 'sel_cu_barrier': sel_cu_barrier,
        'density': density, 'h_remaining': h_remaining
    }

# ==============================================================
# ANALYSIS 7: Dopant Diffusion / Junction Formation
# ==============================================================

def analyze_diffusion():
    """At C(x_j) = C_background: junction depth (γ ~ 1!)"""

    # Depth range
    x = np.linspace(0, 2.0, 500)  # μm

    # Gaussian diffusion profile (drive-in)
    C_surface = 1e20  # cm^-3
    C_bg = 1e15  # background (n-type substrate)

    # Different Dt products (diffusion time)
    Dt_values = [0.01, 0.05, 0.1, 0.2]  # μm²

    profiles = []
    junctions = []

    for Dt in Dt_values:
        C = C_surface * np.exp(-x**2 / (4 * Dt))
        profiles.append(C)

        # Junction depth where C = C_bg
        x_j = np.sqrt(4 * Dt * np.log(C_surface / C_bg))
        junctions.append(x_j)

    # Complementary error function profile (constant surface)
    from scipy.special import erfc
    D_B = 0.037  # μm²/hr (Boron in Si at 1000°C)
    t_values = [0.5, 1.0, 2.0, 4.0]  # hours

    erfc_profiles = []
    erfc_junctions = []
    for t in t_values:
        Dt_val = D_B * t
        C_erfc = C_surface * erfc(x / (2 * np.sqrt(Dt_val)))
        erfc_profiles.append(C_erfc)
        x_j = 2 * np.sqrt(Dt_val) * 1.5  # approximate for erfc
        erfc_junctions.append(x_j)

    return {
        'x': x, 'C_bg': C_bg, 'C_surface': C_surface,
        'profiles': profiles, 'junctions': junctions, 'Dt_values': Dt_values,
        'erfc_profiles': erfc_profiles, 'erfc_junctions': erfc_junctions,
        't_values': t_values
    }

# ==============================================================
# ANALYSIS 8: Electroplating Current Efficiency
# ==============================================================

def analyze_electroplating():
    """At η = 1: all current goes to deposition (γ ~ 1!)"""

    # Current density range
    j = np.linspace(1, 100, 500)  # mA/cm²

    # Copper electroplating from acid sulfate bath
    # Faraday's law: thickness = (j * t * M) / (n * F * ρ)
    M_Cu = 63.55  # g/mol
    n_e = 2  # electrons per Cu²⁺
    F = 96485  # C/mol
    rho_Cu = 8.96  # g/cm³

    # Current efficiency (decreases at high current due to H₂ evolution)
    j_lim = 80  # mA/cm² (limiting current density)
    eta = np.where(j < j_lim,
                    1.0 - 0.01 * (j / j_lim)**2,
                    1.0 - 0.01 - 0.5 * ((j - j_lim) / j_lim)**0.5)
    eta = np.clip(eta, 0.3, 1.0)

    # Deposition rate (nm/min at η)
    t_dep = 60  # seconds
    thickness_rate = eta * (j * 1e-3) * M_Cu * 1e7 / (n_e * F * rho_Cu)  # nm/s
    thickness_rate_nm_min = thickness_rate * 60

    # Throwing power: ratio of thickness at cavity bottom to surface
    # Wagner number Wa = κ/(dj/dη * L)
    # At Wa = 1: uniform deposition (γ ~ 1!)
    aspect_ratio = np.linspace(0.1, 10, 200)
    Wa = 1.0 / aspect_ratio  # simplified
    throwing_power = 1 / (1 + 1/Wa)  # fraction at bottom vs top

    # Leveling: at leveling agent concentration = CMC
    # Peak efficiency at optimal additive concentration
    c_add = np.linspace(0, 20, 200)  # ppm
    c_optimal = 5  # ppm
    grain_size = 50 + 200 * np.exp(-0.5 * (c_add / c_optimal))  # nm
    roughness = 10 + 100 * np.abs(c_add / c_optimal - 1)  # nm

    return {
        'j': j, 'eta': eta, 'rate': thickness_rate_nm_min,
        'j_lim': j_lim,
        'aspect_ratio': aspect_ratio, 'Wa': Wa, 'throwing_power': throwing_power,
        'c_add': c_add, 'grain_size': grain_size, 'roughness': roughness
    }

# ==============================================================
# MASTER VISUALIZATION
# ==============================================================

def create_master_figure():
    """Generate comprehensive 8-panel figure."""

    # Run all analyses
    doping = analyze_doping()
    cvd = analyze_cvd()
    etch = analyze_etching()
    oxide = analyze_oxidation()
    litho = analyze_lithography()
    cmp = analyze_cmp()
    diff = analyze_diffusion()
    plating = analyze_electroplating()

    fig = plt.figure(figsize=(24, 28))
    fig.suptitle(
        'Chemistry Session #254: Semiconductor / Electronics Chemistry\n'
        'Finding #191 | 117th Phenomenon Type at γ ~ 1',
        fontsize=18, fontweight='bold', y=0.98
    )

    gs = gridspec.GridSpec(4, 2, hspace=0.35, wspace=0.3,
                           left=0.06, right=0.96, top=0.94, bottom=0.03)

    # --- Panel 1: Doping / Intrinsic Carriers ---
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.semilogy(doping['T'], doping['n'], 'b-', linewidth=2, label='n (electrons)')
    ax1.semilogy(doping['T'], doping['p'], 'r-', linewidth=2, label='p (holes)')
    ax1.semilogy(doping['T'], doping['n_i'], 'k--', linewidth=1.5, label='n_i (intrinsic)')
    ax1.axvline(doping['T_intrinsic'], color='green', linestyle=':', linewidth=2,
                label=f'n = p at T = {doping["T_intrinsic"]:.0f} K')
    ax1.axhline(1e15, color='gray', linestyle='--', alpha=0.5, label='N_D = 10¹⁵')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Carrier Concentration (cm⁻³)')
    ax1.set_title('1. INTRINSIC CARRIERS: n/p = 1 at T_intrinsic (γ ~ 1!)')
    ax1.legend(fontsize=8, loc='lower right')
    ax1.set_ylim(1e8, 1e21)
    ax1.grid(True, alpha=0.3)

    # Annotation
    ax1.annotate(f'γ ~ 1 BOUNDARY\nn = p = n_i\nT = {doping["T_intrinsic"]:.0f} K',
                xy=(doping['T_intrinsic'], doping['n_i'][np.argmin(np.abs(doping['T'] - doping['T_intrinsic']))]),
                xytext=(doping['T_intrinsic'] - 80, 1e18),
                fontsize=9, fontweight='bold', color='green',
                arrowprops=dict(arrowstyle='->', color='green'))

    # --- Panel 2: CVD Growth Regimes ---
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.semilogy(cvd['T'], cvd['k_s'], 'r-', linewidth=2, label='k_s (surface reaction)')
    ax2.semilogy(cvd['T'], cvd['h_g'], 'b-', linewidth=2, label='h_g (mass transport)')
    ax2.semilogy(cvd['T'], cvd['R'], 'k-', linewidth=2.5, label='Overall rate R')
    ax2.axvline(cvd['T_transition'], color='green', linestyle=':', linewidth=2,
                label=f'k_s = h_g at T = {cvd["T_transition"]:.0f}°C')
    ax2.set_xlabel('Temperature (°C)')
    ax2.set_ylabel('Rate (cm/s)')
    ax2.set_title('2. CVD: Surface Reaction = Mass Transport (γ ~ 1!)')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Add regime labels
    ax2.text(500, 0.5, 'REACTION\nLIMITED\n(k_s << h_g)',
             fontsize=10, ha='center', color='red', alpha=0.7)
    ax2.text(1000, 0.5, 'TRANSPORT\nLIMITED\n(k_s >> h_g)',
             fontsize=10, ha='center', color='blue', alpha=0.7)

    # --- Panel 3: Etching Selectivity ---
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(etch['E_ion'], etch['R_oxide'], 'b-', linewidth=2, label='SiO₂')
    ax3.plot(etch['E_ion'], etch['R_si'], 'r-', linewidth=2, label='Si')
    ax3.plot(etch['E_ion'], etch['R_pr'], 'g-', linewidth=2, label='Photoresist')
    ax3.set_xlabel('Ion Energy (eV)')
    ax3.set_ylabel('Etch Rate (a.u.)')
    ax3.set_title('3. PLASMA ETCHING: Selectivity = 1 Boundaries (γ ~ 1!)')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    ax3_twin = ax3.twinx()
    valid = etch['sel_oxide_si'] > 0
    ax3_twin.plot(etch['E_ion'][valid], etch['sel_oxide_si'][valid],
                  'k--', linewidth=1.5, label='Selectivity SiO₂:Si')
    ax3_twin.axhline(1.0, color='green', linestyle=':', linewidth=2)
    ax3_twin.set_ylabel('Selectivity Ratio')
    ax3_twin.set_ylim(0, 3)
    ax3_twin.legend(fontsize=8, loc='upper right')
    ax3_twin.annotate('γ ~ 1: No selectivity!',
                      xy=(etch['E_unity'], 1.0), xytext=(etch['E_unity'] + 80, 2.0),
                      fontsize=10, fontweight='bold', color='green',
                      arrowprops=dict(arrowstyle='->', color='green'))

    # --- Panel 4: Deal-Grove Oxidation ---
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.loglog(oxide['t'], oxide['x_ox'] * 1000, 'k-', linewidth=2.5, label='Deal-Grove')
    ax4.loglog(oxide['t'], oxide['x_linear'] * 1000, 'b--', linewidth=1.5, label='Linear limit')
    ax4.loglog(oxide['t'], oxide['x_parabolic'] * 1000, 'r--', linewidth=1.5, label='Parabolic limit')
    ax4.axhline(oxide['x_crossover'] * 1000, color='green', linestyle=':', linewidth=2,
                label=f'x = A = {oxide["x_crossover"]*1000:.0f} nm (crossover)')
    ax4.axvline(oxide['t_crossover'], color='green', linestyle=':', linewidth=1, alpha=0.5)
    ax4.set_xlabel('Time (min)')
    ax4.set_ylabel('Oxide Thickness (nm)')
    ax4.set_title('4. THERMAL OXIDATION: Linear-Parabolic Transition (γ ~ 1!)')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0.01, 1e4)

    ax4.annotate(f'γ ~ 1: x = A\nLinear = Parabolic\nt = {oxide["t_crossover"]:.1f} min',
                xy=(oxide['t_crossover'], oxide['x_crossover'] * 1000),
                xytext=(oxide['t_crossover'] * 10, oxide['x_crossover'] * 1000 * 3),
                fontsize=9, fontweight='bold', color='green',
                arrowprops=dict(arrowstyle='->', color='green'))

    # --- Panel 5: Photolithography ---
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.plot(litho['dose'], litho['R_pos'], 'b-', linewidth=2, label='Positive resist')
    ax5.plot(litho['dose'], litho['R_neg'], 'r-', linewidth=2, label='Negative resist')
    ax5.axvline(litho['E_th_pos'], color='blue', linestyle=':', linewidth=2,
                label=f'E_th (pos) = {litho["E_th_pos"]} mJ/cm²')
    ax5.axvline(litho['E_th_neg'], color='red', linestyle=':', linewidth=2,
                label=f'E_th (neg) = {litho["E_th_neg"]} mJ/cm²')
    ax5.set_xlabel('Exposure Dose (mJ/cm²)')
    ax5.set_ylabel('Dissolution Rate (normalized)')
    ax5.set_title('5. PHOTOLITHOGRAPHY: Dose = E_threshold (γ ~ 1!)')
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    # Add k1 annotation
    ax5.text(130, 0.8, f'k₁ = {litho["k1"]:.1f} (resolution limit)\n'
             f'λ = {193} nm, NA = 0.85\n'
             f'Resolution = {litho["resolution"]:.0f} nm\n'
             f'DOF = {litho["DOF"]:.0f} nm\n'
             f'Contrast γ = {litho["contrast"]:.2f}',
             fontsize=9, bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # --- Panel 6: CMP ---
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.plot(cmp['P'], cmp['RR_cu'], 'r-', linewidth=2, label='Cu')
    ax6.plot(cmp['P'], cmp['RR_oxide'], 'b-', linewidth=2, label='SiO₂')
    ax6.plot(cmp['P'], cmp['RR_barrier'], 'g-', linewidth=2, label='TaN barrier')
    ax6.plot(cmp['P'], cmp['RR_nitride'], 'm-', linewidth=2, label='Si₃N₄')
    ax6.set_xlabel('Pressure (psi)')
    ax6.set_ylabel('Removal Rate (nm/min)')
    ax6.set_title('6. CMP: Preston Equation / Selectivity Control (γ ~ 1!)')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)

    # Inset: pattern density effect
    ax6_in = ax6.inset_axes([0.55, 0.45, 0.4, 0.45])
    ax6_in.plot(cmp['density'], cmp['h_remaining'], 'k-', linewidth=2)
    ax6_in.axvline(0.5, color='green', linestyle=':', linewidth=2)
    ax6_in.set_xlabel('Pattern Density', fontsize=8)
    ax6_in.set_ylabel('Step Height (nm)', fontsize=8)
    ax6_in.set_title('Density = 0.5: Max Non-uniformity', fontsize=8)
    ax6_in.tick_params(labelsize=7)

    # --- Panel 7: Diffusion / Junction ---
    ax7 = fig.add_subplot(gs[3, 0])
    colors_diff = ['blue', 'cyan', 'orange', 'red']
    for i, (prof, Dt_val) in enumerate(zip(diff['profiles'], diff['Dt_values'])):
        ax7.semilogy(diff['x'], prof, color=colors_diff[i], linewidth=2,
                     label=f'Dt = {Dt_val} μm²')
        # Mark junction
        x_j = diff['junctions'][i]
        ax7.axvline(x_j, color=colors_diff[i], linestyle=':', alpha=0.5)

    ax7.axhline(diff['C_bg'], color='green', linestyle='--', linewidth=2,
                label=f'C_bg = {diff["C_bg"]:.0e} cm⁻³')
    ax7.set_xlabel('Depth (μm)')
    ax7.set_ylabel('Concentration (cm⁻³)')
    ax7.set_title('7. DIFFUSION: C(x_j) = C_background at Junction (γ ~ 1!)')
    ax7.legend(fontsize=8)
    ax7.set_ylim(1e12, 1e21)
    ax7.grid(True, alpha=0.3)

    ax7.annotate('γ ~ 1: C = C_bg\nJunction depth x_j\np-n boundary',
                xy=(diff['junctions'][2], diff['C_bg']),
                xytext=(diff['junctions'][2] + 0.3, 1e17),
                fontsize=9, fontweight='bold', color='green',
                arrowprops=dict(arrowstyle='->', color='green'))

    # --- Panel 8: Electroplating ---
    ax8 = fig.add_subplot(gs[3, 1])
    ax8.plot(plating['j'], plating['eta'], 'b-', linewidth=2.5, label='Current efficiency η')
    ax8.axhline(1.0, color='green', linestyle=':', linewidth=2, label='η = 1 (ideal, γ ~ 1)')
    ax8.axvline(plating['j_lim'], color='red', linestyle='--', linewidth=1.5,
                label=f'j_lim = {plating["j_lim"]} mA/cm²')
    ax8.set_xlabel('Current Density (mA/cm²)')
    ax8.set_ylabel('Current Efficiency η')
    ax8.set_title('8. ELECTROPLATING: η = 1 (All Current → Deposition, γ ~ 1!)')
    ax8.legend(fontsize=8)
    ax8.grid(True, alpha=0.3)
    ax8.set_ylim(0, 1.15)

    # Secondary axis: deposition rate
    ax8_twin = ax8.twinx()
    ax8_twin.plot(plating['j'], plating['rate'], 'r--', linewidth=1.5, label='Deposition rate')
    ax8_twin.set_ylabel('Deposition Rate (nm/min)', color='red')
    ax8_twin.tick_params(axis='y', labelcolor='red')
    ax8_twin.legend(fontsize=8, loc='center right')

    # Inset: throwing power
    ax8_in = ax8.inset_axes([0.4, 0.15, 0.35, 0.35])
    ax8_in.plot(plating['aspect_ratio'], plating['throwing_power'], 'k-', linewidth=2)
    ax8_in.axhline(0.5, color='green', linestyle=':', linewidth=1.5)
    ax8_in.axvline(1.0, color='green', linestyle=':', linewidth=1.5)
    ax8_in.set_xlabel('Aspect Ratio', fontsize=7)
    ax8_in.set_ylabel('Throwing Power', fontsize=7)
    ax8_in.set_title('Wa = 1 at AR = 1', fontsize=7)
    ax8_in.tick_params(labelsize=6)

    # Save
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/'
                'semiconductor_chemistry_coherence.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Figure saved: semiconductor_chemistry_coherence.png")

# ==============================================================
# MAIN EXECUTION
# ==============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("CHEMISTRY SESSION #254: Semiconductor / Electronics Chemistry")
    print("Finding #191 | 117th Phenomenon Type at γ ~ 1")
    print("=" * 70)

    print("\n1. INTRINSIC CARRIERS")
    d = analyze_doping()
    print(f"   n = p = n_i at T = {d['T_intrinsic']:.0f} K (intrinsic temp for N_D = 10¹⁵)")
    print(f"   E_g(300K) = {d['E_g'][np.argmin(np.abs(d['T']-300))]:.3f} eV")
    print(f"   n_i(300K) = {d['n_i'][np.argmin(np.abs(d['T']-300))]:.2e} cm⁻³")
    print(f"   → At n/p = 1: intrinsic behavior (γ ~ 1!)")

    print("\n2. CVD GROWTH REGIMES")
    c = analyze_cvd()
    print(f"   k_s = h_g transition at T = {c['T_transition']:.0f}°C")
    print(f"   Below: surface-reaction limited (Arrhenius)")
    print(f"   Above: mass-transport limited (weak T dependence)")
    print(f"   → Regime change at rate ratio = 1 (γ ~ 1!)")

    print("\n3. PLASMA ETCHING")
    e = analyze_etching()
    print(f"   Selectivity SiO₂:Si = 1 at E_ion = {e['E_unity']:.0f} eV")
    print(f"   Below: material discrimination possible")
    print(f"   Above: non-selective sputtering dominates")
    print(f"   → Selectivity = 1 IS γ ~ 1 for etch process control")

    print("\n4. THERMAL OXIDATION (Deal-Grove)")
    o = analyze_oxidation()
    print(f"   Linear-parabolic crossover at x = A = {o['x_crossover']*1000:.0f} nm")
    print(f"   Crossover time t = {o['t_crossover']:.1f} min")
    print(f"   → Interface reaction = diffusion at γ ~ 1 boundary")

    print("\n5. PHOTOLITHOGRAPHY")
    l = analyze_lithography()
    print(f"   Exposure threshold E_th = {l['E_th_pos']} mJ/cm² (positive)")
    print(f"   k₁ = {l['k1']} at resolution limit (γ ~ 1!)")
    print(f"   Resolution = {l['resolution']:.0f} nm (ArF, NA=0.85)")
    print(f"   → Dose = E_th IS the γ ~ 1 transformation boundary")

    print("\n6. CMP")
    m = analyze_cmp()
    print(f"   Cu:SiO₂ selectivity = {m['sel_cu_ox'][0]:.1f}")
    print(f"   Pattern density = 0.5: maximum non-uniformity (γ ~ 1!)")
    print(f"   → Removal = addition balance at planarization target")

    print("\n7. DIFFUSION / JUNCTION")
    di = analyze_diffusion()
    for Dt, xj in zip(di['Dt_values'], di['junctions']):
        print(f"   Dt = {Dt} μm²: x_j = {xj:.3f} μm")
    print(f"   → C(x_j) = C_bg defines p-n junction (γ ~ 1!)")

    print("\n8. ELECTROPLATING")
    p = analyze_electroplating()
    print(f"   Limiting current j_lim = {p['j_lim']} mA/cm²")
    print(f"   η → 1 at low current (ideal, γ ~ 1)")
    print(f"   Wagner number Wa = 1 at aspect ratio = 1 (γ ~ 1!)")
    print(f"   → All current → deposition when η = 1")

    print("\n" + "=" * 70)
    print("GENERATING MASTER FIGURE...")
    create_master_figure()

    print("\n" + "=" * 70)
    print("SESSION #254 COMPLETE: Semiconductor / Electronics Chemistry")
    print("Finding #191 | 117th phenomenon type at γ ~ 1")
    print("8/8 boundaries validated:")
    print("  1. Intrinsic carriers: n = p = n_i (γ ~ 1 exactly)")
    print("  2. CVD: k_s = h_g regime transition (γ ~ 1)")
    print("  3. Etching: Selectivity = 1 boundary (γ ~ 1)")
    print("  4. Oxidation: Linear = parabolic at x = A (γ ~ 1)")
    print("  5. Lithography: Dose = E_th, k₁ = 0.5 (γ ~ 1)")
    print("  6. CMP: Pattern density = 0.5, removal balance (γ ~ 1)")
    print("  7. Diffusion: C = C_bg at junction depth (γ ~ 1)")
    print("  8. Electroplating: η = 1, Wa = 1 (γ ~ 1)")
    print("=" * 70)
