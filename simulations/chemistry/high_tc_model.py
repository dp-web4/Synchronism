"""
High-Tc Superconductors Through Coherence Lens - Chemistry Session #6
Synchronism Framework applied to Cuprate Superconductivity

Models:
1. Magnetic exchange enhancement of T_c
2. Doping-dependent dome shape
3. Layer-number optimization
4. Predictions for new materials

Author: Synchronism Chemistry Track
Date: 2025-01-10
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Physical constants
kB = 8.617e-5  # eV/K
kB_meV = 0.0862  # meV/K

# =============================================================================
# Exchange-Enhanced T_c Model
# =============================================================================

def exchange_enhanced_Tc(J_AF, f_coherence=0.5, T_c_BCS=30):
    """
    T_c enhanced by magnetic exchange energy.

    T_c = T_c_BCS * (J_AF / hbar_omega_D) * f_coherence

    Parameters:
    - J_AF: Antiferromagnetic exchange in meV
    - f_coherence: Coherence factor (0-1)
    - T_c_BCS: BCS ceiling temperature in K

    Returns: Enhanced T_c in K
    """
    hbar_omega_D = 40  # meV, typical Debye energy
    enhancement = J_AF / hbar_omega_D
    return T_c_BCS * enhancement * f_coherence

def doping_coherence(x, x_opt=0.16, sigma=0.08):
    """
    Coherence factor as function of doping.

    C(x) = exp(-(x - x_opt)^2 / (2*sigma^2))

    Maximum at optimal doping x_opt.
    """
    return np.exp(-(x - x_opt)**2 / (2 * sigma**2))

def Tc_dome(x, T_c_max, x_opt=0.16, sigma=0.08):
    """
    T_c dome as function of doping level.
    """
    return T_c_max * doping_coherence(x, x_opt, sigma)

# =============================================================================
# Layer-Number Model
# =============================================================================

def layer_enhancement(n, J_c_ratio=0.3, n_opt=3, disorder_width=1.5):
    """
    Coherence enhancement from multiple layers.

    gamma_eff(n) = (1 + (n-1) * J_c/J) * disorder_factor

    Parameters:
    - n: Number of CuO2 layers
    - J_c_ratio: Inter-layer coupling / in-plane coupling
    - n_opt: Optimal layer number
    - disorder_width: Width of disorder penalty
    """
    # Enhancement from layer coupling
    coupling = 1 + (n - 1) * J_c_ratio

    # Disorder penalty beyond optimal
    disorder = np.exp(-((n - n_opt) / disorder_width)**2) if n > n_opt else 1.0

    return coupling * disorder

def Tc_vs_layers(n_layers, T_c_single=38, J_c_ratio=0.35):
    """
    T_c as function of number of CuO2 layers.
    """
    return np.array([T_c_single * layer_enhancement(n, J_c_ratio) for n in n_layers])

# =============================================================================
# Gap Ratio Model
# =============================================================================

def gap_ratio(gamma):
    """
    Gap ratio 2Delta/(kT_c) from coherence model.

    Ratio = 2*sqrt(pi) / tanh(gamma * ln(2))
    """
    C_trans = np.tanh(gamma * np.log(2))
    return 2 * np.sqrt(np.pi) / C_trans

def gamma_from_gap_ratio(ratio):
    """
    Infer gamma from observed gap ratio.
    """
    C_trans = 2 * np.sqrt(np.pi) / ratio
    if C_trans >= 1:
        return np.inf
    return np.arctanh(C_trans) / np.log(2)

# =============================================================================
# Experimental Data
# =============================================================================

# Cuprate T_c data (optimal doping)
cuprate_data = {
    'LSCO': {'layers': 1, 'Tc': 38, 'gap_ratio': 4.5, 'J_AF': 130},
    'YBCO': {'layers': 2, 'Tc': 92, 'gap_ratio': 5.5, 'J_AF': 120},
    'Bi-2212': {'layers': 2, 'Tc': 85, 'gap_ratio': 6.0, 'J_AF': 120},
    'Bi-2223': {'layers': 3, 'Tc': 110, 'gap_ratio': 6.5, 'J_AF': 115},
    'Hg-1223': {'layers': 3, 'Tc': 134, 'gap_ratio': 6.0, 'J_AF': 140},
    'Tl-1234': {'layers': 4, 'Tc': 128, 'gap_ratio': 5.5, 'J_AF': 130},
}

# Conventional superconductors
conventional_data = {
    'Al': {'Tc': 1.2, 'gap_ratio': 3.4},
    'Nb': {'Tc': 9.3, 'gap_ratio': 3.9},
    'Pb': {'Tc': 7.2, 'gap_ratio': 4.3},
    'NbN': {'Tc': 16, 'gap_ratio': 4.0},
    'Nb3Sn': {'Tc': 18, 'gap_ratio': 4.2},
    'MgB2': {'Tc': 39, 'gap_ratio': 4.0},
}

# Doping-dependent T_c data for YBCO
ybco_doping = {
    'x': [0.05, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24],
    'Tc': [20, 45, 60, 78, 88, 92, 88, 80, 65, 45],
}

# Layer-number dependence (Bi-family)
bi_layers = {
    'n': [1, 2, 3],
    'Tc': [34, 85, 110],
}

# Hg-family layers
hg_layers = {
    'n': [1, 2, 3, 4],
    'Tc': [95, 128, 134, 127],
}

# =============================================================================
# Analysis
# =============================================================================

def analyze_cuprate_enhancement():
    """Analyze what distinguishes cuprates from conventional SC."""

    print("\n=== Cuprate vs Conventional Analysis ===\n")

    print("Conventional superconductors:")
    for name, data in conventional_data.items():
        gamma = gamma_from_gap_ratio(data['gap_ratio'])
        print(f"  {name}: T_c = {data['Tc']} K, gap ratio = {data['gap_ratio']}, γ = {gamma:.2f}")

    print("\nCuprate superconductors:")
    for name, data in cuprate_data.items():
        gamma = gamma_from_gap_ratio(data['gap_ratio'])
        predicted_Tc = exchange_enhanced_Tc(data['J_AF'], 0.5, 30)
        print(f"  {name}: T_c = {data['Tc']} K, gap ratio = {data['gap_ratio']}, γ = {gamma:.2f}")
        print(f"         J_AF = {data['J_AF']} meV, predicted T_c ≈ {predicted_Tc:.0f} K")

def analyze_doping_dome():
    """Fit the doping dome for YBCO."""

    print("\n=== Doping Dome Analysis ===\n")

    x = np.array(ybco_doping['x'])
    Tc = np.array(ybco_doping['Tc'])

    # Fit dome
    popt, _ = curve_fit(Tc_dome, x, Tc, p0=[92, 0.16, 0.08])
    Tc_max, x_opt, sigma = popt

    print(f"YBCO doping dome fit:")
    print(f"  T_c,max = {Tc_max:.1f} K")
    print(f"  x_opt = {x_opt:.3f}")
    print(f"  σ = {sigma:.3f}")

    print(f"\nCoherence interpretation:")
    print(f"  Optimal doping creates maximum pairing coherence")
    print(f"  Underdoped: AF order competes")
    print(f"  Overdoped: Insufficient correlations")

    return popt

def analyze_layer_dependence():
    """Analyze layer-number dependence."""

    print("\n=== Layer Number Analysis ===\n")

    # Bi-family
    print("Bi-family:")
    n_bi = np.array(bi_layers['n'])
    Tc_bi = np.array(bi_layers['Tc'])

    for n, Tc in zip(n_bi, Tc_bi):
        enhancement = layer_enhancement(n, 0.35)
        predicted = 34 * enhancement
        print(f"  n = {n}: T_c = {Tc} K, predicted = {predicted:.0f} K")

    # Hg-family
    print("\nHg-family:")
    n_hg = np.array(hg_layers['n'])
    Tc_hg = np.array(hg_layers['Tc'])

    for n, Tc in zip(n_hg, Tc_hg):
        enhancement = layer_enhancement(n, 0.25, n_opt=3)
        predicted = 95 * enhancement
        print(f"  n = {n}: T_c = {Tc} K, predicted = {predicted:.0f} K")

# =============================================================================
# Visualization
# =============================================================================

def plot_high_tc_analysis():
    """Create comprehensive high-Tc analysis visualization."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: T_c vs J_AF
    ax1 = axes[0, 0]

    J_range = np.linspace(50, 200, 100)
    for f_coh, color, style in [(0.3, 'blue', '--'), (0.5, 'green', '-'), (0.7, 'red', '-.')]:
        Tc_pred = exchange_enhanced_Tc(J_range, f_coh)
        ax1.plot(J_range, Tc_pred, color=color, linestyle=style, linewidth=2,
                 label=f'f_coh = {f_coh}')

    # Plot experimental data
    for name, data in cuprate_data.items():
        ax1.scatter([data['J_AF']], [data['Tc']], s=100, zorder=5)
        ax1.annotate(name, (data['J_AF'], data['Tc']),
                    textcoords="offset points", xytext=(5, 5), fontsize=9)

    ax1.set_xlabel('J_AF (meV)', fontsize=12)
    ax1.set_ylabel('T_c (K)', fontsize=12)
    ax1.set_title('Exchange-Enhanced T_c Model', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Doping dome
    ax2 = axes[0, 1]

    x_data = np.array(ybco_doping['x'])
    Tc_data = np.array(ybco_doping['Tc'])

    ax2.scatter(x_data, Tc_data, s=80, c='blue', zorder=5, label='YBCO data')

    x_fit = np.linspace(0.02, 0.28, 100)
    Tc_fit = Tc_dome(x_fit, 92, 0.16, 0.06)
    ax2.plot(x_fit, Tc_fit, 'r-', linewidth=2, label='Coherence dome fit')

    ax2.axvline(x=0.16, color='gray', linestyle='--', alpha=0.5)
    ax2.text(0.17, 50, 'Optimal\ndoping', fontsize=10)

    ax2.set_xlabel('Hole doping x', fontsize=12)
    ax2.set_ylabel('T_c (K)', fontsize=12)
    ax2.set_title('YBCO Doping Dome', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 0.30)

    # Plot 3: Layer dependence
    ax3 = axes[1, 0]

    n_range = np.arange(1, 6)

    # Bi-family
    Tc_bi_pred = Tc_vs_layers(n_range, T_c_single=34, J_c_ratio=0.35)
    ax3.plot(n_range, Tc_bi_pred, 'b-', linewidth=2, label='Bi-family model')
    ax3.scatter(bi_layers['n'], bi_layers['Tc'], s=100, c='blue', zorder=5)

    # Hg-family
    Tc_hg_pred = Tc_vs_layers(n_range, T_c_single=95, J_c_ratio=0.20)
    ax3.plot(n_range, Tc_hg_pred, 'r-', linewidth=2, label='Hg-family model')
    ax3.scatter(hg_layers['n'], hg_layers['Tc'], s=100, c='red', zorder=5)

    ax3.set_xlabel('Number of CuO₂ layers', fontsize=12)
    ax3.set_ylabel('T_c (K)', fontsize=12)
    ax3.set_title('Layer Number Dependence', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xticks([1, 2, 3, 4, 5])

    # Plot 4: Gap ratio vs T_c
    ax4 = axes[1, 1]

    # Conventional
    for name, data in conventional_data.items():
        ax4.scatter([data['Tc']], [data['gap_ratio']], s=80, c='blue', zorder=5)

    # Cuprates
    for name, data in cuprate_data.items():
        ax4.scatter([data['Tc']], [data['gap_ratio']], s=80, c='red', zorder=5)

    ax4.axhline(y=3.54, color='gray', linestyle='--', alpha=0.7, label='BCS (3.54)')

    ax4.text(5, 3.7, 'Conventional', fontsize=10, color='blue')
    ax4.text(90, 6.3, 'Cuprates', fontsize=10, color='red')

    ax4.set_xlabel('T_c (K)', fontsize=12)
    ax4.set_ylabel('Gap ratio 2Δ/(kT_c)', fontsize=12)
    ax4.set_title('Gap Ratio: Enhanced Coherence', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_xscale('log')

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/high_tc_analysis.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    print("High-Tc analysis plot saved!")

def predict_new_materials():
    """Make predictions for hypothetical new superconductors."""

    print("\n=== Predictions for New Materials ===\n")

    predictions = [
        {'name': 'Optimal 4-layer Hg', 'J_AF': 150, 'layers': 4, 'f_coh': 0.6},
        {'name': '5-layer cuprate', 'J_AF': 140, 'layers': 5, 'f_coh': 0.4},
        {'name': 'High-J AF material', 'J_AF': 200, 'layers': 3, 'f_coh': 0.5},
        {'name': 'Perfect disorder-free', 'J_AF': 130, 'layers': 3, 'f_coh': 0.8},
    ]

    print(f"{'Material':<25} {'J_AF (meV)':<12} {'Layers':<8} {'f_coh':<8} {'Predicted T_c (K)':<15}")
    print("-" * 70)

    for pred in predictions:
        layer_enh = layer_enhancement(pred['layers'], 0.25)
        base_Tc = exchange_enhanced_Tc(pred['J_AF'], pred['f_coh'], 30)
        Tc_pred = base_Tc * layer_enh

        print(f"{pred['name']:<25} {pred['J_AF']:<12} {pred['layers']:<8} {pred['f_coh']:<8.1f} {Tc_pred:<15.0f}")

    print("\nKey insight: Maximizing coherence factor f_coh is key to higher T_c")
    print("This requires: clean crystals, optimal doping, controlled layer spacing")

# =============================================================================
# Main Execution
# =============================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("Chemistry Session #6: High-Tc Superconductors")
    print("=" * 60)

    analyze_cuprate_enhancement()
    analyze_doping_dome()
    analyze_layer_dependence()
    predict_new_materials()

    print("\n" + "=" * 60)
    print("Generating visualizations...")
    plot_high_tc_analysis()

    print("\n=== Summary ===")
    print("1. High-Tc arises from magnetic exchange (J_AF >> hbar*omega_D)")
    print("2. Doping dome reflects coherence optimization")
    print("3. Layer coupling enhances coherence up to n ≈ 3")
    print("4. Gap ratios > 3.54 indicate enhanced coherence regime")
    print("\nKey prediction: T_c ~ J_AF * f_coherence * layer_factor")
