"""
Coherence Chemistry Framework - Visual Summary
Chemistry Session #5: Synthesis

Creates visualization summarizing the unified coherence chemistry framework
across all four domains: superconductivity, catalysis, bonding, phase transitions.

Author: Synchronism Chemistry Track
Date: 2025-01-10
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Arrow
import matplotlib.patches as mpatches

# =============================================================================
# Universal Coherence Function
# =============================================================================

def coherence_function(x, gamma=2.0):
    """
    Universal coherence function: C(x) = tanh(gamma * x)
    """
    return np.tanh(gamma * x)

def phase_barrier(delta_phi, E0=1.0):
    """
    Phase barrier: E = E0 * (1 - cos(delta_phi))
    """
    return E0 * (1 - np.cos(delta_phi))

# =============================================================================
# Main Visualization
# =============================================================================

def create_framework_summary():
    """
    Create a comprehensive summary visualization of the Coherence Chemistry Framework.
    """
    fig = plt.figure(figsize=(16, 14))

    # Use GridSpec for complex layout
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)

    # =========================================================================
    # Top row: Universal coherence function
    # =========================================================================

    ax_coherence = fig.add_subplot(gs[0, :])

    x = np.linspace(-3, 3, 200)

    # Plot coherence function for different gamma values
    for gamma, color, label in [(1.0, 'blue', 'γ=1 (catalysis)'),
                                 (1.5, 'green', 'γ=1.5 (bonding)'),
                                 (2.0, 'red', 'γ=2 (superconductivity)'),
                                 (2.5, 'purple', 'γ=2.5 (liquid crystals)')]:
        y = coherence_function(x, gamma)
        ax_coherence.plot(x, y, color=color, linewidth=2, label=label)

    ax_coherence.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax_coherence.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    ax_coherence.set_xlabel('Coherence Argument x', fontsize=12)
    ax_coherence.set_ylabel('Coherence C(x)', fontsize=12)
    ax_coherence.set_title('Universal Coherence Function: C(x) = tanh(γ × x)', fontsize=14, fontweight='bold')
    ax_coherence.legend(loc='lower right', fontsize=10)
    ax_coherence.grid(True, alpha=0.3)
    ax_coherence.set_xlim(-3, 3)
    ax_coherence.set_ylim(-1.1, 1.1)

    # =========================================================================
    # Middle row: Four domains
    # =========================================================================

    # Domain 1: Superconductivity
    ax1 = fig.add_subplot(gs[1, 0])

    T_ratio = np.linspace(0, 1.2, 100)
    delta_ratio = np.where(T_ratio < 1, np.sqrt(1 - T_ratio**3), 0)

    ax1.plot(T_ratio, delta_ratio, 'b-', linewidth=2)
    ax1.axvline(x=1.0, color='red', linestyle='--', label='T_c')
    ax1.fill_between(T_ratio, 0, delta_ratio, alpha=0.3)

    ax1.set_xlabel('T / T_c', fontsize=11)
    ax1.set_ylabel('Δ(T) / Δ₀', fontsize=11)
    ax1.set_title('Superconductivity\n(Cooper Pair Coherence)', fontsize=12, fontweight='bold')
    ax1.set_xlim(0, 1.2)
    ax1.set_ylim(0, 1.1)
    ax1.text(0.5, 0.5, '2Δ₀/kT_c = 2√π', fontsize=10, ha='center',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    ax1.legend(fontsize=9)

    # Domain 2: Catalysis
    ax2 = fig.add_subplot(gs[1, 1])

    # Reaction coordinate with and without catalyst
    x_rxn = np.linspace(0, 1, 100)
    E_uncat = 4 * x_rxn * (1 - x_rxn)  # Parabolic barrier
    E_cat = 2 * x_rxn * (1 - x_rxn)  # Lower barrier

    ax2.plot(x_rxn, E_uncat, 'r-', linewidth=2, label='Uncatalyzed')
    ax2.plot(x_rxn, E_cat, 'g-', linewidth=2, label='Catalyzed')
    ax2.fill_between(x_rxn, E_cat, E_uncat, alpha=0.3, color='green')

    ax2.annotate('', xy=(0.5, 0.5), xytext=(0.5, 1.0),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax2.text(0.55, 0.75, 'ΔE_a', fontsize=10)

    ax2.set_xlabel('Reaction Coordinate', fontsize=11)
    ax2.set_ylabel('Energy', fontsize=11)
    ax2.set_title('Catalysis\n(Phase Barrier Bridging)', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9, loc='upper right')
    ax2.text(0.5, 0.2, 'E_a ∝ (1-cos Δφ)', fontsize=10, ha='center',
             bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    # Domain 3: Bonding
    ax3 = fig.add_subplot(gs[1, 2])

    delta_phi = np.linspace(0, 2*np.pi, 100)
    E_bond = np.cos(delta_phi)

    ax3.plot(delta_phi * 180/np.pi, E_bond, 'b-', linewidth=2)
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

    # Mark bonding and antibonding
    ax3.scatter([0, 360], [1, 1], s=100, c='green', zorder=5, label='Bonding')
    ax3.scatter([180], [-1], s=100, c='red', zorder=5, label='Antibonding')

    ax3.set_xlabel('Phase Difference Δφ (degrees)', fontsize=11)
    ax3.set_ylabel('Bond Energy (normalized)', fontsize=11)
    ax3.set_title('Chemical Bonding\n(Phase-Locked Orbitals)', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.set_xlim(0, 360)
    ax3.set_xticks([0, 90, 180, 270, 360])

    # =========================================================================
    # Bottom row: Phase transitions and summary
    # =========================================================================

    # Domain 4: Phase transitions
    ax4 = fig.add_subplot(gs[2, 0])

    # Correlation length vs T
    T = np.linspace(0.5, 1.5, 100)
    T_c = 1.0

    xi_crystal = np.where(T < T_c, np.inf * np.ones_like(T), np.nan)
    xi_liquid = np.where(T > T_c, 5 / (T - T_c + 0.1), np.nan)

    # Plot symbolically
    ax4.axvline(x=T_c, color='red', linestyle='--', label='T_m')
    ax4.annotate('Crystal\n(ξ → ∞)', xy=(0.7, 0.7), fontsize=11, ha='center')
    ax4.annotate('Liquid\n(ξ ~ 10 Å)', xy=(1.3, 0.7), fontsize=11, ha='center')

    ax4.set_xlabel('Temperature', fontsize=11)
    ax4.set_ylabel('Correlation Length ξ', fontsize=11)
    ax4.set_title('Phase Transitions\n(Coherence States)', fontsize=12, fontweight='bold')
    ax4.set_xlim(0.5, 1.5)
    ax4.set_ylim(0, 1)
    ax4.set_yticks([])
    ax4.legend(fontsize=9)

    # Summary table
    ax_table = fig.add_subplot(gs[2, 1:])
    ax_table.axis('off')

    # Create summary text
    summary_text = """
    ╔══════════════════════════════════════════════════════════════════════════╗
    ║                    COHERENCE CHEMISTRY FRAMEWORK                         ║
    ╠══════════════════════════════════════════════════════════════════════════╣
    ║  Domain           │ Core Insight                  │ Status               ║
    ╠═══════════════════╪═══════════════════════════════╪══════════════════════╣
    ║  Superconductivity│ BCS = coherence theory        │ DERIVED (2√π)        ║
    ║  Catalysis        │ Phase barrier bridging        │ CONSTRAINED          ║
    ║  Bonding          │ Phase-locked orbitals         │ DERIVED (Hückel)     ║
    ║  Phase Transitions│ Coherence states              │ MIXED (concept OK)   ║
    ╠══════════════════════════════════════════════════════════════════════════╣
    ║  Universal equation: C(x) = tanh(γ × g(x))                               ║
    ║  γ = phase space dimensionality ≈ 2 for most chemistry                   ║
    ║  "Chemistry IS phase physics"                                            ║
    ╚══════════════════════════════════════════════════════════════════════════╝
    """

    ax_table.text(0.5, 0.5, summary_text, fontsize=10, family='monospace',
                  ha='center', va='center',
                  bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/coherence_framework_summary.png',
                dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()

    print("Framework summary visualization saved!")

def create_gamma_comparison():
    """
    Create visualization comparing γ values across domains.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    domains = ['Catalysis\n(1D rxn coord)', 'Bonding\n(dipole)', 'Superconductivity\n(2D Fermi)',
               'Liquid Crystals\n(orientation)', 'Gravity\n(3D space)']
    gamma_values = [1.0, 1.5, 2.0, 2.5, 2.0]
    colors = ['green', 'blue', 'red', 'purple', 'orange']
    statuses = ['CONSTRAINED', 'DERIVED', 'DERIVED', 'CONSTRAINED', 'PRIMARY TRACK']

    bars = ax.bar(domains, gamma_values, color=colors, alpha=0.7, edgecolor='black', linewidth=2)

    # Add value labels
    for bar, val, status in zip(bars, gamma_values, statuses):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'γ = {val}\n({status})',
                ha='center', va='bottom', fontsize=9)

    ax.axhline(y=2.0, color='red', linestyle='--', alpha=0.5, label='γ = 2 (theoretical)')
    ax.set_ylabel('Effective γ', fontsize=12)
    ax.set_title('Phase Space Dimensionality Across Coherence Chemistry Domains', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 3.2)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/gamma_comparison.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    print("Gamma comparison plot saved!")

def create_success_failure_chart():
    """
    Create chart showing successes and failures.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    categories = ['BCS Ratio\n(2√π)', 'Hückel\n(4n+2)', 'Catalyst\nBarrier',
                  'Dipole\ntanh', 'Glass\nFragility', 'Melting\nPoint',
                  'Critical\nExponents', 'Period 3\nAngles']

    # 1 = success, 0 = partial, -1 = failure
    scores = [1, 1, 0.7, 0.8, 0.6, -0.5, -0.5, -0.7]

    colors = ['green' if s > 0.5 else 'yellow' if s > 0 else 'red' for s in scores]

    bars = ax.barh(categories, scores, color=colors, alpha=0.7, edgecolor='black')

    ax.axvline(x=0, color='black', linewidth=2)
    ax.axvline(x=0.5, color='green', linestyle='--', alpha=0.5)
    ax.axvline(x=-0.5, color='red', linestyle='--', alpha=0.5)

    ax.set_xlabel('Success Score', fontsize=12)
    ax.set_title('Coherence Chemistry: Successes and Failures', fontsize=14, fontweight='bold')
    ax.set_xlim(-1.2, 1.2)

    # Add legend
    legend_elements = [mpatches.Patch(facecolor='green', alpha=0.7, label='Success'),
                       mpatches.Patch(facecolor='yellow', alpha=0.7, label='Partial'),
                       mpatches.Patch(facecolor='red', alpha=0.7, label='Failure')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/success_failure_chart.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    print("Success/failure chart saved!")

# =============================================================================
# Main Execution
# =============================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("Chemistry Session #5: Coherence Chemistry Framework Synthesis")
    print("=" * 60)

    print("\nGenerating visualizations...")
    create_framework_summary()
    create_gamma_comparison()
    create_success_failure_chart()

    print("\n=== Summary ===")
    print("Created three visualizations:")
    print("1. coherence_framework_summary.png - Overview of all domains")
    print("2. gamma_comparison.png - γ values across chemistry")
    print("3. success_failure_chart.png - What worked and what didn't")
    print("\nCoherence Chemistry Framework synthesis complete!")
