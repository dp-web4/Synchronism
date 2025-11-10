"""
Atomic-Scale Charge-Phase Simulation with STATIC Charges

Version 2: Fixed charges, only phase and intent evolve.

Key fix from v1:
- Don't evolve charge density ρ (keep charges fixed at source positions)
- Only evolve phase φ and intent I
- Measure equilibrium phase configuration
- Extract V(R) from equilibrium phase structure

This tests: "Given fixed atomic charges, what phase/intent structure emerges?"
Rather than: "How do charges and phase co-evolve?" (which was unstable)

Evolution equations (simplified):
  ρ = fixed at source positions
  ∂φ/∂t = α∇²I - β·ρ        (phase driven by intent AND charge sources)
  ∂I/∂t = -γ(I - I_eq)       (intent relaxation toward charge-dependent equilibrium)

The key is that phase is now directly driven by charge density (source term).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from dataclasses import dataclass
import sys
import os


@dataclass
class StaticChargeConfig:
    """Configuration for static charge simulation"""
    Lx: int = 48  # Larger lattice for better separation range
    Ly: int = 48
    dt: float = 0.01

    # Physical parameters
    alpha_sync: float = 0.1  # Intent → phase coupling
    beta_source: float = 0.5  # Charge → phase coupling (new!)
    gamma_relax: float = 0.3  # Intent relaxation

    # Emergent atomic properties
    charge_e: float = 1.0
    mass_ratio: float = 1836.0

    # Simulation
    n_eq: int = 1000  # Equilibration steps (longer for static case)
    n_meas: int = 100   # Measurements
    meas_interval: int = 50


class StaticChargeLattice:
    """
    Lattice with fixed charges, evolving phase and intent.

    Tests: What phase structure emerges from static charge distribution?
    """

    def __init__(self, config: StaticChargeConfig):
        self.cfg = config

        # State variables
        self.rho = np.zeros((config.Lx, config.Ly))  # FIXED charge density
        self.phi = np.zeros((config.Lx, config.Ly))  # Evolving phase
        self.I = np.ones((config.Lx, config.Ly))     # Evolving intent

        # Track which sites have fixed charges
        self.charge_mask = np.zeros((config.Lx, config.Ly), dtype=bool)

        # Energy history
        self.energy_history = []

    def place_charges(self, charges_and_positions):
        """Place fixed charges that won't move"""
        for q, x, y in charges_and_positions:
            self.rho[x, y] = q * self.cfg.charge_e
            self.charge_mask[x, y] = True

            # Set equilibrium intent based on mass
            if q > 0:  # Proton
                self.I[x, y] = self.cfg.mass_ratio
            elif q < 0:  # Electron
                self.I[x, y] = 1.0

    def compute_laplacian(self, field):
        """5-point stencil with periodic BC"""
        lap = (
            np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0) +
            np.roll(field, 1, axis=1) + np.roll(field, -1, axis=1) -
            4.0 * field
        )
        return lap

    def evolve_step(self):
        """
        Single step with FIXED charges.

        Evolution:
          ρ = constant (fixed at sources)
          ∂φ/∂t = α∇²I - β·ρ     (phase driven by intent laplacian AND charge sources)
          ∂I/∂t = -γ(I - I_eq)   (intent relaxes to equilibrium)
        """
        dt = self.cfg.dt

        # 1. Update phase: ∂φ/∂t = α∇²I - β·ρ
        #    The -β·ρ term is NEW: charges directly source phase field
        lap_I = self.compute_laplacian(self.I)
        dphi_dt = self.cfg.alpha_sync * lap_I - self.cfg.beta_source * self.rho
        self.phi += dt * dphi_dt

        # 2. Update intent: ∂I/∂t = -γ(I - I_eq)
        #    Equilibrium intent depends on whether site has charge
        I_eq = np.ones_like(self.I)
        I_eq[self.charge_mask] = self.I[self.charge_mask]  # Charges keep their mass

        dI_dt = -self.cfg.gamma_relax * (self.I - I_eq)
        dI_dt[self.charge_mask] = 0.0  # Don't change intent at charge sites

        self.I += dt * dI_dt

    def compute_total_energy(self):
        """
        Total energy of the configuration.

        For static charges, energy is:
        E = (1/2)∫[|∇φ|² + (I-1)²]d²x + ∫ρ·φ d²x

        The ρ·φ term is the interaction energy.
        """
        # Phase gradient energy (kinetic)
        grad_phi_x = (np.roll(self.phi, -1, axis=0) - np.roll(self.phi, 1, axis=0)) / 2.0
        grad_phi_y = (np.roll(self.phi, -1, axis=1) - np.roll(self.phi, 1, axis=1)) / 2.0
        E_phase = 0.5 * np.sum(grad_phi_x**2 + grad_phi_y**2)

        # Intent deviation energy (potential)
        E_intent = 0.5 * np.sum((self.I - 1.0)**2)

        # Charge-phase interaction
        E_interaction = np.sum(self.rho * self.phi)

        E_total = E_phase + E_intent + E_interaction

        return E_total

    def measure_phase_difference(self, pos1, pos2):
        """
        Measure phase difference between two positions.

        For static charges, V(R) ~ Δφ = φ(R) - φ(0)
        """
        x1, y1 = pos1
        x2, y2 = pos2

        phi1 = self.phi[x1, y1]
        phi2 = self.phi[x2, y2]

        return phi2 - phi1  # Phase difference

    def measure_interaction_energy(self, pos1, pos2):
        """
        Measure interaction energy between charges.

        For static charges with q1, q2:
        V(R) = q1·q2·Δφ(R)

        where Δφ(R) is the phase difference at separation R.
        """
        x1, y1 = pos1
        x2, y2 = pos2

        q1 = self.rho[x1, y1]
        q2 = self.rho[x2, y2]

        # Phase at each charge position
        phi1 = self.phi[x1, y1]
        phi2 = self.phi[x2, y2]

        # Interaction energy
        V = q1 * q2 * np.abs(phi1 - phi2)

        return V


def run_static_charge_measurement(config: StaticChargeConfig, R_values, verbose=True):
    """
    Measure V(R) for fixed charges at various separations.

    Process:
    1. Place charges at separation R
    2. Evolve phase/intent to equilibrium
    3. Measure equilibrium phase structure
    4. Extract V(R) from phase difference
    """
    V_measurements = {R: [] for R in R_values}

    if verbose:
        print(f"\nMeasuring static potential for {len(R_values)} separations...")
        print(f"Lattice: {config.Lx}×{config.Ly}")
        print(f"Equilibration: {config.n_eq} steps")

    for R in R_values:
        if verbose:
            print(f"\n  R = {R}:")

        # Create lattice
        lattice = StaticChargeLattice(config)

        # Place charges
        cx, cy = config.Lx // 2, config.Ly // 2
        pos1 = (cx - R//2, cy)
        pos2 = (cx + R//2, cy)

        lattice.place_charges([
            (+1, pos1[0], pos1[1]),
            (-1, pos2[0], pos2[1]),
        ])

        if verbose:
            print(f"    Charges: +e at {pos1}, -e at {pos2}")
            print(f"    Equilibrating...", end='', flush=True)

        # Equilibrate
        for step in range(config.n_eq):
            lattice.evolve_step()

            # Track energy every 100 steps
            if step % 100 == 0:
                E = lattice.compute_total_energy()
                lattice.energy_history.append(E)

        if verbose:
            print(" done")

            # Check if equilibrated (energy should be stable)
            if len(lattice.energy_history) > 5:
                E_recent = lattice.energy_history[-5:]
                E_std = np.std(E_recent)
                E_mean = np.mean(E_recent)
                if E_std / abs(E_mean) < 0.01:
                    print(f"    ✓ Equilibrated (δE/E < 1%)")
                else:
                    print(f"    ⚠ May not be equilibrated (δE/E = {E_std/abs(E_mean)*100:.1f}%)")

            print(f"    Measuring ({config.n_meas} samples)...", end='', flush=True)

        # Measure
        for _ in range(config.n_meas):
            # Continue evolution
            for _ in range(config.meas_interval):
                lattice.evolve_step()

            # Measure interaction energy
            V = lattice.measure_interaction_energy(pos1, pos2)
            V_measurements[R].append(V)

        if verbose:
            V_mean = np.mean(V_measurements[R])
            V_std = np.std(V_measurements[R])
            print(f" done")
            print(f"    V = {V_mean:.4f} ± {V_std:.4f}")

    # Compute statistics
    R_array = np.array(R_values)
    V_mean = np.array([np.mean(V_measurements[R]) for R in R_values])
    V_err = np.array([np.std(V_measurements[R]) / np.sqrt(len(V_measurements[R])) for R in R_values])

    return {
        'R': R_array,
        'V': V_mean,
        'dV': V_err,
        'measurements': V_measurements
    }


def fit_coulomb_potential(R, V, dV):
    """Fit to V = -α/R + c"""
    def coulomb(R, alpha, c):
        return -alpha / R + c

    try:
        popt, pcov = curve_fit(coulomb, R, V, sigma=dV, absolute_sigma=True, p0=[1.0, 0.0])
        perr = np.sqrt(np.diag(pcov))

        alpha, c = popt
        dalpha, dc = perr

        V_fit = coulomb(R, alpha, c)
        chi2 = np.sum(((V - V_fit) / dV)**2)
        chi2dof = chi2 / (len(R) - 2)

        return {
            'alpha': alpha,
            'c': c,
            'dalpha': dalpha,
            'dc': dc,
            'chi2dof': chi2dof,
            'V_fit': V_fit
        }
    except Exception as e:
        print(f"Fit failed: {e}")
        return None


def plot_results(results, fit_results=None, filename=None):
    """Plot V(R) with optional fit"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    R = results['R']
    V = results['V']
    dV = results['dV']

    # Left: V vs R (linear scale)
    ax1.errorbar(R, V, yerr=dV, fmt='o', capsize=5, markersize=8,
                 color='blue', label='Measured')

    if fit_results:
        ax1.plot(R, fit_results['V_fit'], 'r-', linewidth=2,
                label=f"Fit: V = -{fit_results['alpha']:.3f}/R + {fit_results['c']:.3f}")

        ax1.text(0.05, 0.95,
                f"α = {fit_results['alpha']:.3f} ± {fit_results['dalpha']:.3f}\n"
                f"c = {fit_results['c']:.3f} ± {fit_results['dc']:.3f}\n"
                f"χ²/dof = {fit_results['chi2dof']:.2f}",
                transform=ax1.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax1.set_xlabel('Separation R (lattice units)')
    ax1.set_ylabel('Potential Energy V(R)')
    ax1.set_title('Static Potential: Fixed Charges, Equilibrium Phase')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Right: V vs 1/R (linearized Coulomb check)
    inv_R = 1.0 / R
    ax2.errorbar(inv_R, V, yerr=dV, fmt='o', capsize=5, markersize=8,
                 color='blue', label='Measured')

    if fit_results:
        V_linear = -fit_results['alpha'] * inv_R + fit_results['c']
        ax2.plot(inv_R, V_linear, 'r-', linewidth=2, label='Linear fit')

    ax2.set_xlabel('1/R (inverse lattice units)')
    ax2.set_ylabel('Potential Energy V(R)')
    ax2.set_title('Linearized: Should be straight line if Coulomb')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to {filename}")
    else:
        plt.show()

    plt.close()


if __name__ == "__main__":
    print("=" * 70)
    print("Atomic-Scale Static Charge Simulation (v2)")
    print("Fixed charges, equilibrium phase structure")
    print("=" * 70)

    config = StaticChargeConfig(
        Lx=48,
        Ly=48,
        dt=0.01,
        alpha_sync=0.1,
        beta_source=0.5,
        gamma_relax=0.3,
        n_eq=1000,
        n_meas=100,
        meas_interval=50
    )

    # Test separations
    R_values = [3, 4, 5, 6, 7, 8, 10, 12, 14, 16]

    print("\nRunning static potential measurement...")
    print("(This will take a few minutes)")

    results = run_static_charge_measurement(config, R_values, verbose=True)

    # Fit
    print("\n" + "=" * 70)
    print("Fitting to Coulomb form: V(R) = -α/R + c")
    print("=" * 70)

    fit_results = fit_coulomb_potential(results['R'], results['V'], results['dV'])

    if fit_results:
        print(f"\nFit Results:")
        print(f"  α = {fit_results['alpha']:.4f} ± {fit_results['dalpha']:.4f}")
        print(f"  c = {fit_results['c']:.4f} ± {fit_results['dc']:.4f}")
        print(f"  χ²/dof = {fit_results['chi2dof']:.3f}")

        significance = abs(fit_results['alpha']) / fit_results['dalpha']
        print(f"\nSignificance: |α|/δα = {significance:.2f}")

        if significance > 2:
            print("  ✓ Coulomb coefficient significantly non-zero (>2σ)")
            print(f"  ✓ Evidence for 1/R potential structure!")
        else:
            print("  ✗ Coulomb coefficient consistent with zero")
            print(f"  ✗ No significant 1/R structure detected")

        if fit_results['chi2dof'] < 2:
            print(f"  ✓ Good fit quality (χ²/dof = {fit_results['chi2dof']:.2f})")
        else:
            print(f"  ⚠ Poor fit quality (χ²/dof = {fit_results['chi2dof']:.2f})")

    # Plot
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Research')
    plot_filename = os.path.join(output_dir, 'Session7_Static_Charge_Potential.png')
    plot_results(results, fit_results, filename=plot_filename)

    print("\n" + "=" * 70)
    print("Simulation complete!")
    print("=" * 70)
