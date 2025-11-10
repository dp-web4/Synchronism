"""
Atomic-Scale Charge-Phase Coupled Lattice Simulation

Tests Coulomb potential emergence from charge-phase dynamics at correct abstraction level.

Key differences from Session #6:
- Uses atomic-scale degrees of freedom (charge, phase, intent density)
- Encodes emergent properties from Planck→Atomic (q=±e, m, α_EM)
- Tests Atomic→Molecular emergence (Coulomb from charge correlations)
- Correct scale abstraction per Scale_and_Abstraction.md framework

Physical setup:
- Lattice represents atomic-scale regions (~0.5 Å spacing)
- Each site has inherited emergent properties:
  * Charge density ρ (±e for electron/proton, 0 for vacuum)
  * Coherence phase φ (tracks collective sub-MRH phase)
  * Intent density I (related to mass via m ~ ∫I)

Evolution equations:
  ∂ρ/∂t + ∇·j = 0           (charge conservation)
  ∂φ/∂t = α∇²I              (phase tracking from Synchronism)
  j = -ρ∇φ                   (phase-charge coupling)
  ∂I/∂t = -γ(I - I_eq)       (intent relaxation toward equilibrium)

Measurement:
- Static potential V(R) between fixed charges
- Test: Does V(R) ∝ 1/R emerge naturally?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from dataclasses import dataclass
import sys
import os

# Add parent directory tools for statistics
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'private-context', 'tools', 'lattice-gauge'))
try:
    from stats_utils import jackknife_analysis, calculate_blocking_error
except ImportError:
    print("Warning: stats_utils not available, using simple errors")
    def jackknife_analysis(data, func):
        mean = func(data)
        # Simple standard deviation estimate
        n = len(data)
        samples = [func(np.delete(data, i, axis=0)) for i in range(min(n, 20))]
        error = np.std(samples)
        return mean, error

    def calculate_blocking_error(data):
        return np.std(data) / np.sqrt(len(data))


@dataclass
class AtomicLatticeConfig:
    """Configuration for atomic-scale lattice simulation"""
    Lx: int = 32  # Spatial extent (x)
    Ly: int = 32  # Spatial extent (y)
    dt: float = 0.01  # Time step

    # Physical constants (in atomic units where ℏ = m_e = e = 1)
    alpha_sync: float = 0.1  # Phase tracking strength (∂φ/∂t = α∇²I)
    gamma_relax: float = 0.5  # Intent relaxation rate
    coupling: float = 1.0  # Charge-phase coupling (j = -g·ρ∇φ)

    # Emergent properties (inherited from Planck→Atomic)
    charge_e: float = 1.0  # Elementary charge (in atomic units)
    mass_ratio: float = 1836.0  # m_proton / m_electron
    fine_structure: float = 1.0/137.0  # α_EM (not directly used but available)

    # Simulation parameters
    n_therm: int = 200  # Thermalization steps
    n_meas: int = 500  # Measurement steps
    meas_interval: int = 10  # Steps between measurements

    # Initial conditions
    temperature: float = 0.0  # T=0 for static potential (no thermal fluctuations)


class ChargePhaseAtomicLattice:
    """
    Atomic-scale lattice with charge-phase-intent coupled dynamics.

    This is the CORRECT abstraction for testing Coulomb emergence at atomic scale.
    Session #6 used bare Planck-scale variables - this uses atomic-scale effective DOF.
    """

    def __init__(self, config: AtomicLatticeConfig):
        self.cfg = config

        # State variables at each lattice site
        self.rho = np.zeros((config.Lx, config.Ly))  # Charge density
        self.phi = np.zeros((config.Lx, config.Ly))  # Coherence phase
        self.I = np.ones((config.Lx, config.Ly))     # Intent density (initialized to equilibrium)

        # Current density (derived from ρ and φ)
        self.jx = np.zeros((config.Lx, config.Ly))
        self.jy = np.zeros((config.Lx, config.Ly))

        # Energy tracking
        self.total_energy = []
        self.charge_energy = []
        self.phase_energy = []

    def place_charges(self, charges_and_positions):
        """
        Place fixed charges on lattice.

        Args:
            charges_and_positions: List of (charge, x, y) tuples

        Example:
            place_charges([(+1, 10, 16), (-1, 22, 16)])  # electron-proton pair
        """
        for q, x, y in charges_and_positions:
            self.rho[x, y] = q * self.cfg.charge_e
            # Fixed charges have intent density proportional to mass
            if q > 0:  # Proton
                self.I[x, y] = self.cfg.mass_ratio
            elif q < 0:  # Electron
                self.I[x, y] = 1.0

    def compute_gradients(self, field):
        """Compute ∇field using central differences with periodic BC"""
        grad_x = np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)
        grad_y = np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)
        grad_x /= 2.0
        grad_y /= 2.0
        return grad_x, grad_y

    def compute_laplacian(self, field):
        """Compute ∇²field using 5-point stencil with periodic BC"""
        lap = (
            np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0) +
            np.roll(field, 1, axis=1) + np.roll(field, -1, axis=1) -
            4.0 * field
        )
        return lap

    def evolve_step(self):
        """
        Single time step of charge-phase-intent coupled dynamics.

        Evolution equations:
          ∂ρ/∂t = -∇·j              (charge conservation)
          ∂φ/∂t = α∇²I              (phase tracking from Synchronism)
          ∂I/∂t = -γ(I - I_eq)      (intent relaxation)
          j = -g·ρ∇φ                (charge-phase coupling)
        """
        dt = self.cfg.dt

        # 1. Update current from charge-phase coupling: j = -g·ρ∇φ
        grad_phi_x, grad_phi_y = self.compute_gradients(self.phi)
        self.jx = -self.cfg.coupling * self.rho * grad_phi_x
        self.jy = -self.cfg.coupling * self.rho * grad_phi_y

        # 2. Update charge from continuity: ∂ρ/∂t = -∇·j
        div_j = (
            (np.roll(self.jx, -1, axis=0) - np.roll(self.jx, 1, axis=0)) +
            (np.roll(self.jy, -1, axis=1) - np.roll(self.jy, 1, axis=1))
        ) / 2.0

        # Don't update charge at fixed source positions (keep them fixed)
        # Identify sources: sites with large intent density
        fixed_mask = np.abs(self.I - 1.0) > 0.1
        drho_dt = -div_j
        drho_dt[fixed_mask] = 0.0  # Keep sources fixed

        self.rho += dt * drho_dt

        # 3. Update phase from Synchronism: ∂φ/∂t = α∇²I
        lap_I = self.compute_laplacian(self.I)
        self.phi += dt * self.cfg.alpha_sync * lap_I

        # 4. Update intent from relaxation: ∂I/∂t = -γ(I - I_eq)
        # Equilibrium intent is 1.0 for vacuum, higher for massive particles
        I_eq = np.ones_like(self.I)
        I_eq[fixed_mask] = self.I[fixed_mask]  # Sources maintain their intent

        dI_dt = -self.cfg.gamma_relax * (self.I - I_eq)
        dI_dt[fixed_mask] = 0.0  # Keep source intent fixed

        self.I += dt * dI_dt

    def compute_energy(self):
        """
        Compute system energy components.

        Returns tuple: (total, charge_energy, phase_energy)
        """
        # Charge energy: E_q = (1/2) ∫ ρ V d²x
        # where V is electrostatic potential (related to phase)
        grad_phi_x, grad_phi_y = self.compute_gradients(self.phi)
        E_charge = 0.5 * np.sum(self.rho * self.phi)

        # Phase energy: E_φ = (1/2) ∫ |∇φ|² d²x (kinetic energy of phase field)
        E_phase = 0.5 * np.sum(grad_phi_x**2 + grad_phi_y**2)

        # Intent energy: E_I = (1/2) ∫ (I - I_eq)² d²x (potential energy of intent field)
        I_eq = np.ones_like(self.I)
        E_intent = 0.5 * np.sum((self.I - I_eq)**2)

        E_total = E_charge + E_phase + E_intent

        return E_total, E_charge, E_phase

    def measure_potential_energy(self, pos1, pos2):
        """
        Measure interaction energy between two charge positions.

        This is done by computing energy difference:
        V(R) = E[both charges] - E[charge 1 alone] - E[charge 2 alone]

        Args:
            pos1, pos2: (x, y) positions of the two charges

        Returns:
            Interaction energy V(R)
        """
        # For now, use simpler proxy: phase difference weighted by charge
        x1, y1 = pos1
        x2, y2 = pos2

        # Measure phase at charge positions
        phi1 = self.phi[x1, y1]
        phi2 = self.phi[x2, y2]

        # Interaction energy ~ q1·q2·(φ1 - φ2) / R
        # This is approximation - more sophisticated: integrate energy density
        V = self.rho[x1, y1] * self.rho[x2, y2] * (phi1 - phi2)

        return V


def run_static_potential_measurement(config: AtomicLatticeConfig, R_values, verbose=True):
    """
    Measure static potential V(R) for various separations.

    Args:
        config: Lattice configuration
        R_values: List of separation distances to measure
        verbose: Print progress

    Returns:
        dict with 'R', 'V', 'dV' (mean potential and error vs separation)
    """
    n_R = len(R_values)
    V_measurements = {R: [] for R in R_values}

    if verbose:
        print(f"Measuring static potential for {n_R} separation values...")
        print(f"Lattice: {config.Lx}×{config.Ly}, T={config.temperature}")
        print(f"Thermalization: {config.n_therm}, Measurements: {config.n_meas}")

    for R in R_values:
        if verbose:
            print(f"\n  R = {R}:")

        # Create lattice
        lattice = ChargePhaseAtomicLattice(config)

        # Place charges at separation R (centered on lattice)
        cx, cy = config.Lx // 2, config.Ly // 2
        pos1 = (cx - R//2, cy)
        pos2 = (cx + R//2, cy)

        lattice.place_charges([
            (+1, pos1[0], pos1[1]),  # Proton
            (-1, pos2[0], pos2[1]),  # Electron
        ])

        if verbose:
            print(f"    Placed +e at {pos1}, -e at {pos2}")

        # Thermalization
        if verbose:
            print(f"    Thermalizing ({config.n_therm} steps)...", end='', flush=True)

        for _ in range(config.n_therm):
            lattice.evolve_step()

        if verbose:
            print(" done")
            print(f"    Measuring ({config.n_meas} samples)...", end='', flush=True)

        # Measurements
        for _ in range(config.n_meas):
            # Evolve between measurements
            for _ in range(config.meas_interval):
                lattice.evolve_step()

            # Measure potential
            V = lattice.measure_potential_energy(pos1, pos2)
            V_measurements[R].append(V)

        if verbose:
            V_mean = np.mean(V_measurements[R])
            V_std = np.std(V_measurements[R])
            print(f" done (V = {V_mean:.4f} ± {V_std:.4f})")

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
    """
    Fit V(R) to Coulomb form: V = -α/R + c

    Args:
        R: Separation distances
        V: Measured potentials
        dV: Errors on V

    Returns:
        dict with 'alpha', 'c', 'dalpha', 'dc', 'chi2dof'
    """
    # Coulomb fit function
    def coulomb(R, alpha, c):
        return -alpha / R + c

    try:
        # Weighted fit
        popt, pcov = curve_fit(coulomb, R, V, sigma=dV, absolute_sigma=True, p0=[1.0, 0.0])
        perr = np.sqrt(np.diag(pcov))

        alpha, c = popt
        dalpha, dc = perr

        # Compute chi-squared
        V_fit = coulomb(R, alpha, c)
        chi2 = np.sum(((V - V_fit) / dV)**2)
        chi2dof = chi2 / (len(R) - 2)  # 2 parameters

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


def plot_static_potential(results, fit_results=None, filename=None):
    """Plot V(R) with fit if available"""
    R = results['R']
    V = results['V']
    dV = results['dV']

    plt.figure(figsize=(10, 6))

    # Data
    plt.errorbar(R, V, yerr=dV, fmt='o', label='Measured V(R)',
                 capsize=5, markersize=8, color='blue')

    # Fit
    if fit_results is not None:
        plt.plot(R, fit_results['V_fit'], 'r-', linewidth=2,
                label=f"Fit: V = -{fit_results['alpha']:.3f}/R + {fit_results['c']:.3f}")

        # Add fit quality
        plt.text(0.05, 0.95,
                f"α = {fit_results['alpha']:.3f} ± {fit_results['dalpha']:.3f}\n"
                f"c = {fit_results['c']:.3f} ± {fit_results['dc']:.3f}\n"
                f"χ²/dof = {fit_results['chi2dof']:.2f}",
                transform=plt.gca().transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                fontsize=11)

    plt.xlabel('Separation R (lattice units)', fontsize=12)
    plt.ylabel('Potential Energy V(R)', fontsize=12)
    plt.title('Static Potential from Charge-Phase Dynamics\n(Atomic-Scale Abstraction)', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {filename}")
    else:
        plt.show()

    plt.close()


if __name__ == "__main__":
    print("=" * 70)
    print("Atomic-Scale Charge-Phase Coupled Simulation")
    print("Testing Coulomb Emergence at Correct Abstraction Level")
    print("=" * 70)

    # Configuration
    config = AtomicLatticeConfig(
        Lx=32,
        Ly=32,
        dt=0.01,
        alpha_sync=0.1,
        gamma_relax=0.5,
        coupling=1.0,
        n_therm=200,
        n_meas=500,
        meas_interval=10,
        temperature=0.0
    )

    # Separation distances to probe (in lattice units ~ 0.5 Å)
    R_values = [2, 3, 4, 5, 6, 7, 8, 10, 12]

    print("\nRunning static potential measurement...")
    print("(This will take several minutes)")

    # Run measurement
    results = run_static_potential_measurement(config, R_values, verbose=True)

    # Fit to Coulomb
    print("\nFitting to Coulomb form: V(R) = -α/R + c")
    fit_results = fit_coulomb_potential(results['R'], results['V'], results['dV'])

    if fit_results:
        print(f"\nFit Results:")
        print(f"  α = {fit_results['alpha']:.4f} ± {fit_results['dalpha']:.4f}")
        print(f"  c = {fit_results['c']:.4f} ± {fit_results['dc']:.4f}")
        print(f"  χ²/dof = {fit_results['chi2dof']:.3f}")

        # Statistical significance of α
        significance = abs(fit_results['alpha']) / fit_results['dalpha']
        print(f"\nSignificance: |α|/δα = {significance:.2f}")
        if significance > 2:
            print("  ✓ α significantly different from zero (>2σ)")
        else:
            print("  ✗ α consistent with zero (null result)")

    # Plot
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Research')
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, 'Session7_Atomic_Scale_Potential.png')

    plot_static_potential(results, fit_results, filename=plot_filename)

    print("\n" + "=" * 70)
    print("Simulation complete!")
    print("=" * 70)
