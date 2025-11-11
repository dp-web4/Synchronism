"""
Synchronism Quantum Mechanics Test - Session #10

Testing if (I, φ) dynamics from Synchronism reproduces Schrödinger equation.

Setup:
  - 1D Gaussian wave packet: ψ(x,0) = √I(x,0) · e^{iφ(x,0)}
  - Evolve via Synchronism equations for (I, φ)
  - Compare to analytical Schrödinger evolution

Synchronism dynamics (Madelung-like form):
  ∂I/∂t = -∂/∂x(I · ∂φ/∂x)  [continuity equation]
  ∂φ/∂t = -(1/2m)(∂φ/∂x)² + Q[I]  [Hamilton-Jacobi with quantum potential]

Where Q[I] is quantum potential from intent gradients:
  Q[I] = (ℏ²/2m) · ∂²√I/∂x² / √I

This is the KEY TEST: Does Synchronism give Schrödinger?
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import os


class SynchronismQuantum1D:
    """
    1D quantum evolution via Synchronism (I, φ) dynamics.

    Tests if Synchronism reproduces Schrödinger equation.
    """

    def __init__(self, L, N, m=1.0, hbar=1.0):
        """
        Initialize 1D quantum system.

        Args:
            L: System length
            N: Number of lattice points
            m: Particle mass (natural units)
            hbar: Reduced Planck constant
        """
        self.L = L
        self.N = N
        self.dx = L / N
        self.m = m
        self.hbar = hbar

        # Spatial grid
        self.x = np.linspace(-L/2, L/2, N, endpoint=False)

        # Fields
        self.I = np.zeros(N)  # Intent density (= |ψ|²)
        self.phi = np.zeros(N)  # Phase field (= arg(ψ))

    def initialize_gaussian_wavepacket(self, x0, k0, sigma):
        """
        Initialize Gaussian wave packet.

        ψ(x,0) = (2πσ²)^{-1/4} · exp(-(x-x0)²/(4σ²) + ik0·x)

        Args:
            x0: Initial position
            k0: Initial momentum (p = ℏk0)
            sigma: Wave packet width
        """
        # Intent density I = |ψ|²
        self.I = (2*np.pi*sigma**2)**(-0.5) * np.exp(-(self.x - x0)**2 / (2*sigma**2))

        # Phase field φ
        self.phi = k0 * self.x

        # Normalize
        norm = np.sum(self.I) * self.dx
        self.I /= norm

        print(f"  Initialized Gaussian wave packet:")
        print(f"    x0 = {x0}, k0 = {k0}, σ = {sigma}")
        print(f"    ∫I dx = {np.sum(self.I)*self.dx:.6f} (should be 1.0)")

    def get_psi(self):
        """Get wave function ψ = √I · e^{iφ}"""
        return np.sqrt(self.I) * np.exp(1j * self.phi)

    def compute_quantum_potential(self):
        """
        Compute quantum potential from intent density.

        Q = (ℏ²/2m) · ∂²√I/∂x² / √I
          = (ℏ²/2m) · [∂²I/(4I^{3/2}) - (∂I)²/(4I^{5/2})]

        This is the key Synchronism prediction!
        """
        # Compute ∂I/∂x (central difference)
        dI_dx = np.gradient(self.I, self.dx)

        # Compute ∂²I/∂x² (second derivative)
        d2I_dx2 = np.gradient(dI_dx, self.dx)

        # Quantum potential (avoid division by zero)
        sqrt_I = np.sqrt(self.I + 1e-12)
        Q = (self.hbar**2 / (2*self.m)) * (
            d2I_dx2 / (2*sqrt_I**3) - dI_dx**2 / (4*sqrt_I**5)
        )

        return Q

    def evolve_synchronism(self, dt, steps):
        """
        Evolve (I, φ) via Synchronism dynamics.

        ∂I/∂t = -∂/∂x(I · ∂φ/∂x)
        ∂φ/∂t = -(1/2m)(∂φ/∂x)² + Q[I]

        Using finite difference with periodic BC.
        """
        for step in range(steps):
            # Compute spatial derivatives
            dphi_dx = np.gradient(self.phi, self.dx)
            Q = self.compute_quantum_potential()

            # dI/dt = -∂/∂x(I · ∂φ/∂x)  [continuity]
            flux = self.I * dphi_dx
            dflux_dx = np.gradient(flux, self.dx)
            dI_dt = -dflux_dx

            # dφ/dt = -(1/2m)(∂φ/∂x)² + Q  [Hamilton-Jacobi + quantum potential]
            dphi_dt = -(1/(2*self.m)) * dphi_dx**2 + Q

            # Euler step (simple but stable for small dt)
            self.I += dt * dI_dt
            self.phi += dt * dphi_dt

            # Prevent negative density
            self.I = np.maximum(self.I, 0.0)

            # Renormalize (conserve probability)
            norm = np.sum(self.I) * self.dx
            if norm > 0:
                self.I /= norm


def analytical_schrodinger_1d(x, t, x0, k0, sigma, m=1.0, hbar=1.0):
    """
    Analytical solution for free particle Gaussian wave packet.

    ψ(x,t) = (2πσ²(t))^{-1/4} · exp(-(x-x0-v0·t)²/(4σ²(t)) + i·k0·(x-v0·t/2))

    Where:
      σ²(t) = σ²·(1 + (ℏt/(mσ²))²)  [spreading]
      v0 = ℏk0/m  [group velocity]
    """
    v0 = hbar * k0 / m
    sigma_t_sq = sigma**2 * (1 + (hbar*t / (m*sigma**2))**2)

    # Gaussian envelope (spreading)
    envelope = (2*np.pi*sigma_t_sq)**(-0.25) * np.exp(-(x - x0 - v0*t)**2 / (4*sigma_t_sq))

    # Phase (moving wave with time-dependent correction)
    phase = k0 * (x - v0*t/2) - (hbar*k0**2*t)/(2*m)

    psi = envelope * np.exp(1j * phase)

    return psi


def run_quantum_test(x0=0.0, k0=5.0, sigma=1.0, T=2.0, dt=0.001, L=20.0, N=256):
    """
    Run Synchronism vs Schrödinger comparison test.

    Args:
        x0: Initial position
        k0: Initial momentum
        sigma: Wave packet width
        T: Total evolution time
        dt: Time step
        L: System length
        N: Spatial grid points

    Returns:
        dict with results
    """
    print(f"\n{'='*70}")
    print(f"Synchronism Quantum Mechanics Test - Session #10")
    print(f"{'='*70}")
    print(f"\nParameters:")
    print(f"  Initial position: x0 = {x0}")
    print(f"  Initial momentum: k0 = {k0} (p = ℏk0)")
    print(f"  Wave packet width: σ = {sigma}")
    print(f"  Evolution time: T = {T}")
    print(f"  Time step: dt = {dt}")
    print(f"  Grid: {N} points over L = {L}")

    # Initialize Synchronism
    sync = SynchronismQuantum1D(L, N, m=1.0, hbar=1.0)
    sync.initialize_gaussian_wavepacket(x0, k0, sigma)

    # Store initial state
    psi_sync_initial = sync.get_psi()

    # Evolve via Synchronism
    print(f"\n  Evolving via Synchronism dynamics...")
    steps = int(T / dt)
    sync.evolve_synchronism(dt, steps)

    # Get final Synchronism state
    psi_sync_final = sync.get_psi()
    I_sync = sync.I
    phi_sync = sync.phi

    print(f"    Evolution complete ({steps} steps)")
    print(f"    Final norm: ∫I dx = {np.sum(I_sync)*sync.dx:.6f}")

    # Analytical Schrödinger solution
    print(f"\n  Computing analytical Schrödinger solution...")
    psi_qm_initial = analytical_schrodinger_1d(sync.x, 0.0, x0, k0, sigma)
    psi_qm_final = analytical_schrodinger_1d(sync.x, T, x0, k0, sigma)

    # Compute errors
    error_initial = np.sqrt(np.sum(np.abs(psi_sync_initial - psi_qm_initial)**2) * sync.dx)
    error_final = np.sqrt(np.sum(np.abs(psi_sync_final - psi_qm_final)**2) * sync.dx)

    # Density errors
    I_qm_final = np.abs(psi_qm_final)**2
    density_error = np.sqrt(np.sum((I_sync - I_qm_final)**2) * sync.dx)
    density_error_rel = density_error / np.sqrt(np.sum(I_qm_final**2) * sync.dx)

    print(f"\n{'='*70}")
    print(f"RESULTS")
    print(f"{'='*70}")
    print(f"  ψ error (initial): {error_initial:.6f}")
    print(f"  ψ error (final):   {error_final:.6f}")
    print(f"  Density error (absolute): {density_error:.6f}")
    print(f"  Density error (relative): {density_error_rel*100:.2f}%")

    # Assessment
    if error_final < 0.1:
        print(f"\n  ✅ SUCCESS: Synchronism reproduces Schrödinger (error < 10%)")
    elif error_final < 0.5:
        print(f"\n  ⚠️  PARTIAL: Synchronism approximates Schrödinger (error < 50%)")
    else:
        print(f"\n  ❌ FAILURE: Synchronism differs from Schrödinger (error > 50%)")

    return {
        'x': sync.x,
        't': T,
        'psi_sync': psi_sync_final,
        'psi_qm': psi_qm_final,
        'I_sync': I_sync,
        'I_qm': I_qm_final,
        'phi_sync': phi_sync,
        'error': error_final,
        'density_error_rel': density_error_rel
    }


def plot_comparison(results, filename=None):
    """Plot Synchronism vs Schrödinger comparison"""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    x = results['x']
    t = results['t']

    # Density |ψ|²
    ax = axes[0]
    ax.plot(x, results['I_sync'], 'b-', linewidth=2, label='Synchronism I(x,t)')
    ax.plot(x, results['I_qm'], 'r--', linewidth=2, label='Schrödinger |ψ(x,t)|²')
    ax.set_xlabel('Position x', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)
    ax.set_title(f'Quantum Evolution Comparison (t = {t:.2f})', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Phase
    ax = axes[1]
    phi_qm = np.angle(results['psi_qm'])
    ax.plot(x, results['phi_sync'], 'b-', linewidth=2, label='Synchronism φ(x,t)')
    ax.plot(x, phi_qm, 'r--', linewidth=2, label='Schrödinger arg(ψ)')
    ax.set_xlabel('Position x', fontsize=12)
    ax.set_ylabel('Phase (radians)', fontsize=12)
    ax.set_title('Phase Field Comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Real part of ψ
    ax = axes[2]
    psi_sync_real = np.real(results['psi_sync'])
    psi_qm_real = np.real(results['psi_qm'])
    ax.plot(x, psi_sync_real, 'b-', linewidth=2, label='Synchronism Re(ψ)')
    ax.plot(x, psi_qm_real, 'r--', linewidth=2, label='Schrödinger Re(ψ)')
    ax.set_xlabel('Position x', fontsize=12)
    ax.set_ylabel('Re(ψ)', fontsize=12)
    ax.set_title('Wave Function Real Part', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add error info
    error_text = f"L² error: {results['error']:.4f}\nDensity error: {results['density_error_rel']*100:.2f}%"
    fig.text(0.15, 0.02, error_text, fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {filename}")
    else:
        plt.show()

    plt.close()


if __name__ == "__main__":
    print("\n" + "="*70)
    print("SYNCHRONISM QUANTUM MECHANICS TEST - Session #10")
    print("="*70)
    print("\nTesting if Synchronism (I, φ) dynamics reproduces Schrödinger equation")
    print("\nMethod:")
    print("  1. Initialize Gaussian wave packet: ψ = √I · e^{iφ}")
    print("  2. Evolve via Synchronism dynamics (Madelung form)")
    print("  3. Compare to analytical Schrödinger solution")
    print("\nKey test: Does quantum potential Q[I] = (ℏ²/2m)·∂²√I/∂x²/√I work?")

    # Run test
    results = run_quantum_test(
        x0=0.0,      # Center at origin
        k0=5.0,      # Moderate momentum
        sigma=1.0,   # Unit width
        T=2.0,       # Evolve for 2 time units
        dt=0.001,    # Small time step for stability
        L=40.0,      # Large box to avoid boundary effects
        N=512        # High resolution
    )

    # Plot
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Research')
    plot_filename = os.path.join(output_dir, 'Session10_Quantum_Emergence.png')
    plot_comparison(results, filename=plot_filename)

    print("\n" + "="*70)
    print("SESSION #10 QUANTUM TEST COMPLETE")
    print("="*70)
    print("\nKey findings:")
    print("  1. Synchronism (I, φ) dynamics implemented")
    print("  2. Quantum potential Q[I] computed from intent gradients")
    print("  3. Compared to analytical Schrödinger solution")
    print("  4. Tests if quantum mechanics emerges from Synchronism")
    print("\nComparison to Sessions #8-9:")
    print("  - Session #8: Coulomb (classical electrostatics) ✓")
    print("  - Session #9: Magnetism (classical magnetostatics) ✓")
    print("  - Session #10: Schrödinger (quantum mechanics) ?")
    print("\n" + "="*70)
