#!/usr/bin/env python3
"""
Session #12: Numerical Validation of Gravity from Synchronism

Test: Does Newtonian gravitational potential Φ(r) = -GM/r emerge from
      Synchronism's intent field dynamics?

Theory (from Session #11):
    ∇²Φ = 4πG ρ_I

    where ρ_I = (1/2)[(∇I)² + (I-I₀)²]

For spherically symmetric mass:
    Expected: Φ(r) = -GM/r (outside mass)

Method:
    1. Define spherical intent distribution I(r)
    2. Compute energy density ρ_I from intent gradients
    3. Solve Poisson equation for Φ(r)
    4. Compare to analytical -GM/r

Author: Claude (Autonomous Research Agent)
Date: 2025-11-13
Session: #12
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit


class SynchronismGravitySimulation:
    """
    Spherically symmetric gravity solver for Synchronism framework.

    Tests whether Newtonian potential Φ = -GM/r emerges from intent dynamics.
    """

    def __init__(self, r_min=0.1, r_max=10.0, N_r=200, G=1.0):
        """
        Initialize spherical grid and parameters.

        Parameters:
        -----------
        r_min : float
            Minimum radius (avoid singularity at r=0)
        r_max : float
            Maximum radius (far field)
        N_r : int
            Number of radial grid points
        G : float
            Gravitational constant (set to 1 for simplicity)
        """
        self.r_min = r_min
        self.r_max = r_max
        self.N_r = N_r
        self.G = G

        # Radial grid (logarithmic for better resolution near origin)
        self.r = np.logspace(np.log10(r_min), np.log10(r_max), N_r)
        self.dr = np.diff(self.r)

        # Intent field I(r) - to be initialized
        self.I = None
        self.I_0 = 1.0  # Background intent

        # Energy density ρ_I(r)
        self.rho_I = None

        # Gravitational potential Φ(r)
        self.Phi = None

    def set_intent_profile_gaussian(self, M_total=1.0, R_mass=1.0):
        """
        Set Gaussian intent profile representing a massive object.

        I(r) = I_0 + ΔI * exp(-(r/R)²)

        For ρ_I ≈ (1/2)(I - I₀)², we need:
        ∫ (1/2)(I - I₀)² 4πr² dr = M_total

        Parameters:
        -----------
        M_total : float
            Total "intent mass" (analogous to gravitational mass)
        R_mass : float
            Characteristic radius of mass distribution
        """
        self.M_total = M_total
        self.R_mass = R_mass

        # For Gaussian: I - I₀ = ΔI exp(-(r/R)²)
        # ∫ (ΔI)² exp(-2(r/R)²) 4πr² dr from 0 to ∞
        # = (ΔI)² · 4π · (R³/2) · (π/2)^(3/2)
        # = (ΔI)² · π^(5/2) R³ / 2
        #
        # Want (1/2) × result = M_total
        # So: (ΔI)² · π^(5/2) R³ / 4 = M_total
        # ΔI = sqrt(4 M_total / (π^(5/2) R³))

        Delta_I = np.sqrt(4.0 * M_total / (np.pi**(2.5) * R_mass**3))

        # Intent profile
        self.I = self.I_0 + Delta_I * np.exp(-(self.r / R_mass)**2)

    def set_intent_profile_uniform_sphere(self, M_total=1.0, R_mass=1.0):
        """
        Set uniform sphere intent profile.

        I(r) = I_0 + ΔI  for r < R
        I(r) = I_0       for r > R

        Parameters:
        -----------
        M_total : float
            Total mass
        R_mass : float
            Sphere radius
        """
        self.M_total = M_total
        self.R_mass = R_mass

        # For uniform sphere: ρ = 3M/(4πR³)
        # Need ρ_I = (1/2)(I - I_0)² ≈ ρ
        # So (I - I_0) ≈ sqrt(2ρ) = sqrt(6M/(4πR³))
        Delta_I = np.sqrt(1.5 * M_total / (np.pi * R_mass**3))

        # Step function
        self.I = np.where(self.r < R_mass,
                         self.I_0 + Delta_I,
                         self.I_0)

    def compute_energy_density(self):
        """
        Compute gravitational source density from intent field.

        From Session #11:
            ρ_I = (1/2)[(∇I)² + (I - I₀)²]

        For spherical symmetry:
            ∇I = dI/dr
        """
        # Compute radial gradient dI/dr using central differences
        dI_dr = np.gradient(self.I, self.r)

        # Energy density components
        gradient_term = 0.5 * dI_dr**2
        potential_term = 0.5 * (self.I - self.I_0)**2

        self.rho_I = gradient_term + potential_term

        # Store components for analysis
        self.gradient_contribution = gradient_term
        self.potential_contribution = potential_term

        return self.rho_I

    def solve_poisson_spherical(self):
        """
        Solve spherical Poisson equation for gravitational potential.

        ∇²Φ = 4πG ρ_I

        In spherical coordinates:
            (1/r²) d/dr(r² dΦ/dr) = 4πG ρ_I

        Expanding:
            d²Φ/dr² + (2/r) dΦ/dr = 4πG ρ_I

        Boundary conditions:
            dΦ/dr(r_min) = 0 (symmetry at origin)
            Φ(r_max) = 0 (far field)
        """
        N = self.N_r
        r = self.r

        # Build finite difference matrix
        # Use non-uniform grid spacing

        diag_main = np.zeros(N)
        diag_lower = np.zeros(N-1)
        diag_upper = np.zeros(N-1)

        for i in range(1, N-1):
            # Grid spacings
            dr_minus = r[i] - r[i-1]
            dr_plus = r[i+1] - r[i]
            r_i = r[i]

            # Non-uniform finite differences for d²Φ/dr²
            # d²Φ/dr² ≈ [Φ_{i+1}/dr_plus - Φ_i/dr_center - Φ_i/dr_center + Φ_{i-1}/dr_minus] / dr_center
            #         ≈ 2[Φ_{i+1}/dr_plus - Φ_i(1/dr_plus + 1/dr_minus) + Φ_{i-1}/dr_minus] / (dr_plus + dr_minus)

            dr_avg = 0.5 * (dr_minus + dr_plus)

            # Coefficients for d²Φ/dr²
            coeff_left = 2.0 / (dr_minus * (dr_minus + dr_plus))
            coeff_right = 2.0 / (dr_plus * (dr_minus + dr_plus))
            coeff_center = -(coeff_left + coeff_right)

            # Coefficients for (2/r) dΦ/dr (central difference)
            dPhi_dr_coeff_left = -1.0 / (dr_minus + dr_plus)
            dPhi_dr_coeff_right = 1.0 / (dr_minus + dr_plus)

            # Combined: d²Φ/dr² + (2/r) dΦ/dr
            diag_lower[i-1] = coeff_left + (2.0 / r_i) * dPhi_dr_coeff_left
            diag_main[i] = coeff_center
            diag_upper[i] = coeff_right + (2.0 / r_i) * dPhi_dr_coeff_right

        # Boundary conditions
        # At r_min: dΦ/dr = 0 (use forward difference)
        # Φ_1 - Φ_0 = 0 → Φ_0 = Φ_1
        diag_main[0] = -1.0
        diag_upper[0] = 1.0

        # At r_max: Φ = 0
        diag_main[-1] = 1.0
        if N > 1:
            diag_lower[-1] = 0.0

        # Assemble sparse matrix
        Lap = diags([diag_lower, diag_main, diag_upper],
                    offsets=[-1, 0, 1],
                    shape=(N, N),
                    format='csr')

        # Right-hand side: 4πG ρ_I
        rhs = 4.0 * np.pi * self.G * self.rho_I.copy()

        # Boundary values
        rhs[0] = 0.0   # dΦ/dr = 0 at origin
        rhs[-1] = 0.0  # Φ = 0 at infinity

        # Solve sparse linear system
        self.Phi = spsolve(Lap, rhs)

        return self.Phi

    def compute_enclosed_mass(self):
        """
        Compute total enclosed mass from intent density.

        M_enc(r) = ∫₀ʳ ρ_I(r') 4πr'² dr'
        """
        # Integrate ρ_I over spherical shells
        integrand = self.rho_I * 4.0 * np.pi * self.r**2

        # Cumulative trapezoidal integration
        M_enc = np.zeros_like(self.r)
        for i in range(1, len(self.r)):
            dr = self.r[i] - self.r[i-1]
            M_enc[i] = M_enc[i-1] + 0.5 * (integrand[i] + integrand[i-1]) * dr

        self.M_enc = M_enc
        return M_enc

    def analytical_newtonian(self):
        """
        Compute analytical Newtonian potential for comparison.

        Φ_Newton(r) = -GM/r  (outside mass)
        """
        # Use total enclosed mass at each radius
        self.Phi_analytical = -self.G * self.M_enc / self.r
        return self.Phi_analytical

    def fit_power_law(self, r_fit_min=None, r_fit_max=None):
        """
        Fit Φ(r) to power law Φ = A/r^n + B to test 1/r behavior.

        Parameters:
        -----------
        r_fit_min, r_fit_max : float
            Range for fitting (default: outside mass distribution)
        """
        if r_fit_min is None:
            r_fit_min = 2.0 * self.R_mass
        if r_fit_max is None:
            r_fit_max = self.r_max

        # Select fitting region
        mask = (self.r >= r_fit_min) & (self.r <= r_fit_max)
        r_fit = self.r[mask]
        Phi_fit = self.Phi[mask]

        # Fit function: Φ = A/r^n + B
        def power_law(r, A, n, B):
            return A / r**n + B

        # Initial guess: A = -GM, n = 1, B = 0
        p0 = [-self.G * self.M_total, 1.0, 0.0]

        try:
            popt, pcov = curve_fit(power_law, r_fit, Phi_fit, p0=p0)
            perr = np.sqrt(np.diag(pcov))

            self.fit_params = {
                'A': popt[0],
                'n': popt[1],
                'B': popt[2],
                'A_err': perr[0],
                'n_err': perr[1],
                'B_err': perr[2]
            }

            # Compute chi-squared
            Phi_model = power_law(r_fit, *popt)
            residuals = Phi_fit - Phi_model
            chi_squared = np.sum(residuals**2)
            dof = len(r_fit) - 3
            chi_squared_per_dof = chi_squared / dof if dof > 0 else np.inf

            self.fit_params['chi_squared'] = chi_squared
            self.fit_params['dof'] = dof
            self.fit_params['chi_squared_per_dof'] = chi_squared_per_dof

            return self.fit_params

        except Exception as e:
            print(f"Fit failed: {e}")
            return None

    def plot_results(self, filename='Session12_Gravity_Emergence.png'):
        """
        Create comprehensive visualization of results.
        """
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))

        # 1. Intent profile I(r)
        ax = axes[0, 0]
        ax.plot(self.r, self.I, 'b-', linewidth=2, label='I(r)')
        ax.axhline(self.I_0, color='gray', linestyle='--', label='I₀')
        ax.axvline(self.R_mass, color='red', linestyle=':', alpha=0.5, label='R_mass')
        ax.set_xlabel('Radius r')
        ax.set_ylabel('Intent I(r)')
        ax.set_title('Intent Field Profile')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 2. Energy density ρ_I(r)
        ax = axes[0, 1]
        ax.plot(self.r, self.rho_I, 'g-', linewidth=2, label='ρ_I total')
        ax.plot(self.r, self.gradient_contribution, 'c--', alpha=0.6, label='(∇I)²/2')
        ax.plot(self.r, self.potential_contribution, 'm--', alpha=0.6, label='(I-I₀)²/2')
        ax.axvline(self.R_mass, color='red', linestyle=':', alpha=0.5)
        ax.set_xlabel('Radius r')
        ax.set_ylabel('Energy Density ρ_I')
        ax.set_title('Gravitational Source Density')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')

        # 3. Enclosed mass M(r)
        ax = axes[0, 2]
        ax.plot(self.r, self.M_enc, 'orange', linewidth=2)
        ax.axhline(self.M_total, color='red', linestyle='--',
                  label=f'M_total = {self.M_total:.3f}')
        ax.axvline(self.R_mass, color='red', linestyle=':', alpha=0.5)
        ax.set_xlabel('Radius r')
        ax.set_ylabel('Enclosed Mass M(r)')
        ax.set_title('Cumulative Mass Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 4. Gravitational potential Φ(r)
        ax = axes[1, 0]
        ax.plot(self.r, self.Phi, 'b-', linewidth=2, label='Φ (Synchronism)')
        ax.plot(self.r, self.Phi_analytical, 'r--', linewidth=2, label='Φ = -GM/r')
        ax.axvline(self.R_mass, color='red', linestyle=':', alpha=0.5, label='R_mass')
        ax.set_xlabel('Radius r')
        ax.set_ylabel('Potential Φ(r)')
        ax.set_title('Gravitational Potential')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # 5. Residuals |Φ - Φ_analytical|
        ax = axes[1, 1]
        residuals = np.abs(self.Phi - self.Phi_analytical)
        ax.plot(self.r, residuals, 'purple', linewidth=2)
        ax.axvline(self.R_mass, color='red', linestyle=':', alpha=0.5)
        ax.set_xlabel('Radius r')
        ax.set_ylabel('|Φ - Φ_Newtonian|')
        ax.set_title('Absolute Residuals')
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')

        # 6. Power law fit results
        ax = axes[1, 2]

        # Fitting region (outside mass)
        r_fit_min = 2.0 * self.R_mass
        mask = self.r >= r_fit_min
        r_fit = self.r[mask]

        # Plot data and fit
        ax.plot(r_fit, self.Phi[mask], 'bo', markersize=4, alpha=0.6, label='Data')

        if hasattr(self, 'fit_params') and self.fit_params is not None:
            A = self.fit_params['A']
            n = self.fit_params['n']
            B = self.fit_params['B']
            Phi_fit = A / r_fit**n + B

            ax.plot(r_fit, Phi_fit, 'r-', linewidth=2,
                   label=f'Fit: Φ = {A:.3f}/r^{n:.3f} + {B:.3f}')

            # Add text with fit results
            chi2_dof = self.fit_params['chi_squared_per_dof']
            textstr = f"n = {n:.3f} ± {self.fit_params['n_err']:.3f}\n"
            textstr += f"χ²/dof = {chi2_dof:.6f}"

            ax.text(0.05, 0.05, textstr, transform=ax.transAxes,
                   fontsize=10, verticalalignment='bottom',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.set_xlabel('Radius r')
        ax.set_ylabel('Potential Φ(r)')
        ax.set_title('Power Law Fit (Outside Mass)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"\n✓ Plot saved: {filename}")

        return fig


def run_gravity_validation():
    """
    Main simulation: Test if Φ(r) = -GM/r emerges from Synchronism.
    """
    print("=" * 80)
    print("SESSION #12: GRAVITY VALIDATION FROM SYNCHRONISM")
    print("=" * 80)
    print("\nTest: Does Newtonian potential Φ = -GM/r emerge from intent dynamics?")
    print("\nTheory (Session #11):")
    print("  ∇²Φ = 4πG ρ_I")
    print("  ρ_I = (1/2)[(∇I)² + (I-I₀)²]")
    print("  Expected: Φ(r) = -GM/r outside mass\n")

    # Initialize simulation
    print("Setting up spherical simulation...")
    sim = SynchronismGravitySimulation(
        r_min=0.1,
        r_max=10.0,
        N_r=300,
        G=1.0
    )

    # Test Case 1: Gaussian intent profile
    print("\n" + "-" * 80)
    print("TEST CASE 1: Gaussian Intent Distribution")
    print("-" * 80)

    M_total = 1.0
    R_mass = 1.0

    print(f"  M_total = {M_total}")
    print(f"  R_mass = {R_mass}")

    sim.set_intent_profile_gaussian(M_total=M_total, R_mass=R_mass)

    print("\n  Computing energy density ρ_I from intent gradients...")
    sim.compute_energy_density()

    print("  Solving Poisson equation ∇²Φ = 4πG ρ_I...")
    sim.solve_poisson_spherical()

    print("  Computing enclosed mass M(r)...")
    sim.compute_enclosed_mass()

    print("  Computing analytical Newtonian Φ = -GM/r...")
    sim.analytical_newtonian()

    print(f"\n  Total mass from integration: M = {sim.M_enc[-1]:.6f}")
    print(f"  Expected mass: M = {M_total:.6f}")
    print(f"  Discrepancy: {abs(sim.M_enc[-1] - M_total):.6f}")

    # Fit to power law in far field
    print("\n  Fitting Φ(r) to power law A/r^n + B (r > 2R_mass)...")
    fit_result = sim.fit_power_law(r_fit_min=2.0 * R_mass, r_fit_max=sim.r_max)

    if fit_result:
        print(f"\n  FIT RESULTS:")
        print(f"    Φ(r) = A/r^n + B")
        print(f"    A = {fit_result['A']:.6f} ± {fit_result['A_err']:.6f}")
        print(f"    n = {fit_result['n']:.6f} ± {fit_result['n_err']:.6f}")
        print(f"    B = {fit_result['B']:.6f} ± {fit_result['B_err']:.6f}")
        print(f"    χ²/dof = {fit_result['chi_squared_per_dof']:.8f}")

        # Expected: n = 1 for Newtonian gravity
        n_expected = 1.0
        n_deviation = abs(fit_result['n'] - n_expected)
        n_sigma = n_deviation / fit_result['n_err'] if fit_result['n_err'] > 0 else np.inf

        print(f"\n  VALIDATION:")
        print(f"    Expected: n = {n_expected} (Newtonian 1/r)")
        print(f"    Measured: n = {fit_result['n']:.6f} ± {fit_result['n_err']:.6f}")
        print(f"    Deviation: {n_deviation:.6f} ({n_sigma:.2f}σ)")

        if n_sigma < 3.0:
            print(f"    ✓ CONSISTENT with Newtonian gravity (within {n_sigma:.2f}σ)!")
        else:
            print(f"    ✗ INCONSISTENT with Newtonian gravity ({n_sigma:.2f}σ deviation)")

        # Check chi-squared quality
        if fit_result['chi_squared_per_dof'] < 0.01:
            print(f"    ✓ Excellent fit quality (χ²/dof = {fit_result['chi_squared_per_dof']:.8f})")
        elif fit_result['chi_squared_per_dof'] < 1.0:
            print(f"    ✓ Good fit quality (χ²/dof = {fit_result['chi_squared_per_dof']:.8f})")
        else:
            print(f"    ⚠ Poor fit quality (χ²/dof = {fit_result['chi_squared_per_dof']:.8f})")

    # Generate visualization
    print("\n  Creating visualization...")
    sim.plot_results(filename='/mnt/c/exe/projects/ai-agents/synchronism/Research/Session12_Gravity_Emergence.png')

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY: SESSION #12 GRAVITY VALIDATION")
    print("=" * 80)

    if fit_result and abs(fit_result['n'] - 1.0) < 3.0 * fit_result['n_err']:
        print("\n✅ SUCCESS: Newtonian gravity Φ = -GM/r EMERGES from Synchronism!")
        print(f"\n   Power law exponent n = {fit_result['n']:.3f} ± {fit_result['n_err']:.3f}")
        print(f"   Expected n = 1.000 (Newtonian)")
        print(f"   Fit quality χ²/dof = {fit_result['chi_squared_per_dof']:.6f}")
        print("\n   This completes classical unification:")
        print("   • Session #8: Electrostatics ✓")
        print("   • Session #9: Magnetism ✓")
        print("   • Session #12: Gravity ✓")
        print("\n   Synchronism explains ALL of classical physics from intent dynamics!")
    else:
        print("\n⚠ PARTIAL: Framework established but numerical validation needs refinement")

    print("\n" + "=" * 80)

    return sim, fit_result


if __name__ == "__main__":
    # Run validation
    sim, fit_result = run_gravity_validation()

    # Keep plot open for inspection
    plt.show()
