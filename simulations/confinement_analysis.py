#!/usr/bin/env python3
"""
Confinement Analysis Tools for SU(3) Lattice Gauge Theory

This module provides analysis tools for extracting confinement physics from
SU(3) lattice gauge simulations, specifically designed for Synchronism
Session #35 physics extraction.

Key observables:
- Wilson loops W(R,T) for multiple spatial (R) and temporal (T) extents
- Static quark-antiquark potential V(R) from Wilson loop ratios
- String tension σ from linear fit to V(R) = σR - α/R + C
- Comparison to QCD string tension σ_QCD ≈ 0.9 GeV/fm

Author: Autonomous Synchronism Research
Session: #33 (preparation for Session #35)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2
import sys


class ConfinementAnalyzer:
    """Analyzes Wilson loops to extract confinement physics from SU(3) lattice data."""

    def __init__(self, beta, lattice_spacing_fm=0.1):
        """
        Initialize confinement analyzer.

        Parameters:
        -----------
        beta : float
            Lattice coupling β = 6/g² for SU(3)
        lattice_spacing_fm : float
            Lattice spacing in femtometers (default: 0.1 fm)
        """
        self.beta = beta
        self.a = lattice_spacing_fm  # Lattice spacing in fm
        self.sigma_QCD_GeV_fm = 0.9  # QCD string tension in GeV/fm

    def measure_wilson_loop(self, U, R, T, mu=0, nu=1):
        """
        Measure Wilson loop W(R,T) for given spatial (R) and temporal (T) extent.

        A Wilson loop is a closed rectangular path of SU(3) links:
        - Spatial extent R in direction mu
        - Temporal extent T in direction nu

        For confinement measurement:
        - R = quark-antiquark separation
        - T = temporal extent (take T → ∞ limit)

        Parameters:
        -----------
        U : ndarray
            SU(3) link variables, shape (Nx, Ny, Nz, Nt, 4, 3, 3)
        R : int
            Spatial extent of Wilson loop
        T : int
            Temporal extent of Wilson loop
        mu : int
            Spatial direction (default: 0 for x-direction)
        nu : int
            Temporal direction (default: 1 for y or t-direction)

        Returns:
        --------
        W : float
            Wilson loop expectation value (real part of trace)
        """
        Nx, Ny, Nz, Nt = U.shape[:4]

        # Average over all starting positions
        wilson_sum = 0.0
        count = 0

        # Loop over all starting positions
        for x in range(Nx):
            for y in range(Ny):
                for z in range(Nz):
                    for t in range(Nt):
                        # Construct R×T Wilson loop starting at (x,y,z,t)
                        W_loop = np.eye(3, dtype=complex)

                        # Forward path in mu direction (R steps)
                        pos = [x, y, z, t]
                        for step in range(R):
                            W_loop = W_loop @ U[tuple(pos + [mu])]
                            pos[mu] = (pos[mu] + 1) % [Nx, Ny, Nz, Nt][mu]

                        # Forward path in nu direction (T steps)
                        for step in range(T):
                            W_loop = W_loop @ U[tuple(pos + [nu])]
                            pos[nu] = (pos[nu] + 1) % [Nx, Ny, Nz, Nt][nu]

                        # Backward path in mu direction (R steps)
                        for step in range(R):
                            pos[mu] = (pos[mu] - 1) % [Nx, Ny, Nz, Nt][mu]
                            W_loop = W_loop @ U[tuple(pos + [mu])].conj().T

                        # Backward path in nu direction (T steps)
                        for step in range(T):
                            pos[nu] = (pos[nu] - 1) % [Nx, Ny, Nz, Nt][nu]
                            W_loop = W_loop @ U[tuple(pos + [nu])].conj().T

                        # Take trace and accumulate
                        wilson_sum += np.real(np.trace(W_loop)) / 3.0  # Normalize by 3 for SU(3)
                        count += 1

        return wilson_sum / count

    def extract_potential_from_wilson_loops(self, wilson_loops):
        """
        Extract static quark-antiquark potential V(R) from Wilson loop measurements.

        The potential is extracted from:
            W(R,T) ∝ exp(-V(R) * T)

        Therefore:
            V(R) = -(1/a) * (1/T) * ln(W(R,T) / W(R,T-1))

        where a is the lattice spacing.

        Parameters:
        -----------
        wilson_loops : dict
            Dictionary of Wilson loop measurements
            Keys: (R, T) tuples
            Values: W(R,T) expectation values

        Returns:
        --------
        R_values : ndarray
            Array of separation distances R
        V_values : ndarray
            Array of potential values V(R) in GeV
        V_errors : ndarray
            Statistical errors on V(R)
        """
        # Organize data by R
        R_dict = {}
        for (R, T), W in wilson_loops.items():
            if R not in R_dict:
                R_dict[R] = []
            R_dict[R].append((T, W))

        R_values = []
        V_values = []
        V_errors = []

        # Extract V(R) for each R using adjacent T values
        for R in sorted(R_dict.keys()):
            measurements = sorted(R_dict[R])  # Sort by T

            if len(measurements) < 2:
                continue

            # Use largest T values for best signal (ground state dominance)
            T1, W1 = measurements[-2]
            T2, W2 = measurements[-1]

            if W1 > 0 and W2 > 0:
                # V(R) = -(1/a) * ln(W(R,T2) / W(R,T1)) / (T2 - T1)
                V = -(1.0 / self.a) * np.log(W2 / W1) / (T2 - T1)

                # Error estimate (simplified, assumes Gaussian)
                # δV ≈ (1/a) * sqrt((δW/W)² + (δW'/W')²) / ΔT
                # For now, use 10% uncertainty placeholder
                V_err = 0.1 * abs(V)

                R_values.append(R)
                V_values.append(V)
                V_errors.append(V_err)

        return np.array(R_values), np.array(V_values), np.array(V_errors)

    def fit_linear_plus_coulomb(self, R, V, V_err):
        """
        Fit potential to confinement form: V(R) = σR - α/R + C

        This is the expected form for QCD confinement:
        - σR: Linear confining term (string tension)
        - -α/R: Coulombic term (one-gluon exchange)
        - C: Constant offset

        Parameters:
        -----------
        R : ndarray
            Separation distances (in lattice units)
        V : ndarray
            Potential values (in GeV)
        V_err : ndarray
            Errors on potential values

        Returns:
        --------
        params : tuple
            Best-fit parameters (σ, α, C)
        param_errors : tuple
            Errors on parameters
        chi2_dof : float
            Chi-squared per degree of freedom
        """
        def confinement_potential(R, sigma, alpha, C):
            return sigma * R - alpha / R + C

        # Initial guess: σ ≈ 0.9 GeV/fm, α ≈ 0.5, C ≈ 0
        p0 = [0.9, 0.5, 0.0]

        try:
            # Perform weighted fit
            popt, pcov = curve_fit(
                confinement_potential,
                R * self.a,  # Convert to physical units (fm)
                V,
                p0=p0,
                sigma=V_err,
                absolute_sigma=True,
                maxfev=10000
            )

            # Extract parameters and errors
            sigma, alpha, C = popt
            sigma_err, alpha_err, C_err = np.sqrt(np.diag(pcov))

            # Calculate chi-squared
            V_fit = confinement_potential(R * self.a, sigma, alpha, C)
            chi2_value = np.sum(((V - V_fit) / V_err) ** 2)
            dof = len(R) - 3  # 3 parameters
            chi2_dof = chi2_value / dof if dof > 0 else np.inf

            return (sigma, alpha, C), (sigma_err, alpha_err, C_err), chi2_dof

        except Exception as e:
            print(f"Fit failed: {e}")
            return (np.nan, np.nan, np.nan), (np.nan, np.nan, np.nan), np.inf

    def fit_pure_linear(self, R, V, V_err):
        """
        Fit potential to pure linear form: V(R) = σR + C

        This tests whether confinement is present without Coulomb term.

        Parameters:
        -----------
        R : ndarray
            Separation distances (in lattice units)
        V : ndarray
            Potential values (in GeV)
        V_err : ndarray
            Errors on potential values

        Returns:
        --------
        params : tuple
            Best-fit parameters (σ, C)
        param_errors : tuple
            Errors on parameters
        chi2_dof : float
            Chi-squared per degree of freedom
        """
        def linear_potential(R, sigma, C):
            return sigma * R + C

        p0 = [0.9, 0.0]

        try:
            popt, pcov = curve_fit(
                linear_potential,
                R * self.a,
                V,
                p0=p0,
                sigma=V_err,
                absolute_sigma=True
            )

            sigma, C = popt
            sigma_err, C_err = np.sqrt(np.diag(pcov))

            V_fit = linear_potential(R * self.a, sigma, C)
            chi2_value = np.sum(((V - V_fit) / V_err) ** 2)
            dof = len(R) - 2
            chi2_dof = chi2_value / dof if dof > 0 else np.inf

            return (sigma, C), (sigma_err, C_err), chi2_dof

        except Exception as e:
            print(f"Fit failed: {e}")
            return (np.nan, np.nan), (np.nan, np.nan), np.inf

    def validate_confinement(self, sigma, sigma_err):
        """
        Validate whether measured string tension indicates confinement.

        Validation criteria:
        1. σ > 0 (necessary for confinement)
        2. σ within factor 2-3 of QCD value (0.9 GeV/fm)
        3. Statistical significance: σ/σ_err > 2 (2-sigma detection)

        Parameters:
        -----------
        sigma : float
            Measured string tension (GeV/fm)
        sigma_err : float
            Error on string tension

        Returns:
        --------
        validated : bool
            True if confinement validated
        message : str
            Validation message
        """
        # Check if σ > 0
        if sigma <= 0:
            return False, f"String tension σ = {sigma:.3f} ≤ 0 (no confinement)"

        # Check statistical significance
        if sigma / sigma_err < 2.0:
            return False, f"String tension σ/σ_err = {sigma/sigma_err:.2f} < 2 (low significance)"

        # Check agreement with QCD
        ratio = sigma / self.sigma_QCD_GeV_fm
        if ratio < 0.33 or ratio > 3.0:
            return False, f"String tension ratio σ/σ_QCD = {ratio:.2f} outside [0.33, 3.0]"

        # All criteria met!
        return True, f"✅ Confinement validated: σ = {sigma:.3f}±{sigma_err:.3f} GeV/fm (QCD: {self.sigma_QCD_GeV_fm} GeV/fm)"

    def plot_potential(self, R, V, V_err, fit_params_full, fit_params_linear, output_file):
        """
        Create publication-quality plot of V(R) with fits.

        Parameters:
        -----------
        R : ndarray
            Separation distances (lattice units)
        V : ndarray
            Potential values (GeV)
        V_err : ndarray
            Errors on potential
        fit_params_full : tuple
            Parameters from full fit (σ, α, C)
        fit_params_linear : tuple
            Parameters from linear fit (σ, C)
        output_file : str
            Path to save figure
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        # Physical separation in fm
        R_fm = R * self.a

        # --- Left panel: V(R) with both fits ---
        ax1.errorbar(R_fm, V, yerr=V_err, fmt='o', color='blue',
                     markersize=8, capsize=5, label='Synchronism Data')

        # Plot fits
        R_fine = np.linspace(R_fm.min(), R_fm.max(), 200)

        # Full fit: σR - α/R + C
        sigma, alpha, C = fit_params_full[0]
        if not np.isnan(sigma):
            V_full = sigma * R_fine - alpha / R_fine + C
            ax1.plot(R_fine, V_full, '-', color='red', linewidth=2,
                     label=f'Linear+Coulomb: σ={sigma:.3f} GeV/fm')

        # Linear fit: σR + C
        sigma_lin, C_lin = fit_params_linear[0]
        if not np.isnan(sigma_lin):
            V_linear = sigma_lin * R_fine + C_lin
            ax1.plot(R_fine, V_linear, '--', color='green', linewidth=2,
                     label=f'Pure Linear: σ={sigma_lin:.3f} GeV/fm')

        # QCD reference
        ax1.axhline(y=0, color='black', linestyle=':', linewidth=1, alpha=0.5)

        ax1.set_xlabel('Quark Separation R (fm)', fontsize=14)
        ax1.set_ylabel('Potential V(R) (GeV)', fontsize=14)
        ax1.set_title('SU(3) Static Quark-Antiquark Potential', fontsize=16)
        ax1.legend(fontsize=11, loc='best')
        ax1.grid(alpha=0.3)

        # --- Right panel: Linearized plot R*V(R) vs R ---
        # For pure linear potential, R*V should be quadratic
        # For Coulomb, R*V should be linear
        ax2.errorbar(R_fm, R_fm * V, yerr=R_fm * V_err, fmt='s', color='blue',
                     markersize=8, capsize=5, label='R × V(R) Data')

        if not np.isnan(sigma):
            RV_full = R_fine * (sigma * R_fine - alpha / R_fine + C)
            ax2.plot(R_fine, RV_full, '-', color='red', linewidth=2,
                     label='Linear+Coulomb Fit')

        if not np.isnan(sigma_lin):
            RV_linear = R_fine * (sigma_lin * R_fine + C_lin)
            ax2.plot(R_fine, RV_linear, '--', color='green', linewidth=2,
                     label='Pure Linear Fit')

        ax2.set_xlabel('Quark Separation R (fm)', fontsize=14)
        ax2.set_ylabel('R × V(R) (GeV·fm)', fontsize=14)
        ax2.set_title('Linearized Potential (Tests Confinement)', fontsize=16)
        ax2.legend(fontsize=11, loc='best')
        ax2.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved plot to {output_file}")
        plt.close()

    def generate_report(self, R, V, V_err, fit_full, fit_linear):
        """
        Generate comprehensive text report of confinement analysis.

        Parameters:
        -----------
        R : ndarray
            Separation distances
        V : ndarray
            Potential values
        V_err : ndarray
            Errors
        fit_full : tuple
            (params, errors, chi2_dof) from full fit
        fit_linear : tuple
            (params, errors, chi2_dof) from linear fit

        Returns:
        --------
        report : str
            Formatted report text
        """
        report = []
        report.append("=" * 80)
        report.append("SU(3) CONFINEMENT ANALYSIS REPORT")
        report.append("Synchronism Research Session #35")
        report.append("=" * 80)
        report.append("")

        # Lattice parameters
        report.append("LATTICE PARAMETERS:")
        report.append(f"  β = {self.beta:.4f}")
        report.append(f"  Lattice spacing a = {self.a:.4f} fm")
        report.append(f"  QCD reference: σ_QCD ≈ {self.sigma_QCD_GeV_fm} GeV/fm")
        report.append("")

        # Data summary
        report.append("POTENTIAL MEASUREMENTS:")
        report.append("  R (fm)    V(R) (GeV)    Error")
        report.append("  " + "-" * 40)
        for i in range(len(R)):
            report.append(f"  {R[i]*self.a:6.3f}    {V[i]:8.4f}      {V_err[i]:6.4f}")
        report.append("")

        # Full fit results
        report.append("FIT 1: LINEAR + COULOMB  V(R) = σR - α/R + C")
        report.append("  " + "-" * 60)
        params, errors, chi2_dof = fit_full
        sigma, alpha, C = params
        sigma_err, alpha_err, C_err = errors

        report.append(f"  String tension σ  = {sigma:.4f} ± {sigma_err:.4f} GeV/fm")
        report.append(f"  Coulomb term α    = {alpha:.4f} ± {alpha_err:.4f} GeV")
        report.append(f"  Constant C        = {C:.4f} ± {C_err:.4f} GeV")
        report.append(f"  χ²/dof            = {chi2_dof:.4f}")

        if not np.isnan(sigma):
            ratio = sigma / self.sigma_QCD_GeV_fm
            report.append(f"  σ/σ_QCD           = {ratio:.3f}")
        report.append("")

        # Linear fit results
        report.append("FIT 2: PURE LINEAR  V(R) = σR + C")
        report.append("  " + "-" * 60)
        params_lin, errors_lin, chi2_dof_lin = fit_linear
        sigma_lin, C_lin = params_lin
        sigma_lin_err, C_lin_err = errors_lin

        report.append(f"  String tension σ  = {sigma_lin:.4f} ± {sigma_lin_err:.4f} GeV/fm")
        report.append(f"  Constant C        = {C_lin:.4f} ± {C_lin_err:.4f} GeV")
        report.append(f"  χ²/dof            = {chi2_dof_lin:.4f}")

        if not np.isnan(sigma_lin):
            ratio_lin = sigma_lin / self.sigma_QCD_GeV_fm
            report.append(f"  σ/σ_QCD           = {ratio_lin:.3f}")
        report.append("")

        # Validation
        report.append("CONFINEMENT VALIDATION:")
        report.append("  " + "-" * 60)

        # Check full fit
        validated, message = self.validate_confinement(sigma, sigma_err)
        report.append(f"  Full fit: {message}")

        # Check linear fit
        validated_lin, message_lin = self.validate_confinement(sigma_lin, sigma_lin_err)
        report.append(f"  Linear fit: {message_lin}")
        report.append("")

        # Overall conclusion
        report.append("CONCLUSION:")
        report.append("  " + "-" * 60)
        if validated or validated_lin:
            report.append("  ✅ CONFINEMENT VALIDATED")
            report.append("  Linear potential V(R) ∝ R observed with σ > 0")
            report.append("  String tension consistent with QCD")
            report.append("  → Synchronism successfully derives strong force confinement!")
        else:
            report.append("  ❌ CONFINEMENT NOT VALIDATED")
            report.append("  Either σ ≤ 0, low significance, or inconsistent with QCD")
            report.append("  → Synchronism requires refinement for strong force")

        report.append("=" * 80)

        return "\n".join(report)


def example_analysis():
    """
    Example usage of ConfinementAnalyzer with synthetic data.

    This demonstrates the complete analysis workflow for Session #35.
    """
    print("Confinement Analysis Example")
    print("=" * 80)

    # Create analyzer
    analyzer = ConfinementAnalyzer(beta=5.7, lattice_spacing_fm=0.1)

    # Synthetic Wilson loop data (for demonstration)
    # In Session #35, this will come from actual SU(3) simulation
    print("Generating synthetic Wilson loop data...")

    # Expected behavior: W(R,T) ≈ exp(-V(R)*T) where V(R) = σR - α/R + C
    sigma_true = 0.85  # GeV/fm (close to QCD)
    alpha_true = 0.3   # GeV
    C_true = 0.1       # GeV

    wilson_loops = {}
    for R in range(2, 8):  # R = 2 to 7 lattice units
        for T in range(3, 7):  # T = 3 to 6 lattice units
            R_fm = R * analyzer.a
            V_R = sigma_true * R_fm - alpha_true / R_fm + C_true
            # Add noise
            noise = np.random.normal(0, 0.05)
            W = np.exp(-V_R * T * analyzer.a) * (1 + noise)
            wilson_loops[(R, T)] = max(W, 1e-10)  # Avoid log(0)

    print(f"Generated {len(wilson_loops)} Wilson loop measurements")
    print()

    # Extract potential
    print("Extracting potential from Wilson loops...")
    R, V, V_err = analyzer.extract_potential_from_wilson_loops(wilson_loops)
    print(f"Extracted V(R) for {len(R)} separation distances")
    print()

    # Fit potentials
    print("Fitting confinement models...")
    fit_full = analyzer.fit_linear_plus_coulomb(R, V, V_err)
    fit_linear = analyzer.fit_pure_linear(R, V, V_err)
    print()

    # Generate report
    report = analyzer.generate_report(R, V, V_err, fit_full, fit_linear)
    print(report)
    print()

    # Create plot
    output_file = "confinement_analysis_example.png"
    analyzer.plot_potential(R, V, V_err, fit_full, fit_linear, output_file)

    print("Example analysis complete!")
    print()
    print("For Session #35:")
    print("1. Run SU(3) simulation to generate wilson_loops dictionary")
    print("2. Use this analyzer to extract V(R)")
    print("3. Fit to confinement form")
    print("4. Validate against QCD")


if __name__ == "__main__":
    example_analysis()
