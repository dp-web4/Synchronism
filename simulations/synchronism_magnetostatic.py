"""
Synchronism Magnetostatic Potential - Session #9

Testing if magnetic interactions emerge from Synchronism.

From Session #9 derivation:
  Static vector potential equation: ∇²A = -j

This is analogous to Session #8 electrostatic: ∇²φ = -ρ

Expected result:
- A(r) ∝ ln(r) in 2D (like φ in Session #8)
- Force between currents F ∝ 1/R (magnetic Coulomb)

Method:
1. Solve vector Poisson equation for A_z(x,y)
2. Two parallel current sources (in z direction)
3. Measure A_z vs separation R
4. Check if F ∝ 1/R emerges
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import os


class SynchronismMagnetostaticSolver:
    """
    Solve magnetostatic Synchronism equation on 2D lattice.

    Equation (from Session #9 derivation):
      ∇²A = -j

    For currents in z direction with (x,y) symmetry:
      ∇²A_z(x,y) = -j_z(x,y)

    This is identical structure to Session #8 electrostatic!
    """

    def __init__(self, Lx, Ly, alpha=0.0, beta=0.0):
        self.Lx = Lx
        self.Ly = Ly
        self.alpha = alpha  # Intent-scalar coupling (not used in magnetostatic)
        self.beta = beta    # Intent-vector coupling (for future)

        # State
        self.A_z = np.zeros((Lx, Ly))  # Vector potential (z component)
        self.I = np.ones((Lx, Ly))     # Intent field (for future intent coupling)
        self.j_z = np.zeros((Lx, Ly))  # Current density (z component)

    def set_current_sources(self, currents_and_positions):
        """Place current sources on lattice

        Args:
            currents_and_positions: List of (I, x, y) tuples
                I: current strength
                x, y: position on lattice
        """
        self.j_z.fill(0.0)
        for I_current, x, y in currents_and_positions:
            # Gaussian current distribution (smoother than delta function)
            for i in range(self.Lx):
                for j in range(self.Ly):
                    r2 = (i-x)**2 + (j-y)**2
                    self.j_z[i,j] += I_current * np.exp(-r2/2.0) / (2*np.pi)

    def build_laplacian_matrix(self):
        """
        Build discrete Laplacian matrix for 2D lattice with periodic BC.

        Returns sparse matrix such that: (Lap @ u.flatten()) gives ∇²u

        Identical to Session #8 implementation.
        """
        N = self.Lx * self.Ly

        # Main diagonal: -4
        # Off-diagonals: +1 for neighbors
        main = -4 * np.ones(N)
        off_x = np.ones(N)
        off_y = np.ones(N)

        # Periodic boundary handling
        for i in range(self.Lx):
            for j in range(self.Ly):
                idx = i * self.Ly + j

                # x boundaries
                if i == 0 or i == self.Lx - 1:
                    off_x[idx] = 1  # Periodic

                # y boundaries
                if j == 0 or j == self.Ly - 1:
                    off_y[idx] = 1  # Periodic

        # Build matrix
        Lap = diags([main, off_x, off_x, off_y, off_y],
                    [0, -self.Ly, self.Ly, -1, 1],
                    shape=(N, N), format='csr')

        return Lap

    def solve_magnetostatic(self):
        """
        Solve magnetostatic equation:
          ∇²A_z = -j_z

        This is a direct solve (no iteration needed for pure magnetostatic).
        """
        Lap = self.build_laplacian_matrix()
        N = self.Lx * self.Ly

        # Flatten fields
        j_flat = self.j_z.flatten()

        # Solve: ∇²A = -j
        # Lap @ A = -j
        A_flat = spsolve(Lap, -j_flat)

        # Reshape back
        self.A_z = A_flat.reshape((self.Lx, self.Ly))

        print(f"  Solved magnetostatic equation")

    def measure_vector_potential_difference(self, pos1, pos2):
        """Measure vector potential difference between two positions

        For two parallel currents I₁, I₂ at positions r₁, r₂:
          Interaction energy: U = I₁ I₂ (A(r₁) - A(r₂))
          Force: F = I₁ I₂ dA/dr

        We measure A difference to get interaction strength.
        """
        x1, y1 = pos1
        x2, y2 = pos2

        # Vector potential difference
        A_diff = self.A_z[x1,y1] - self.A_z[x2,y2]

        # Get currents at those positions (from j_z)
        I1 = np.sum(self.j_z[max(0,x1-1):x1+2, max(0,y1-1):y1+2])
        I2 = np.sum(self.j_z[max(0,x2-1):x2+2, max(0,y2-1):y2+2])

        # Interaction energy
        U = I1 * I2 * abs(A_diff)

        return U


def run_magnetostatic_scan(R_values, Lx=64, Ly=64, verbose=True):
    """
    Measure magnetic interaction U(R) for various separations.

    Args:
        R_values: Separation distances to test
        Lx, Ly: Lattice size

    Returns:
        dict with R, U (interaction energy)
    """
    U_results = []

    if verbose:
        print(f"\n{'='*70}")
        print(f"Synchronism Magnetostatic Test - Session #9")
        print(f"{'='*70}")
        print(f"Lattice: {Lx}×{Ly}")
        print(f"Testing {len(R_values)} separation values")
        print(f"Expected: U(R) ∝ ln(R) in 2D, or U(R) ∝ 1/R in 3D")

    for R in R_values:
        if verbose:
            print(f"\n  R = {R}:")

        # Create solver
        solver = SynchronismMagnetostaticSolver(Lx, Ly)

        # Place currents (parallel, in z direction)
        cx, cy = Lx // 2, Ly // 2
        pos1 = (cx - R//2, cy)
        pos2 = (cx + R//2, cy)

        solver.set_current_sources([
            (+1.0, pos1[0], pos1[1]),  # Current +I
            (-1.0, pos2[0], pos2[1]),  # Current -I (antiparallel for interaction)
        ])

        if verbose:
            print(f"    Currents: +I at {pos1}, -I at {pos2}")
            print(f"    Solving ∇²A = -j...", end='', flush=True)

        # Solve magnetostatic equation
        solver.solve_magnetostatic()

        if verbose:
            print(" done")

        # Measure interaction energy
        U = solver.measure_vector_potential_difference(pos1, pos2)
        U_results.append(U)

        if verbose:
            print(f"    U(R={R}) = {U:.6f}")

    return {
        'R': np.array(R_values),
        'U': np.array(U_results),
    }


def fit_magnetic_interaction(R, U):
    """Fit to U = A*ln(R) + B (2D) or U = -A/R + B (3D-like)"""

    # Try logarithmic (2D expected)
    def log_form(R, A, B):
        return A * np.log(R) + B

    # Try inverse (3D-like, or if simulation has 3D-like behavior)
    def inverse_form(R, A, B):
        return -A / R + B

    results = {}

    # Fit logarithmic
    try:
        mask = R > 2
        R_fit = R[mask]
        U_fit = U[mask]

        popt_log, pcov_log = curve_fit(log_form, R_fit, U_fit, p0=[1.0, 0.0])
        perr_log = np.sqrt(np.diag(pcov_log))

        A_log, B_log = popt_log
        dA_log, dB_log = perr_log

        U_model_log = log_form(R_fit, A_log, B_log)
        chi2_log = np.sum((U_fit - U_model_log)**2) / len(R_fit)

        results['log'] = {
            'A': A_log, 'dA': dA_log,
            'B': B_log, 'dB': dB_log,
            'chi2': chi2_log,
            'U_fit': log_form(R, A_log, B_log)
        }
    except Exception as e:
        print(f"Log fit failed: {e}")
        results['log'] = None

    # Fit inverse
    try:
        mask = R > 2
        R_fit = R[mask]
        U_fit = U[mask]

        popt_inv, pcov_inv = curve_fit(inverse_form, R_fit, U_fit, p0=[1.0, 0.0])
        perr_inv = np.sqrt(np.diag(pcov_inv))

        A_inv, B_inv = popt_inv
        dA_inv, dB_inv = perr_inv

        U_model_inv = inverse_form(R_fit, A_inv, B_inv)
        chi2_inv = np.sum((U_fit - U_model_inv)**2) / len(R_fit)

        results['inverse'] = {
            'A': A_inv, 'dA': dA_inv,
            'B': B_inv, 'dB': dB_inv,
            'chi2': chi2_inv,
            'U_fit': inverse_form(R, A_inv, B_inv)
        }
    except Exception as e:
        print(f"Inverse fit failed: {e}")
        results['inverse'] = None

    return results


def plot_magnetostatic_results(results, fit_results, filename=None):
    """Plot U(R) with both logarithmic and inverse fits"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    R = results['R']
    U = results['U']

    # Left: U vs R
    ax1.plot(R, U, 'o-', markersize=8, linewidth=2,
             label='Synchronism magnetostatic', color='blue')

    if fit_results.get('log'):
        ax1.plot(R, fit_results['log']['U_fit'], '--', linewidth=2,
                label=f"Log fit: A·ln(R)+B", color='red')

    if fit_results.get('inverse'):
        ax1.plot(R, fit_results['inverse']['U_fit'], '-.', linewidth=2,
                label=f"Inverse fit: -A/R+B", color='green')

    # Add fit info
    info_text = ""
    if fit_results.get('log'):
        info_text += f"Log: χ² = {fit_results['log']['chi2']:.4f}\n"
    if fit_results.get('inverse'):
        info_text += f"Inv: χ² = {fit_results['inverse']['chi2']:.4f}\n"

    if info_text:
        ax1.text(0.05, 0.95, info_text.strip(),
                transform=ax1.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                fontsize=10)

    ax1.set_xlabel('Separation R (lattice units)', fontsize=12)
    ax1.set_ylabel('Interaction Energy U(R)', fontsize=12)
    ax1.set_title('Magnetic Interaction from Synchronism', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Right: Test both linearizations
    # Top half: U vs ln(R)
    # Bottom half: U vs 1/R
    ax2_top = ax2

    ln_R = np.log(R[R > 1])
    U_ln = U[R > 1]

    ax2_top.plot(ln_R, U_ln, 'o-', markersize=8, linewidth=2, color='blue')
    if fit_results.get('log'):
        ln_R_full = np.log(R)
        ax2_top.plot(ln_R_full, fit_results['log']['U_fit'], '--',
                    linewidth=2, color='red', label='Log fit')

    ax2_top.set_xlabel('ln(R)', fontsize=12)
    ax2_top.set_ylabel('U(R)', fontsize=12)
    ax2_top.set_title('Linearized: U vs ln(R) [2D expected]', fontsize=12)
    ax2_top.legend(fontsize=10)
    ax2_top.grid(True, alpha=0.3)

    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {filename}")
    else:
        plt.show()

    plt.close()


if __name__ == "__main__":
    print("\n" + "="*70)
    print("SYNCHRONISM MAGNETOSTATIC TEST - Session #9")
    print("="*70)
    print("\nBased on Session #9 magnetism derivation:")
    print("  ∇²A = -j  (magnetostatic equation)")
    print("\nExpected in 2D:")
    print("  A(r) ∝ ln(r)")
    print("  U(R) ∝ ln(R) for current interaction")
    print("\nThis is magnetic analog of Session #8 Coulomb test!")

    # Parameters
    R_values = [3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20]

    print(f"\nTesting {len(R_values)} separations: {R_values}")

    # Run measurement
    results = run_magnetostatic_scan(R_values, Lx=64, Ly=64, verbose=True)

    # Fit to both forms
    print(f"\n{'='*70}")
    print("FITTING TO BOTH FORMS")
    print("="*70)

    fit_results = fit_magnetic_interaction(results['R'], results['U'])

    if fit_results.get('log'):
        print(f"\nLogarithmic fit: U(R) = A·ln(R) + B")
        print(f"  A = {fit_results['log']['A']:.4f} ± {fit_results['log']['dA']:.4f}")
        print(f"  B = {fit_results['log']['B']:.4f} ± {fit_results['log']['dB']:.4f}")
        print(f"  χ² = {fit_results['log']['chi2']:.6f}")

    if fit_results.get('inverse'):
        print(f"\nInverse fit: U(R) = -A/R + B")
        print(f"  A = {fit_results['inverse']['A']:.4f} ± {fit_results['inverse']['dA']:.4f}")
        print(f"  B = {fit_results['inverse']['B']:.4f} ± {fit_results['inverse']['dB']:.4f}")
        print(f"  χ² = {fit_results['inverse']['chi2']:.6f}")

    # Determine which fits better
    if fit_results.get('log') and fit_results.get('inverse'):
        chi2_log = fit_results['log']['chi2']
        chi2_inv = fit_results['inverse']['chi2']

        print(f"\n{'='*70}")
        if chi2_log < chi2_inv:
            print("RESULT: Logarithmic fit is better (2D behavior confirmed!)")
        else:
            print("RESULT: Inverse fit is better (3D-like behavior)")
        print(f"{'='*70}")

    # Plot
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'Research')
    plot_filename = os.path.join(output_dir, 'Session9_Magnetostatic_Emergence.png')
    plot_magnetostatic_results(results, fit_results, filename=plot_filename)

    print("\n" + "="*70)
    print("SESSION #9 MAGNETOSTATIC TEST COMPLETE")
    print("="*70)
    print("\nKey findings:")
    print("  1. Magnetostatic equation ∇²A = -j solved successfully")
    print("  2. Interaction U(R) measured between current sources")
    print("  3. Fitted to both ln(R) [2D] and 1/R [3D-like] forms")
    print("  4. Tests if magnetic forces emerge from Synchronism")
    print("\nComparison to Session #8:")
    print("  - Session #8: ∇²φ = -ρ → V ∝ 1/R (Coulomb)")
    print("  - Session #9: ∇²A = -j → U ∝ ? (magnetic analog)")
    print("\n" + "="*70)
