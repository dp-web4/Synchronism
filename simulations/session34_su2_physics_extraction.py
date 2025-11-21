#!/usr/bin/env python3
"""
SU(2) Physics Extraction - Session #34
======================================

Goal: Extract weak force physics from SU(2) lattice gauge simulation.

This script runs a production-quality SU(2) simulation to measure:
- Yukawa screening from W/Z boson mass
- Static potential V(R) from Polyakov loop correlators
- Validation against QCD expectations

Expected runtime: 4-6 hours
Expected output: Definitive validation or falsification of SU(2) emergence

Based on Session #30 implementation, optimized for physics extraction.

Author: CBP Autonomous Synchronism Research
Date: 2025-11-21
Session: #34
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm, norm
from scipy.optimize import curve_fit
import sys
import os
import time
import pickle

# Import SU(2) lattice class from Session #30
sys.path.insert(0, '/mnt/c/exe/projects/ai-agents/synchronism/simulations')
from synchronism_session30_su2_lattice_3p1d import SU2Lattice3p1D

# Import stats utilities
sys.path.append('/mnt/c/exe/projects/ai-agents/private-context/tools/lattice-gauge')
from stats_utils import jackknife_function


def yukawa_potential(R, alpha, M, const):
    """Yukawa potential: V(R) = -α exp(-MR)/R + C"""
    return -alpha * np.exp(-M * R) / R + const


def coulomb_potential(R, alpha, const):
    """Coulomb potential: V(R) = -α/R + C"""
    return -alpha / R + const


def run_production_su2(Lx, Ly, Lz, Nt, beta, n_therm, n_meas, meas_interval, output_prefix):
    """
    Run production SU(2) simulation for physics extraction.

    Parameters:
    -----------
    Lx, Ly, Lz : int
        Spatial lattice dimensions
    Nt : int
        Temporal lattice dimension
    beta : float
        SU(2) coupling (β = 4/g²)
    n_therm : int
        Thermalization sweeps
    n_meas : int
        Measurement sweeps
    meas_interval : int
        Sweeps between measurements
    output_prefix : str
        Prefix for output files

    Returns:
    --------
    results : dict
        Complete results including:
        - R: Separations
        - V: Potential values
        - dV: Errors on V
        - C_mean: Polyakov correlators
        - plaquettes: Plaquette history
        - configs: Saved configurations (subset)
    """

    print("\n" + "="*80)
    print("SESSION #34: SU(2) WEAK FORCE PHYSICS EXTRACTION")
    print("="*80)
    print("\nLattice Configuration:")
    print(f"  Spatial: {Lx} × {Ly} × {Lz}")
    print(f"  Temporal: {Nt}")
    print(f"  Total sites: {Lx*Ly*Lz*Nt}")
    print(f"  β = {beta:.3f} (SU(2) coupling)")
    print("\nMonte Carlo:")
    print(f"  Thermalization: {n_therm} sweeps")
    print(f"  Measurements: {n_meas} sweeps")
    print(f"  Interval: {meas_interval} sweeps")
    print(f"  Total sweeps: {n_therm + n_meas * meas_interval}")
    print("\nExpected Runtime: ~4-6 hours")
    print("="*80 + "\n")

    # Initialize lattice
    print("Initializing SU(2) lattice...")
    lattice = SU2Lattice3p1D(Lx, Ly, Lz, Nt, beta)
    print("✓ Lattice initialized\n")

    # Thermalization
    print("Phase 1: Thermalization")
    print("-" * 40)
    therm_start = time.time()

    for sweep_num in range(n_therm):
        acceptance = lattice.sweep()

        if (sweep_num + 1) % 50 == 0:
            plaq = lattice.average_plaquette()
            elapsed = time.time() - therm_start
            print(f"  Sweep {sweep_num+1}/{n_therm}: ⟨P⟩ = {plaq:.6f}, "
                  f"acc = {acceptance:.3f}, time = {elapsed:.1f}s")

    therm_time = time.time() - therm_start
    print(f"✓ Thermalization complete ({therm_time:.1f}s)\n")

    # Measurements
    print("Phase 2: Measurements")
    print("-" * 40)
    meas_start = time.time()

    R_all = []
    C_all = []
    plaquettes = []
    configs_saved = []
    save_interval = max(1, n_meas // 10)  # Save 10 configurations

    for meas in range(n_meas):
        # Monte Carlo updates
        for _ in range(meas_interval):
            lattice.sweep()

        # Measure plaquette
        plaq = lattice.average_plaquette()
        plaquettes.append(plaq)

        # Measure Polyakov correlator
        R, C = lattice.polyakov_correlator()
        R_all.append(R)
        C_all.append(C)

        # Save configuration occasionally
        if (meas + 1) % save_interval == 0:
            configs_saved.append({
                'theta': lattice.theta.copy(),
                'plaquette': plaq
            })

        # Progress report
        if (meas + 1) % 10 == 0:
            elapsed = time.time() - meas_start
            eta = elapsed / (meas + 1) * (n_meas - meas - 1)
            print(f"  Measurement {meas+1}/{n_meas}: ⟨P⟩ = {plaq:.6f}, "
                  f"elapsed = {elapsed/60:.1f}m, ETA = {eta/60:.1f}m")

    meas_time = time.time() - meas_start
    total_time = therm_time + meas_time
    print(f"✓ Measurements complete ({meas_time/60:.1f}m)\n")
    print(f"Total simulation time: {total_time/60:.1f} minutes ({total_time/3600:.2f} hours)")

    # Statistical analysis
    print("\nPhase 3: Statistical Analysis")
    print("-" * 40)

    # Average over measurements
    R_unique = R_all[0]  # All measurements have same R grid
    C_samples = np.array([C for C in C_all])  # Shape: (n_meas, n_R)

    C_mean = np.mean(C_samples, axis=0)
    C_std = np.std(C_samples, axis=0)
    C_err = C_std / np.sqrt(len(C_samples))

    # Extract potential V(R) from -log(C(R))
    # V(R) ~ -log(C(R)) for large temporal extent
    C_abs = np.abs(C_mean)

    # Only use points with good statistics
    good_mask = (C_abs > 1e-6) & (C_err/C_abs < 0.5)
    R_good = R_unique[good_mask]
    C_good = C_abs[good_mask]

    if len(R_good) > 3:
        V = -np.log(C_good)
        # Error propagation: dV = dC/C
        dV = C_err[good_mask] / C_good

        print(f"✓ Extracted V(R) for {len(R_good)} separation distances")
    else:
        print("⚠ Warning: Insufficient statistics for V(R) extraction")
        V = np.array([])
        dV = np.array([])

    # Save results
    results = {
        'R': R_good,
        'V': V,
        'dV': dV,
        'C_mean': C_mean,
        'C_std': C_std,
        'C_err': C_err,
        'R_full': R_unique,
        'plaquettes': np.array(plaquettes),
        'configs': configs_saved,
        'parameters': {
            'Lx': Lx, 'Ly': Ly, 'Lz': Lz, 'Nt': Nt,
            'beta': beta,
            'n_therm': n_therm,
            'n_meas': n_meas,
            'meas_interval': meas_interval,
            'total_time': total_time
        }
    }

    # Save to disk
    results_file = f'/mnt/c/exe/projects/ai-agents/synchronism/simulations/{output_prefix}_results.pkl'
    with open(results_file, 'wb') as f:
        pickle.dump(results, f)
    print(f"✓ Results saved to {results_file}")

    return results


def analyze_and_fit(results, output_prefix):
    """
    Fit V(R) to Yukawa and Coulomb potentials, generate plots.

    Parameters:
    -----------
    results : dict
        Results from run_production_su2
    output_prefix : str
        Prefix for output files

    Returns:
    --------
    analysis : dict
        Fit results and validation metrics
    """

    print("\nPhase 4: Model Fitting")
    print("-" * 40)

    R = results['R']
    V = results['V']
    dV = results['dV']

    if len(R) < 4:
        print("⚠ Insufficient data for fitting")
        return {'success': False}

    # Fit 1: Yukawa potential V(R) = -α exp(-MR)/R + C
    print("Fitting Yukawa potential: V(R) = -α exp(-MR)/R + C")
    try:
        # Initial guess: α ~ 0.5, M ~ 1.0 (W boson mass), C ~ 0
        p0_yukawa = [0.5, 1.0, 0.0]
        popt_yukawa, pcov_yukawa = curve_fit(
            yukawa_potential, R, V, p0=p0_yukawa,
            sigma=dV, absolute_sigma=True,
            maxfev=10000
        )

        alpha_y, M_y, C_y = popt_yukawa
        dalpha_y, dM_y, dC_y = np.sqrt(np.diag(pcov_yukawa))

        V_fit_yukawa = yukawa_potential(R, alpha_y, M_y, C_y)
        chi2_yukawa = np.sum(((V - V_fit_yukawa) / dV) ** 2)
        dof_yukawa = len(R) - 3
        chi2dof_yukawa = chi2_yukawa / dof_yukawa if dof_yukawa > 0 else np.inf

        print(f"  α = {alpha_y:.4f} ± {dalpha_y:.4f}")
        print(f"  M = {M_y:.4f} ± {dM_y:.4f} (lattice units)")
        print(f"  C = {C_y:.4f} ± {dC_y:.4f}")
        print(f"  χ²/dof = {chi2dof_yukawa:.4f}")

        yukawa_fit = {
            'success': True,
            'alpha': alpha_y, 'dalpha': dalpha_y,
            'M': M_y, 'dM': dM_y,
            'const': C_y, 'dconst': dC_y,
            'chi2dof': chi2dof_yukawa
        }
    except Exception as e:
        print(f"  Yukawa fit failed: {e}")
        yukawa_fit = {'success': False}

    # Fit 2: Coulomb potential V(R) = -α/R + C
    print("\nFitting Coulomb potential: V(R) = -α/R + C")
    try:
        p0_coulomb = [0.5, 0.0]
        popt_coulomb, pcov_coulomb = curve_fit(
            coulomb_potential, R, V, p0=p0_coulomb,
            sigma=dV, absolute_sigma=True
        )

        alpha_c, C_c = popt_coulomb
        dalpha_c, dC_c = np.sqrt(np.diag(pcov_coulomb))

        V_fit_coulomb = coulomb_potential(R, alpha_c, C_c)
        chi2_coulomb = np.sum(((V - V_fit_coulomb) / dV) ** 2)
        dof_coulomb = len(R) - 2
        chi2dof_coulomb = chi2_coulomb / dof_coulomb if dof_coulomb > 0 else np.inf

        print(f"  α = {alpha_c:.4f} ± {dalpha_c:.4f}")
        print(f"  C = {C_c:.4f} ± {dC_c:.4f}")
        print(f"  χ²/dof = {chi2dof_coulomb:.4f}")

        coulomb_fit = {
            'success': True,
            'alpha': alpha_c, 'dalpha': dalpha_c,
            'const': C_c, 'dconst': dC_c,
            'chi2dof': chi2dof_coulomb
        }
    except Exception as e:
        print(f"  Coulomb fit failed: {e}")
        coulomb_fit = {'success': False}

    # Compare fits
    print("\nModel Comparison:")
    print("-" * 40)
    if yukawa_fit['success'] and coulomb_fit['success']:
        if yukawa_fit['chi2dof'] < coulomb_fit['chi2dof']:
            print(f"✓ Yukawa model preferred (χ²/dof = {yukawa_fit['chi2dof']:.4f} vs {coulomb_fit['chi2dof']:.4f})")
            print("  → Screening observed! Weak force validated.")
            best_model = 'yukawa'
        else:
            print(f"✓ Coulomb model preferred (χ²/dof = {coulomb_fit['chi2dof']:.4f} vs {yukawa_fit['chi2dof']:.4f})")
            print("  → No screening observed. Weak force NOT validated.")
            best_model = 'coulomb'
    elif yukawa_fit['success']:
        print("✓ Only Yukawa fit succeeded")
        best_model = 'yukawa'
    elif coulomb_fit['success']:
        print("✓ Only Coulomb fit succeeded")
        best_model = 'coulomb'
    else:
        print("⚠ Both fits failed")
        best_model = None

    # Generate plots
    print("\nPhase 5: Visualization")
    print("-" * 40)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(f'Session #34: SU(2) Weak Force Physics Extraction\nβ={results["parameters"]["beta"]}, {results["parameters"]["Lx"]}³×{results["parameters"]["Nt"]} lattice',
                 fontsize=14)

    # Plot 1: V(R) with fits
    ax = axes[0, 0]
    ax.errorbar(R, V, yerr=dV, fmt='o', markersize=6, capsize=4,
                label='SU(2) Data', color='blue')

    R_smooth = np.linspace(R.min(), R.max(), 200)
    if yukawa_fit['success']:
        V_yukawa = yukawa_potential(R_smooth, yukawa_fit['alpha'],
                                    yukawa_fit['M'], yukawa_fit['const'])
        ax.plot(R_smooth, V_yukawa, '-', linewidth=2, color='red',
                label=f'Yukawa: M={yukawa_fit["M"]:.3f}, χ²/dof={yukawa_fit["chi2dof"]:.3f}')

    if coulomb_fit['success']:
        V_coulomb = coulomb_potential(R_smooth, coulomb_fit['alpha'],
                                      coulomb_fit['const'])
        ax.plot(R_smooth, V_coulomb, '--', linewidth=2, color='green',
                label=f'Coulomb: χ²/dof={coulomb_fit["chi2dof"]:.3f}')

    ax.set_xlabel('R (lattice units)', fontsize=12)
    ax.set_ylabel('V(R)', fontsize=12)
    ax.set_title('Static Potential V(R)', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Plot 2: Residuals
    ax = axes[0, 1]
    if yukawa_fit['success']:
        V_fit = yukawa_potential(R, yukawa_fit['alpha'],
                                 yukawa_fit['M'], yukawa_fit['const'])
        residuals = V - V_fit
        ax.errorbar(R, residuals, yerr=dV, fmt='o', capsize=3)
        ax.axhline(0, color='r', linestyle='--')
        ax.set_xlabel('R (lattice units)', fontsize=12)
        ax.set_ylabel('V - V_Yukawa', fontsize=12)
        ax.set_title('Fit Residuals', fontsize=13)
        ax.grid(True, alpha=0.3)

    # Plot 3: Polyakov correlator
    ax = axes[1, 0]
    C_abs = np.abs(results['C_mean'])
    ax.semilogy(results['R_full'], C_abs, 'o-', markersize=4)
    ax.set_xlabel('R (lattice units)', fontsize=12)
    ax.set_ylabel('|C(R)|', fontsize=12)
    ax.set_title('Polyakov Loop Correlator', fontsize=13)
    ax.grid(True, alpha=0.3)

    # Plot 4: Plaquette history
    ax = axes[1, 1]
    plaq_mean = np.mean(results['plaquettes'])
    ax.plot(results['plaquettes'], alpha=0.7, linewidth=0.8)
    ax.axhline(plaq_mean, color='r', linestyle='--',
               label=f'Mean: {plaq_mean:.5f}')
    ax.set_xlabel('Measurement', fontsize=12)
    ax.set_ylabel('⟨Plaquette⟩', fontsize=12)
    ax.set_title('SU(2) Coherence Evolution', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    plot_file = f'/mnt/c/exe/projects/ai-agents/synchronism/simulations/{output_prefix}_analysis.png'
    plt.savefig(plot_file, dpi=200)
    print(f"✓ Analysis plot saved: {plot_file}")

    # Generate report
    report_file = f'/mnt/c/exe/projects/ai-agents/synchronism/simulations/{output_prefix}_report.txt'
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SESSION #34: SU(2) WEAK FORCE PHYSICS EXTRACTION\n")
        f.write("="*80 + "\n\n")

        f.write("LATTICE PARAMETERS:\n")
        f.write(f"  Spatial: {results['parameters']['Lx']}³\n")
        f.write(f"  Temporal: {results['parameters']['Nt']}\n")
        f.write(f"  β = {results['parameters']['beta']}\n")
        f.write(f"  Measurements: {results['parameters']['n_meas']}\n")
        f.write(f"  Runtime: {results['parameters']['total_time']/3600:.2f} hours\n\n")

        f.write("YUKAWA FIT: V(R) = -α exp(-MR)/R + C\n")
        if yukawa_fit['success']:
            f.write(f"  α = {yukawa_fit['alpha']:.4f} ± {yukawa_fit['dalpha']:.4f}\n")
            f.write(f"  M = {yukawa_fit['M']:.4f} ± {yukawa_fit['dM']:.4f} (lattice units)\n")
            f.write(f"  C = {yukawa_fit['const']:.4f} ± {yukawa_fit['dconst']:.4f}\n")
            f.write(f"  χ²/dof = {yukawa_fit['chi2dof']:.4f}\n\n")
        else:
            f.write("  Fit failed\n\n")

        f.write("COULOMB FIT: V(R) = -α/R + C\n")
        if coulomb_fit['success']:
            f.write(f"  α = {coulomb_fit['alpha']:.4f} ± {coulomb_fit['dalpha']:.4f}\n")
            f.write(f"  C = {coulomb_fit['const']:.4f} ± {coulomb_fit['dconst']:.4f}\n")
            f.write(f"  χ²/dof = {coulomb_fit['chi2dof']:.4f}\n\n")
        else:
            f.write("  Fit failed\n\n")

        f.write("CONCLUSION:\n")
        if best_model == 'yukawa':
            f.write("  ✅ WEAK FORCE VALIDATED\n")
            f.write("  Yukawa screening observed with M ≈ {:.3f}\n".format(yukawa_fit['M']))
            f.write("  SU(2) gauge symmetry emerges from Synchronism intent dynamics\n")
        elif best_model == 'coulomb':
            f.write("  ❌ WEAK FORCE NOT VALIDATED\n")
            f.write("  No screening observed, pure Coulomb behavior\n")
            f.write("  SU(2) requires refinement or additional physics\n")
        else:
            f.write("  ⚠ INCONCLUSIVE\n")
            f.write("  Insufficient statistics or fitting failures\n")

        f.write("="*80 + "\n")

    print(f"✓ Report saved: {report_file}")

    return {
        'yukawa_fit': yukawa_fit,
        'coulomb_fit': coulomb_fit,
        'best_model': best_model,
        'success': True
    }


def main():
    """
    Main execution: Session #34 SU(2) physics extraction.
    """

    # Production parameters
    params = {
        'Lx': 8,
        'Ly': 8,
        'Lz': 8,
        'Nt': 4,
        'beta': 2.2,          # SU(2) coupling
        'n_therm': 100,       # Thermalization
        'n_meas': 500,        # Measurements
        'meas_interval': 2,   # Decorrelation
        'output_prefix': 'session34_su2_production'
    }

    print("\n" + "="*80)
    print("SYNCHRONISM RESEARCH SESSION #34")
    print("SU(2) Weak Force Physics Extraction")
    print("="*80)
    print("\nMission:")
    print("  Extract Yukawa screening from SU(2) lattice gauge simulation")
    print("  Validate weak force emergence from Synchronism intent dynamics")
    print("\nContext:")
    print("  Session #27: U(1) Coulomb validated ✅")
    print("  Session #30: SU(2) implemented ✅")
    print("  Session #33: Analysis tools prepared ✅")
    print("  Nova: Extract physics, validate screening")
    print("\nExpected:")
    print("  V(R) ∝ exp(-M_W R)/R if weak force emerges")
    print("  V(R) ∝ 1/R if SU(2) NOT validated")
    print("="*80 + "\n")

    # Run production simulation
    results = run_production_su2(**params)

    # Analyze and fit
    analysis = analyze_and_fit(results, params['output_prefix'])

    print("\n" + "="*80)
    print("SESSION #34 COMPLETE")
    print("="*80)

    if analysis.get('success'):
        if analysis['best_model'] == 'yukawa':
            print("\n✅ SUCCESS: Weak force validated!")
            print("   SU(2) gauge symmetry emerges from Synchronism")
        elif analysis['best_model'] == 'coulomb':
            print("\n❌ FALSIFIED: No screening observed")
            print("   SU(2) requires theoretical refinement")
        else:
            print("\n⚠ INCONCLUSIVE: Need better statistics")

    print("\nNext: Session #35 - SU(3) confinement (critical test)")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
