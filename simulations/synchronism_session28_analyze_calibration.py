"""
Coupling Calibration Analysis - Session #28
Synchronism Research: Determining β_target for α = 1/137.036

Analyzes all β calibration runs to establish the mapping:
    Synchronism Intent Coherence (β) → Physical QED Coupling (α)

Determines β_target where α_eff = 1/137.036 ≈ 0.00730

Usage:
    python3 synchronism_session28_analyze_calibration.py

Requires: All session28_calibration_beta_*.npz files from calibration runs

Author: CBP Autonomous Synchronism Research
Date: 2025-11-19
Session: #28
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
import os

def load_calibration_data():
    """
    Load all calibration data points from individual β runs.

    Returns:
    --------
    data : dict
        Contains beta, alpha_eff, alpha_err, chi2_dof arrays
    """
    # Find all calibration files
    pattern = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/session28_calibration_beta_*.npz'
    files = sorted(glob.glob(pattern))

    if len(files) == 0:
        print(f"Error: No calibration files found matching {pattern}")
        print("Run synchronism_session28_coupling_calibration.py for multiple β values first")
        return None

    print(f"Found {len(files)} calibration data points")

    # Load data
    beta_values = []
    alpha_values = []
    alpha_errors = []
    chi2_values = []

    for f in files:
        data = np.load(f)
        beta_values.append(data['beta'])
        alpha_values.append(data['alpha_eff'])
        alpha_errors.append(data['alpha_err'])
        chi2_values.append(data['chi2_dof'])

        print(f"  β = {data['beta']:.4f}: α = {data['alpha_eff']:.6f} ± {data['alpha_err']:.6f}, χ²/dof = {data['chi2_dof']:.2f}")

    return {
        'beta': np.array(beta_values),
        'alpha_eff': np.array(alpha_values),
        'alpha_err': np.array(alpha_errors),
        'chi2_dof': np.array(chi2_values)
    }


def fit_scaling_relation(data, model='power'):
    """
    Fit α_eff(β) to theoretical scaling form.

    Models:
    -------
    'power': α_eff = A / β^B
    'log': α_eff = A + B * log(β)
    'linear': α_eff = A / β

    Returns:
    --------
    fit_params : dict
        Fit parameters and statistics
    """
    beta = data['beta']
    alpha = data['alpha_eff']
    alpha_err = data['alpha_err']

    if model == 'power':
        # α_eff = A / β^B
        def func(b, A, B):
            return A / b**B

        # Initial guess: A ~ 0.1, B ~ 1
        p0 = [0.1, 1.0]

        popt, pcov = curve_fit(func, beta, alpha, p0=p0, sigma=alpha_err, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))

        return {
            'model': 'power',
            'func': func,
            'A': popt[0],
            'A_err': perr[0],
            'B': popt[1],
            'B_err': perr[1],
            'params': popt,
            'formula': f'α = {popt[0]:.4f}/β^{popt[1]:.3f}'
        }

    elif model == 'linear':
        # α_eff = A / β
        def func(b, A):
            return A / b

        p0 = [0.1]

        popt, pcov = curve_fit(func, beta, alpha, p0=p0, sigma=alpha_err, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))

        return {
            'model': 'linear',
            'func': func,
            'A': popt[0],
            'A_err': perr[0],
            'params': popt,
            'formula': f'α = {popt[0]:.4f}/β'
        }

    elif model == 'log':
        # α_eff = A + B * log(β)
        def func(b, A, B):
            return A + B * np.log(b)

        p0 = [0.1, -0.05]

        popt, pcov = curve_fit(func, beta, alpha, p0=p0, sigma=alpha_err, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))

        return {
            'model': 'log',
            'func': func,
            'A': popt[0],
            'A_err': perr[0],
            'B': popt[1],
            'B_err': perr[1],
            'params': popt,
            'formula': f'α = {popt[0]:.4f} + {popt[1]:.4f}*log(β)'
        }


def find_beta_target(data, fit_params, alpha_target=1/137.036):
    """
    Determine β value that produces α_eff = alpha_target.

    Uses fitted scaling relation to interpolate or extrapolate.

    Returns:
    --------
    beta_target : float
        β value for target α
    beta_err : float
        Estimated uncertainty (rough)
    """
    func = fit_params['func']
    params = fit_params['params']

    # Solve numerically: func(β, *params) = alpha_target
    from scipy.optimize import fsolve

    # Initial guess based on data range
    if alpha_target < data['alpha_eff'].min():
        # Extrapolate to higher β
        beta_guess = data['beta'].max() * 2
    elif alpha_target > data['alpha_eff'].max():
        # Extrapolate to lower β
        beta_guess = data['beta'].min() / 2
    else:
        # Interpolate
        beta_guess = np.interp(alpha_target, data['alpha_eff'][::-1], data['beta'][::-1])

    # Solve
    def equation(b):
        return func(b, *params) - alpha_target

    beta_target = fsolve(equation, beta_guess)[0]

    # Estimate uncertainty (propagate fit parameter errors)
    # Rough estimate: δβ ~ β * (δA/A) for power law
    if 'A_err' in fit_params:
        beta_err = beta_target * (fit_params['A_err'] / fit_params['A'])
    else:
        beta_err = beta_target * 0.2  # 20% uncertainty if can't estimate

    return beta_target, beta_err


def plot_calibration(data, fits, beta_target, alpha_target=1/137.036):
    """
    Create comprehensive calibration plots.

    Panels:
    1. α_eff vs β with fits
    2. 1/α_eff vs β (check linearity)
    3. Residuals
    4. χ²/dof quality check
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel 1: α_eff(β) with fits
    ax = axes[0, 0]

    ax.errorbar(data['beta'], data['alpha_eff'], yerr=data['alpha_err'],
                fmt='o', label='Measured', capsize=3, markersize=8)

    # Plot fits
    beta_smooth = np.linspace(data['beta'].min()*0.5, data['beta'].max()*2, 200)

    colors = ['red', 'green', 'blue']
    for i, (name, fit) in enumerate(fits.items()):
        alpha_fit = fit['func'](beta_smooth, *fit['params'])
        ax.plot(beta_smooth, alpha_fit, '--', color=colors[i], label=f"{name}: {fit['formula']}", alpha=0.7)

    # Mark target
    ax.axhline(alpha_target, color='orange', linestyle=':', linewidth=2, label=f'Target: α = 1/137 = {alpha_target:.6f}')
    ax.axvline(beta_target, color='purple', linestyle=':', linewidth=2, label=f'β_target = {beta_target:.2f}')

    ax.set_xlabel('β (coherence coupling)', fontsize=12)
    ax.set_ylabel('α_eff (effective fine-structure)', fontsize=12)
    ax.set_title('Coupling Calibration Curve', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Panel 2: 1/α vs β (check if linear)
    ax = axes[0, 1]

    alpha_inv = 1.0 / data['alpha_eff']
    alpha_inv_err = data['alpha_err'] / data['alpha_eff']**2

    ax.errorbar(data['beta'], alpha_inv, yerr=alpha_inv_err,
                fmt='o', capsize=3, markersize=8)
    ax.axhline(137.036, color='orange', linestyle=':', linewidth=2, label='Target: 1/α = 137')

    ax.set_xlabel('β (coherence coupling)', fontsize=12)
    ax.set_ylabel('1/α_eff', fontsize=12)
    ax.set_title('Inverse Coupling (Check Linearity)', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 3: Fit residuals (best fit)
    ax = axes[1, 0]

    # Use power law fit
    best_fit = fits['power']
    alpha_fit_data = best_fit['func'](data['beta'], *best_fit['params'])
    residuals = (data['alpha_eff'] - alpha_fit_data) / data['alpha_err']

    ax.errorbar(data['beta'], residuals, yerr=np.ones_like(residuals),
                fmt='o', capsize=3, markersize=8)
    ax.axhline(0, color='red', linestyle='--')
    ax.axhline(2, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(-2, color='gray', linestyle=':', alpha=0.5)

    ax.set_xlabel('β (coherence coupling)', fontsize=12)
    ax.set_ylabel('Residuals (σ)', fontsize=12)
    ax.set_title('Fit Quality (Power Law Model)', fontsize=13)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')

    # Panel 4: χ²/dof check
    ax = axes[1, 1]

    ax.plot(data['beta'], data['chi2_dof'], 'o-', markersize=8)
    ax.axhline(1, color='green', linestyle='--', label='Perfect fit')
    ax.axhline(2, color='orange', linestyle=':', label='Acceptable')

    ax.set_xlabel('β (coherence coupling)', fontsize=12)
    ax.set_ylabel('χ²/dof', fontsize=12)
    ax.set_title('V(R) Fit Quality at Each β', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')

    plt.tight_layout()

    plot_file = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/session28_calibration_analysis.png'
    plt.savefig(plot_file, dpi=150)
    print(f"\n✓ Saved plot: {plot_file}")


def main():
    """
    Main analysis: Load calibration data, fit scaling relation, determine β_target.
    """
    print("="*70)
    print("SESSION #28: COUPLING CALIBRATION ANALYSIS")
    print("="*70)
    print()

    # Load data
    data = load_calibration_data()
    if data is None:
        return

    print("\n" + "="*70)
    print("FITTING SCALING RELATIONS")
    print("="*70)

    # Try multiple models
    fits = {}

    for model in ['power', 'linear', 'log']:
        try:
            fit = fit_scaling_relation(data, model=model)
            fits[model] = fit
            print(f"\n{model.capitalize()} law fit:")
            print(f"  {fit['formula']}")

            if 'A_err' in fit:
                print(f"  A = {fit['A']:.6f} ± {fit['A_err']:.6f}")
            if 'B_err' in fit:
                print(f"  B = {fit['B']:.3f} ± {fit['B_err']:.3f}")
        except Exception as e:
            print(f"\n{model.capitalize()} fit failed: {e}")

    # Determine β_target (use power law as primary)
    print("\n" + "="*70)
    print("DETERMINING β_TARGET FOR α = 1/137.036")
    print("="*70)

    alpha_target = 1.0 / 137.036
    print(f"\nTarget: α = {alpha_target:.8f}")

    if 'power' in fits:
        beta_target, beta_err = find_beta_target(data, fits['power'], alpha_target)
        print(f"\nUsing power law model:")
        print(f"  β_target = {beta_target:.2f} ± {beta_err:.2f}")

        # Check if extrapolation needed
        if beta_target < data['beta'].min():
            print(f"  ⚠ WARNING: Extrapolating to LOWER β (outside measured range)")
        elif beta_target > data['beta'].max():
            print(f"  ⚠ WARNING: Extrapolating to HIGHER β (outside measured range)")
        else:
            print(f"  ✓ Interpolating within measured range")

        # Extrapolation uncertainty
        extrap_factor = max(beta_target / data['beta'].max(), data['beta'].min() / beta_target)
        if extrap_factor > 1.2:
            print(f"  ⚠ Extrapolation factor: {extrap_factor:.2f}x (uncertainty may be larger)")

    else:
        print("  ⚠ Power law fit failed, cannot determine β_target")
        beta_target = None

    # Create plots
    print("\n" + "="*70)
    print("GENERATING CALIBRATION PLOTS")
    print("="*70)

    if beta_target is not None:
        plot_calibration(data, fits, beta_target, alpha_target)

    # Summary
    print("\n" + "="*70)
    print("SESSION #28 CALIBRATION SUMMARY")
    print("="*70)
    print(f"\nData points: {len(data['beta'])}")
    print(f"β range: {data['beta'].min():.2f} - {data['beta'].max():.2f}")
    print(f"α_eff range: {data['alpha_eff'].min():.6f} - {data['alpha_eff'].max():.6f}")

    if beta_target is not None:
        print(f"\n🎯 RESULT: β_target = {beta_target:.2f} ± {beta_err:.2f}")
        print(f"   Produces: α_eff ≈ {alpha_target:.6f} (1/137.036)")
        print(f"\n   Synchronism Interpretation:")
        print(f"   Intent coherence strength β = {beta_target:.2f}")
        print(f"   generates observed electromagnetic coupling α = 1/137")

        if extrap_factor < 1.2:
            print(f"\n   ✓ High confidence (interpolated within data)")
        else:
            print(f"\n   ⚠ Moderate confidence (extrapolated {extrap_factor:.2f}x)")
            print(f"   Recommend: Run simulation at β = {beta_target:.2f} to validate")

    print("\n" + "="*70)
    print("DONE!")
    print("="*70)


if __name__ == '__main__':
    main()
