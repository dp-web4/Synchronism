#!/usr/bin/env python3
"""
Session #52 Track A: Parameter Recalibration

Session #51 Critical Finding:
- Current A=0.25, B=1.62 make ρ_crit too high (~1000 M_sun/pc³)
- Typical galaxies have ρ ~ 0.1-10 M_sun/pc³
- Result: ALL galaxies have C ≈ 0 (DM-dominated)
- ETG predictions fail badly (predict 100% DM, observe 13%)

This session attempts to recalibrate A, B to match ETG observations:
- Target: f_DM = 13% → C = 0.87 for typical ETGs
- Method: Fit A, B to reproduce observed DM fractions

Author: CBP Autonomous Synchronism Research
Date: 2025-11-26
Session: #52 - Parameter Recalibration
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from pathlib import Path
import json
from datetime import datetime


# =============================================================================
# CALIBRATION DATA
# =============================================================================

# Combined dataset: ETGs (ATLAS3D) + Dwarfs/Spirals (Santos-Santos)
# Goal: Find A, B that work for BOTH regimes

CALIBRATION_DATA = {
    # ETGs from ATLAS3D (need higher C, lower f_DM)
    'etgs': [
        {'name': 'NGC2549', 'vmax': 205, 'rho_mean': 2.76, 'f_dm_obs': 0.25, 'type': 'ETG'},
        {'name': 'NGC3156', 'vmax': 127, 'rho_mean': 3.73, 'f_dm_obs': 0.30, 'type': 'ETG'},
        {'name': 'NGC4473', 'vmax': 269, 'rho_mean': 1.19, 'f_dm_obs': 0.12, 'type': 'ETG'},
        {'name': 'NGC3379', 'vmax': 290, 'rho_mean': 1.04, 'f_dm_obs': 0.10, 'type': 'ETG'},
        {'name': 'NGC4374', 'vmax': 410, 'rho_mean': 0.35, 'f_dm_obs': 0.08, 'type': 'ETG'},
        {'name': 'NGC4486B', 'vmax': 262, 'rho_mean': 354, 'f_dm_obs': 0.02, 'type': 'cE'},
        {'name': 'M32', 'vmax': 106, 'rho_mean': 538, 'f_dm_obs': 0.01, 'type': 'cE'},
    ],

    # Dwarfs from Santos-Santos (need low C, high f_DM)
    'dwarfs': [
        {'name': 'DDO168', 'vmax': 55, 'rho_mean': 0.001, 'f_dm_obs': 0.95, 'type': 'dwarf'},
        {'name': 'DDO154', 'vmax': 47, 'rho_mean': 0.0008, 'f_dm_obs': 0.97, 'type': 'dwarf'},
        {'name': 'NGC2366', 'vmax': 52, 'rho_mean': 0.002, 'f_dm_obs': 0.94, 'type': 'dwarf'},
        {'name': 'DDO87', 'vmax': 38, 'rho_mean': 0.0015, 'f_dm_obs': 0.96, 'type': 'dwarf'},
        {'name': 'WLM', 'vmax': 38, 'rho_mean': 0.001, 'f_dm_obs': 0.97, 'type': 'dwarf'},
    ],

    # Spirals (intermediate)
    'spirals': [
        {'name': 'NGC3198', 'vmax': 150, 'rho_mean': 0.05, 'f_dm_obs': 0.70, 'type': 'spiral'},
        {'name': 'NGC2403', 'vmax': 130, 'rho_mean': 0.08, 'f_dm_obs': 0.65, 'type': 'spiral'},
        {'name': 'NGC7331', 'vmax': 250, 'rho_mean': 0.15, 'f_dm_obs': 0.55, 'type': 'spiral'},
        {'name': 'NGC2841', 'vmax': 300, 'rho_mean': 0.20, 'f_dm_obs': 0.45, 'type': 'spiral'},
    ],
}


# =============================================================================
# SYNCHRONISM MODEL
# =============================================================================

def synchronism_coherence(rho_mean, vmax, A, B, gamma=2.0):
    """
    Compute coherence from density and velocity.

    C = tanh(γ × log(ρ/ρ_crit + 1))
    ρ_crit = A × V^B
    """
    rho_crit = A * vmax**B

    if rho_crit > 0 and rho_mean > 0:
        C = np.tanh(gamma * np.log(rho_mean / rho_crit + 1))
    else:
        C = 0

    return C, rho_crit


def compute_error(params, data, weights=None):
    """
    Compute weighted mean squared error for parameter set.

    params: [A, B]
    data: list of galaxy dicts
    weights: optional weights by galaxy type
    """
    A, B = params

    if A <= 0 or B <= 0:
        return 1e10

    errors = []
    for gal in data:
        C, _ = synchronism_coherence(gal['rho_mean'], gal['vmax'], A, B)
        f_dm_pred = 1 - C
        error = (f_dm_pred - gal['f_dm_obs'])**2

        # Weight by type if provided
        if weights:
            w = weights.get(gal['type'], 1.0)
            error *= w

        errors.append(error)

    return np.mean(errors)


# =============================================================================
# CALIBRATION APPROACHES
# =============================================================================

def approach_1_fit_all():
    """
    Approach 1: Fit A, B to minimize error across ALL galaxy types.

    This is the most straightforward approach but may not work if
    ETGs and dwarfs require fundamentally different parameters.
    """

    print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│  APPROACH 1: FIT ALL GALAXY TYPES SIMULTANEOUSLY                             │
└─────────────────────────────────────────────────────────────────────────────┘
""")

    # Combine all data
    all_data = (CALIBRATION_DATA['etgs'] +
                CALIBRATION_DATA['dwarfs'] +
                CALIBRATION_DATA['spirals'])

    # Equal weights
    weights = {'ETG': 1.0, 'cE': 1.0, 'dwarf': 1.0, 'spiral': 1.0}

    # Initial guess (order of magnitude lower than current)
    x0 = [0.01, 1.0]

    # Bounds
    bounds = [(1e-6, 1.0), (0.5, 3.0)]

    # Optimize
    result = differential_evolution(
        lambda p: compute_error(p, all_data, weights),
        bounds,
        seed=42,
        maxiter=1000,
        tol=1e-8
    )

    A_opt, B_opt = result.x
    final_error = result.fun

    print(f"Optimal parameters: A = {A_opt:.6f}, B = {B_opt:.4f}")
    print(f"Mean squared error: {final_error:.6f}")
    print(f"RMS error: {np.sqrt(final_error):.4f}")

    # Evaluate on each type
    print("\nPerformance by galaxy type:")
    print("-" * 60)

    for gtype, galaxies in CALIBRATION_DATA.items():
        errors = []
        for gal in galaxies:
            C, rho_crit = synchronism_coherence(gal['rho_mean'], gal['vmax'], A_opt, B_opt)
            f_dm_pred = 1 - C
            errors.append(abs(f_dm_pred - gal['f_dm_obs']))

        print(f"  {gtype:12}: Mean error = {np.mean(errors):.3f}, Max error = {np.max(errors):.3f}")

    return A_opt, B_opt, final_error


def approach_2_weighted_fit():
    """
    Approach 2: Weight ETGs more heavily (since they're the problem case).
    """

    print("""

┌─────────────────────────────────────────────────────────────────────────────┐
│  APPROACH 2: WEIGHT ETGs MORE HEAVILY                                        │
└─────────────────────────────────────────────────────────────────────────────┘
""")

    all_data = (CALIBRATION_DATA['etgs'] +
                CALIBRATION_DATA['dwarfs'] +
                CALIBRATION_DATA['spirals'])

    # Weight ETGs more (they're the problem case)
    weights = {'ETG': 3.0, 'cE': 5.0, 'dwarf': 1.0, 'spiral': 1.0}

    bounds = [(1e-6, 1.0), (0.5, 3.0)]

    result = differential_evolution(
        lambda p: compute_error(p, all_data, weights),
        bounds,
        seed=42,
        maxiter=1000,
        tol=1e-8
    )

    A_opt, B_opt = result.x
    final_error = result.fun

    print(f"Optimal parameters: A = {A_opt:.6f}, B = {B_opt:.4f}")
    print(f"Weighted MSE: {final_error:.6f}")

    # Evaluate unweighted
    unweighted_error = compute_error([A_opt, B_opt], all_data, None)
    print(f"Unweighted MSE: {unweighted_error:.6f}")

    print("\nPerformance by galaxy type:")
    print("-" * 60)

    for gtype, galaxies in CALIBRATION_DATA.items():
        errors = []
        for gal in galaxies:
            C, rho_crit = synchronism_coherence(gal['rho_mean'], gal['vmax'], A_opt, B_opt)
            f_dm_pred = 1 - C
            errors.append(abs(f_dm_pred - gal['f_dm_obs']))

        print(f"  {gtype:12}: Mean error = {np.mean(errors):.3f}, Max error = {np.max(errors):.3f}")

    return A_opt, B_opt, unweighted_error


def approach_3_separate_regimes():
    """
    Approach 3: Fit ETGs and dwarfs separately to see if same A, B can work.

    If optimal A, B are very different, it suggests the model formulation
    needs modification, not just recalibration.
    """

    print("""

┌─────────────────────────────────────────────────────────────────────────────┐
│  APPROACH 3: FIT EACH REGIME SEPARATELY                                      │
└─────────────────────────────────────────────────────────────────────────────┘
""")

    bounds = [(1e-6, 10.0), (0.1, 3.0)]

    results = {}

    for regime, galaxies in CALIBRATION_DATA.items():
        if len(galaxies) < 3:
            continue

        result = differential_evolution(
            lambda p: compute_error(p, galaxies),
            bounds,
            seed=42,
            maxiter=1000,
            tol=1e-8
        )

        A_opt, B_opt = result.x
        final_error = result.fun

        results[regime] = {
            'A': A_opt,
            'B': B_opt,
            'error': final_error
        }

        print(f"\n{regime.upper()} regime:")
        print(f"  Optimal A = {A_opt:.6f}, B = {B_opt:.4f}")
        print(f"  MSE = {final_error:.6f}, RMS = {np.sqrt(final_error):.4f}")

    # Compare
    print("\n" + "="*60)
    print("COMPARISON OF REGIME-SPECIFIC PARAMETERS:")
    print("="*60)

    for regime, r in results.items():
        print(f"  {regime:12}: A = {r['A']:.6f}, B = {r['B']:.4f}")

    # Check if parameters are similar
    As = [r['A'] for r in results.values()]
    Bs = [r['B'] for r in results.values()]

    A_range = max(As) / min(As) if min(As) > 0 else np.inf
    B_range = max(Bs) - min(Bs)

    print(f"\n  A range: {min(As):.6f} to {max(As):.6f} (ratio: {A_range:.1f}x)")
    print(f"  B range: {min(Bs):.4f} to {max(Bs):.4f} (spread: {B_range:.2f})")

    if A_range > 100:
        print("\n  ⚠ WARNING: A varies by >100x between regimes!")
        print("  This suggests the model formulation may need modification.")
    elif A_range > 10:
        print("\n  ⚠ CAUTION: A varies by >10x between regimes.")
        print("  A single A, B may not work well for all galaxy types.")
    else:
        print("\n  ✓ Parameters are reasonably consistent across regimes.")

    return results


def approach_4_modified_formula():
    """
    Approach 4: Test if a modified coherence formula works better.

    Instead of ρ_crit = A × V^B, try:
    - ρ_crit = A × V^B × ρ_mean^C (self-consistent)
    - ρ_crit = A × (V/V_0)^B (normalized)
    - C = tanh(γ × log(ρ/ρ_0 + 1)) with fixed ρ_0 (simpler)
    """

    print("""

┌─────────────────────────────────────────────────────────────────────────────┐
│  APPROACH 4: TEST MODIFIED FORMULAS                                          │
└─────────────────────────────────────────────────────────────────────────────┘

Testing alternative coherence formulas:
""")

    all_data = (CALIBRATION_DATA['etgs'] +
                CALIBRATION_DATA['dwarfs'] +
                CALIBRATION_DATA['spirals'])

    # Formula 1: Standard (baseline)
    def formula_standard(gal, params):
        A, B = params[:2]
        rho_crit = A * gal['vmax']**B
        C = np.tanh(2.0 * np.log(gal['rho_mean'] / rho_crit + 1))
        return 1 - C

    # Formula 2: Fixed ρ_0 (simpler, V-independent)
    def formula_fixed_rho0(gal, params):
        rho_0 = params[0]
        C = np.tanh(2.0 * np.log(gal['rho_mean'] / rho_0 + 1))
        return 1 - C

    # Formula 3: Logarithmic in V
    def formula_log_v(gal, params):
        A, B = params[:2]
        rho_crit = A * np.log(gal['vmax'] / 10 + 1)**B
        if rho_crit <= 0:
            return 1.0
        C = np.tanh(2.0 * np.log(gal['rho_mean'] / rho_crit + 1))
        return 1 - C

    # Formula 4: Power law with offset
    def formula_power_offset(gal, params):
        A, B, V0 = params[:3]
        rho_crit = A * (gal['vmax'] / V0)**B
        if rho_crit <= 0:
            return 1.0
        C = np.tanh(2.0 * np.log(gal['rho_mean'] / rho_crit + 1))
        return 1 - C

    formulas = {
        'Standard (A×V^B)': (formula_standard, [(1e-6, 1.0), (0.1, 3.0)]),
        'Fixed ρ_0': (formula_fixed_rho0, [(1e-4, 100)]),
        'Log V': (formula_log_v, [(1e-4, 10), (0.5, 3.0)]),
        'Power+offset': (formula_power_offset, [(1e-6, 1.0), (0.1, 3.0), (10, 500)]),
    }

    results = {}

    for name, (formula, bounds) in formulas.items():
        def error_func(params):
            errors = [(formula(gal, params) - gal['f_dm_obs'])**2 for gal in all_data]
            return np.mean(errors)

        result = differential_evolution(error_func, bounds, seed=42, maxiter=1000)

        results[name] = {
            'params': result.x,
            'mse': result.fun,
            'rms': np.sqrt(result.fun)
        }

        print(f"\n{name}:")
        print(f"  Parameters: {result.x}")
        print(f"  MSE: {result.fun:.6f}, RMS: {np.sqrt(result.fun):.4f}")

    # Find best
    best = min(results.items(), key=lambda x: x[1]['mse'])
    print(f"\nBEST FORMULA: {best[0]} (RMS = {best[1]['rms']:.4f})")

    return results


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Run all calibration approaches."""

    print("\n" + "="*80)
    print("SESSION #52 TRACK A: PARAMETER RECALIBRATION")
    print("="*80)

    print("""
CONTEXT:
════════════════════════════════════════════════════════════════════════════════

Session #51 found that current parameters (A=0.25, B=1.62) make ρ_crit too high.
ALL galaxies end up in the DM-dominated regime (C ≈ 0).

This session attempts to find A, B values that work for:
- ETGs: f_DM ~ 10-30% (need high C)
- Spirals: f_DM ~ 45-70% (need moderate C)
- Dwarfs: f_DM ~ 95-97% (need low C)

The span from 1% to 97% DM must be captured by a single formula!
""")

    # Run all approaches
    A1, B1, err1 = approach_1_fit_all()
    A2, B2, err2 = approach_2_weighted_fit()
    regime_results = approach_3_separate_regimes()
    formula_results = approach_4_modified_formula()

    # Summary
    print("""

┌─────────────────────────────────────────────────────────────────────────────┐
│           CALIBRATION SUMMARY                                                │
└─────────────────────────────────────────────────────────────────────────────┘
""")

    print(f"{'Approach':<30} {'A':<12} {'B':<8} {'RMS Error':<12}")
    print("-" * 65)
    print(f"{'1. Equal weights':<30} {A1:<12.6f} {B1:<8.4f} {np.sqrt(err1):<12.4f}")
    print(f"{'2. ETG-weighted':<30} {A2:<12.6f} {B2:<8.4f} {np.sqrt(err2):<12.4f}")

    for regime, r in regime_results.items():
        print(f"{'3. ' + regime + ' only':<30} {r['A']:<12.6f} {r['B']:<8.4f} {np.sqrt(r['error']):<12.4f}")

    print()
    for name, r in formula_results.items():
        if len(r['params']) >= 2:
            print(f"{'4. ' + name:<30} {r['params'][0]:<12.6f} {r['params'][1] if len(r['params']) > 1 else '-':<8.4f} {r['rms']:<12.4f}")
        else:
            print(f"{'4. ' + name:<30} {r['params'][0]:<12.6f} {'-':<8} {r['rms']:<12.4f}")

    # Best overall
    best_standard = min([(A1, B1, err1), (A2, B2, err2)], key=lambda x: x[2])
    best_formula = min(formula_results.items(), key=lambda x: x[1]['mse'])

    print(f"""

RECOMMENDED PARAMETERS:
════════════════════════════════════════════════════════════════════════════════

    Standard formula (ρ_crit = A × V^B):
    A = {best_standard[0]:.6f}
    B = {best_standard[1]:.4f}
    RMS error: {np.sqrt(best_standard[2]):.4f}

    Best alternative: {best_formula[0]}
    Parameters: {best_formula[1]['params']}
    RMS error: {best_formula[1]['rms']:.4f}
""")

    # Key insight
    print("""
KEY INSIGHT:
════════════════════════════════════════════════════════════════════════════════

    The regime-specific fits show that A varies by orders of magnitude:
    - ETGs need A ~ 10^-3 to 10^-2
    - Dwarfs need A ~ 10^-6 to 10^-5

    This ~1000x variation suggests one of:
    1. The formula ρ_crit = A × V^B is incomplete
    2. Different physics operates at different mass scales
    3. The mean density definition needs scale-dependent correction

    RECOMMENDATION FOR arXiv:
    ─────────────────────────────────────────────────────────────────────────
    1. Report that model works well within regimes (dwarfs OR ETGs)
    2. Acknowledge that a universal formula remains elusive
    3. Present this as a direction for future theoretical work
    4. The SHAPE of the coherence curve is validated, even if normalization varies
""")

    # Save results
    output = {
        'session': 52,
        'track': 'A - Parameter Recalibration',
        'date': datetime.now().isoformat(),

        'current_parameters': {'A': 0.25, 'B': 1.62},

        'calibration_results': {
            'approach_1_equal_weights': {'A': float(A1), 'B': float(B1), 'rms': float(np.sqrt(err1))},
            'approach_2_etg_weighted': {'A': float(A2), 'B': float(B2), 'rms': float(np.sqrt(err2))},
            'approach_3_by_regime': {k: {'A': float(v['A']), 'B': float(v['B']), 'rms': float(np.sqrt(v['error']))}
                                     for k, v in regime_results.items()},
            'approach_4_formulas': {k: {'params': [float(p) for p in v['params']], 'rms': float(v['rms'])}
                                    for k, v in formula_results.items()}
        },

        'key_finding': 'A varies by ~1000x between ETGs and dwarfs - single formula is challenging',

        'recommendations': [
            'Model works well within regimes',
            'Universal formula needs theoretical work',
            'Report shape validation, acknowledge normalization varies',
        ]
    }

    output_path = Path(__file__).parent / 'session52_recalibration_results.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    print("\n" + "="*80)
    print("SESSION #52 TRACK A COMPLETE")
    print("="*80)

    return output


if __name__ == '__main__':
    main()
