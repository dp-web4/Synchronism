#!/usr/bin/env python3
"""
Session #42: Alternative Coherence Functions - Addressing 30.3% Bound-Hitting

Session #41 Results:
- Success rate: 56.0% (χ² < 5)
- 30.3% still at upper bound (ρ_crit = 100000)
- Discovered ρ_crit ∝ v_max^1.74 (virial scaling)

Session #41 Recommendations:
- Test alternative coherence functions
- Multi-parameter optimization
- Leverage v_max correlation for predictive modeling

This session explores 4 alternative coherence function forms to improve
performance and reduce bound-hitting.

Alternative Forms Tested:
1. **Current (baseline)**: C = 1 - exp(-(ρ/ρ_c)^γ)  [γ = 0.30 fixed]
2. **Stretched exponential**: C = 1 - exp(-(ρ/ρ_c)^γ)  [γ FREE parameter]
3. **Power-law cutoff**: C = (ρ/ρ_c)^γ / [1 + (ρ/ρ_c)^γ]
4. **Tanh (smooth saturation)**: C = tanh(γ × log(ρ/ρ_c + 1))
5. **Double exponential**: C = 1 - exp(-ρ/ρ_c1) × exp(-(ρ/ρ_c2)^γ)

Author: CBP Autonomous Synchronism Research
Date: 2025-11-23
Session: #42 - Alternative Coherence Functions
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy.optimize import minimize, minimize_scalar
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple

# Import SPARC validation utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    SynchronismPredictor,
    RealSPARCLoader,
    SPARCGalaxy
)


@dataclass
class CoherenceFunction:
    """Define a coherence function with its properties."""
    name: str
    formula: str
    compute: Callable  # Function signature: (rho_vis, params) -> C_vis
    param_names: List[str]
    param_bounds: List[Tuple[float, float]]
    description: str


class AlternativeCoherencePredictor(SynchronismPredictor):
    """
    Predictor with pluggable coherence functions.

    Allows testing different functional forms beyond exponential.
    """

    def __init__(self, coherence_func: CoherenceFunction,
                 params: Dict[str, float],
                 beta: float = 0.30):
        """
        Initialize with custom coherence function.

        Args:
            coherence_func: CoherenceFunction instance
            params: Dictionary of parameter values
            beta: Dark matter scaling exponent (fixed at 0.30)
        """
        super().__init__(gamma=0.30, beta=beta)  # gamma not used for alternatives
        self.coherence_func = coherence_func
        self.params = params

    def compute_coherence(self, rho_vis: np.ndarray, rho_0=None) -> np.ndarray:
        """Compute coherence using alternative function."""
        C_vis = self.coherence_func.compute(rho_vis, self.params)

        # Ensure 0 ≤ C < 1
        C_vis = np.clip(C_vis, 0.0, 0.99999)

        return C_vis


# ============================================================================
# Alternative Coherence Functions
# ============================================================================

def baseline_exponential(rho_vis: np.ndarray, params: Dict) -> np.ndarray:
    """
    BASELINE: Current Session #38-41 formula

    C = 1 - exp(-(ρ/ρ_crit)^γ)

    Properties:
    - Smooth transition from 0 to 1
    - Never fully saturates
    - γ = 0.30 fixed (empirical from Session #34)
    """
    rho_c = params['rho_crit']
    gamma = params.get('gamma', 0.30)  # Default to 0.30 if not provided

    return 1.0 - np.exp(-(rho_vis / rho_c) ** gamma)


def stretched_exponential(rho_vis: np.ndarray, params: Dict) -> np.ndarray:
    """
    VARIANT 1: Stretched exponential with FREE γ

    C = 1 - exp(-(ρ/ρ_crit)^γ)

    Difference from baseline: γ is optimized, not fixed at 0.30

    Test: Does optimal γ vary by galaxy? Or is 0.30 universal?
    """
    rho_c = params['rho_crit']
    gamma = params['gamma']  # FREE parameter

    return 1.0 - np.exp(-(rho_vis / rho_c) ** gamma)


def power_law_cutoff(rho_vis: np.ndarray, params: Dict) -> np.ndarray:
    """
    VARIANT 2: Power-law with cutoff (similar to Fermi-Dirac)

    C = (ρ/ρ_crit)^γ / [1 + (ρ/ρ_crit)^γ]

    Properties:
    - C ~ (ρ/ρ_c)^γ for ρ << ρ_c (power-law rise)
    - C → 1 - (ρ_c/ρ)^γ for ρ >> ρ_c (smooth saturation)
    - γ controls transition sharpness

    Physical motivation: Occupation probability (quantum analogy)
    """
    rho_c = params['rho_crit']
    gamma = params['gamma']

    ratio = (rho_vis / rho_c) ** gamma
    return ratio / (1.0 + ratio)


def tanh_log(rho_vis: np.ndarray, params: Dict) -> np.ndarray:
    """
    VARIANT 3: Tanh of logarithm (smooth symmetric transition)

    C = tanh(γ × log(ρ/ρ_crit + 1))

    Properties:
    - Symmetric S-curve in log-space
    - C → 0 as ρ → 0
    - C → 1 as ρ → ∞ (but slowly)
    - γ controls transition width

    Physical motivation: Logistic growth / phase transition
    """
    rho_c = params['rho_crit']
    gamma = params['gamma']

    # Add 1 to avoid log(0)
    return np.tanh(gamma * np.log(rho_vis / rho_c + 1.0))


def double_exponential(rho_vis: np.ndarray, params: Dict) -> np.ndarray:
    """
    VARIANT 4: Double exponential (two-scale model)

    C = 1 - exp(-ρ/ρ_c1) × exp(-(ρ/ρ_c2)^γ)

    Properties:
    - Two characteristic densities: ρ_c1 (linear), ρ_c2 (nonlinear)
    - Low ρ: C ~ ρ/ρ_c1 (linear rise)
    - High ρ: C ~ 1 - exp(-(ρ/ρ_c2)^γ) (standard form)

    Physical motivation: Two-stage coherence (local vs global)
    """
    rho_c1 = params['rho_crit']
    rho_c2 = params.get('rho_crit2', rho_c1)  # Second scale
    gamma = params.get('gamma', 0.30)

    term1 = np.exp(-rho_vis / rho_c1)
    term2 = np.exp(-(rho_vis / rho_c2) ** gamma)

    return 1.0 - term1 * term2


# ============================================================================
# Define Coherence Function Library
# ============================================================================

COHERENCE_FUNCTIONS = {
    'baseline': CoherenceFunction(
        name='Baseline Exponential',
        formula='C = 1 - exp(-(ρ/ρ_c)^0.30)',
        compute=baseline_exponential,
        param_names=['rho_crit'],
        param_bounds=[(0.01, 100000)],
        description='Current Session #38-41 model (γ = 0.30 fixed)'
    ),

    'stretched': CoherenceFunction(
        name='Stretched Exponential',
        formula='C = 1 - exp(-(ρ/ρ_c)^γ)',
        compute=stretched_exponential,
        param_names=['rho_crit', 'gamma'],
        param_bounds=[(0.01, 100000), (0.10, 1.0)],
        description='FREE γ exponent (test universality of 0.30)'
    ),

    'powerlaw': CoherenceFunction(
        name='Power-Law Cutoff',
        formula='C = (ρ/ρ_c)^γ / [1 + (ρ/ρ_c)^γ]',
        compute=power_law_cutoff,
        param_names=['rho_crit', 'gamma'],
        param_bounds=[(0.01, 100000), (0.10, 1.0)],
        description='Fermi-Dirac-like occupation probability'
    ),

    'tanh': CoherenceFunction(
        name='Tanh Log',
        formula='C = tanh(γ × log(ρ/ρ_c + 1))',
        compute=tanh_log,
        param_names=['rho_crit', 'gamma'],
        param_bounds=[(0.01, 100000), (0.10, 2.0)],
        description='Smooth symmetric S-curve transition'
    ),
}


# ============================================================================
# Optimization Functions
# ============================================================================

def fit_galaxy_alternative(galaxy: SPARCGalaxy,
                          coherence_func: CoherenceFunction,
                          method='continuous') -> Dict:
    """
    Fit galaxy with alternative coherence function.

    Args:
        galaxy: SPARC galaxy data
        coherence_func: CoherenceFunction to test
        method: 'continuous' (default) or 'grid'

    Returns:
        Dictionary with fit results and parameters
    """

    def objective(params_array):
        """Objective function for optimization."""
        # Convert array to parameter dictionary
        params_dict = {
            name: val
            for name, val in zip(coherence_func.param_names, params_array)
        }

        # Create predictor with these parameters
        predictor = AlternativeCoherencePredictor(
            coherence_func=coherence_func,
            params=params_dict,
            beta=0.30
        )

        # Fit alpha and compute χ²
        result = predictor.fit_alpha(galaxy)
        return result['chi2_red']

    # Initial guess (center of log-space)
    x0 = []
    for bounds in coherence_func.param_bounds:
        if bounds[0] > 0:
            # Log-space midpoint
            x0.append(np.sqrt(bounds[0] * bounds[1]))
        else:
            # Linear midpoint
            x0.append((bounds[0] + bounds[1]) / 2)

    # Optimize
    if len(coherence_func.param_names) == 1:
        # 1D optimization (faster)
        result = minimize_scalar(
            lambda x: objective([x]),
            bounds=coherence_func.param_bounds[0],
            method='bounded',
            options={'xatol': 1e-3}
        )
        params_opt = [result.x]
        chi2_opt = result.fun
        nfev = result.nfev
        success = result.success
    else:
        # Multi-dimensional optimization
        result = minimize(
            objective,
            x0=x0,
            bounds=coherence_func.param_bounds,
            method='L-BFGS-B',
            options={'ftol': 1e-6, 'maxiter': 100}
        )
        params_opt = result.x
        chi2_opt = result.fun
        nfev = result.nfev
        success = result.success

    # Get full fit at optimal parameters
    params_dict = {
        name: val
        for name, val in zip(coherence_func.param_names, params_opt)
    }

    predictor_opt = AlternativeCoherencePredictor(
        coherence_func=coherence_func,
        params=params_dict,
        beta=0.30
    )

    fit_result = predictor_opt.fit_alpha(galaxy)

    # Check bounds
    at_bounds = {}
    for i, (name, val) in enumerate(zip(coherence_func.param_names, params_opt)):
        bounds = coherence_func.param_bounds[i]
        at_bounds[f'{name}_at_lower'] = abs(val - bounds[0]) < 0.01 * bounds[0]
        at_bounds[f'{name}_at_upper'] = abs(val - bounds[1]) < 0.01 * bounds[1]

    # Compile results (convert numpy types to native Python for JSON)
    return {
        'name': galaxy.name,
        'coherence_function': coherence_func.name,
        'params': {k: float(v) for k, v in params_dict.items()},
        'chi2_red': float(chi2_opt),
        'alpha': float(fit_result['alpha_best']),
        'nfev': int(nfev),
        'success': bool(success),
        **{k: bool(v) for k, v in at_bounds.items()}
    }


def compare_all_functions(galaxy: SPARCGalaxy) -> Dict:
    """Compare all coherence functions on a single galaxy."""
    results = {}

    for func_key, func in COHERENCE_FUNCTIONS.items():
        results[func_key] = fit_galaxy_alternative(galaxy, func)

    return results


def run_full_comparison_sparc(limit=None, output_file='session42_alternative_coherence_results.json'):
    """
    Test all alternative coherence functions on SPARC dataset.

    Compares:
    - Baseline (γ = 0.30 fixed)
    - Stretched (γ free)
    - Power-law cutoff
    - Tanh log

    For each galaxy, finds best function by AIC/BIC.
    """
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=limit)

    print("\n" + "="*80)
    print("SESSION #42: ALTERNATIVE COHERENCE FUNCTIONS")
    print("="*80)
    print(f"\nTesting {len(COHERENCE_FUNCTIONS)} coherence function forms")
    print(f"on {len(galaxies)} SPARC galaxies\n")

    for key, func in COHERENCE_FUNCTIONS.items():
        print(f"  {key:10s}: {func.formula}")

    print("\n" + "-"*80)

    all_results = []

    for i, galaxy in enumerate(galaxies, 1):
        if i % 25 == 0:
            print(f"Processed {i}/{len(galaxies)} galaxies...")

        galaxy_results = compare_all_functions(galaxy)
        all_results.append({
            'galaxy': galaxy.name,
            'results': galaxy_results
        })

    # Analyze results
    print("\n" + "="*80)
    print("SESSION #42 ALTERNATIVE COHERENCE ANALYSIS")
    print("="*80)

    # Success rates by function
    print("\n### Success Rates (χ² < 5):\n")

    for func_key in COHERENCE_FUNCTIONS.keys():
        successes = sum(
            1 for r in all_results
            if r['results'][func_key]['chi2_red'] < 5.0
        )
        success_rate = 100.0 * successes / len(all_results)
        print(f"  {COHERENCE_FUNCTIONS[func_key].name:25s}: {success_rate:5.1f}% ({successes}/{len(all_results)})")

    # Median χ² by function
    print("\n### Median χ²:\n")

    for func_key in COHERENCE_FUNCTIONS.keys():
        chi2_values = [r['results'][func_key]['chi2_red'] for r in all_results]
        median_chi2 = np.median(chi2_values)
        print(f"  {COHERENCE_FUNCTIONS[func_key].name:25s}: {median_chi2:6.2f}")

    # Bound-hitting analysis (for multi-parameter functions)
    print("\n### Parameter Bound-Hitting:\n")

    for func_key in COHERENCE_FUNCTIONS.keys():
        func = COHERENCE_FUNCTIONS[func_key]

        for param_name in func.param_names:
            at_upper = sum(
                1 for r in all_results
                if r['results'][func_key].get(f'{param_name}_at_upper', False)
            )
            at_lower = sum(
                1 for r in all_results
                if r['results'][func_key].get(f'{param_name}_at_lower', False)
            )

            if at_upper > 0 or at_lower > 0:
                print(f"  {func.name:25s} ({param_name:10s}): {at_upper:3d} at upper, {at_lower:3d} at lower")

    # Best function per galaxy
    print("\n### Best Function Distribution:\n")

    best_counts = {k: 0 for k in COHERENCE_FUNCTIONS.keys()}

    for result in all_results:
        best_func = min(
            COHERENCE_FUNCTIONS.keys(),
            key=lambda k: result['results'][k]['chi2_red']
        )
        best_counts[best_func] += 1

    for func_key, count in best_counts.items():
        pct = 100.0 * count / len(all_results)
        print(f"  {COHERENCE_FUNCTIONS[func_key].name:25s}: {count:3d} galaxies ({pct:5.1f}%)")

    # Save results
    output_path = Path(__file__).parent / output_file

    save_data = {
        'session': 42,
        'date': datetime.now().isoformat(),
        'n_galaxies': len(all_results),
        'functions_tested': {
            k: {
                'name': v.name,
                'formula': v.formula,
                'n_params': len(v.param_names)
            }
            for k, v in COHERENCE_FUNCTIONS.items()
        },
        'results': all_results
    }

    with open(output_path, 'w') as f:
        json.dump(save_data, f, indent=2)

    print(f"\nResults saved to: {output_path}")
    print("\n" + "="*80)
    print("SESSION #42 COMPLETE")
    print("="*80)

    return all_results


if __name__ == '__main__':
    # Run full SPARC comparison
    results = run_full_comparison_sparc(limit=None)  # All 175 galaxies
