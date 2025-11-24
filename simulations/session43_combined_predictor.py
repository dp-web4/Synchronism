#!/usr/bin/env python3
"""
Session #43 Track C: Combined Virial + γ Predictor

Session #43 Findings:
- Track A: 88% of galaxies at γ bounds [0.10, 2.0] in Session #42
- Track B: Extended bounds [0.001, 20.0] → 66.3% success, still 50% at bounds
- Conclusion: γ optimization is difficult; many galaxies want extreme values

Alternative Strategy:
Instead of optimizing γ per galaxy, test:
1. **Fixed γ**: Use optimal median γ from non-bound galaxies
2. **ρ_crit predicted from v_max**: Use Session #42 virial scaling

This creates a FULLY PREDICTIVE model:
- 0 galaxy-specific parameters (both ρ_crit and γ predicted)
- Compare to Session #42 virial predictor (44.6% with optimized γ)

Models to Test:
A. Baseline tanh: ρ_crit predicted + γ fixed at median
B. Stretched exp: ρ_crit predicted + γ fixed at 0.30
C. Simple exponential: ρ_crit predicted + γ = 0.30 (Session #41 model with prediction)

Author: CBP Autonomous Synchronism Research
Date: 2025-11-24
Session: #43 - Combined Predictor
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime
from scipy import stats

# Import utilities
sys.path.append(str(Path(__file__).parent))
from session38_sparc_refined_coherence import (
    SynchronismPredictor,
    RefinedCoherencePredictor,
    RealSPARCLoader,
    SPARCGalaxy
)


class TanhCoherencePredictor(SynchronismPredictor):
    """Tanh coherence with fixed γ."""

    def __init__(self, rho_crit: float = 1.0, gamma: float = 0.30, beta: float = 0.30):
        super().__init__(gamma=gamma, beta=beta)
        self.rho_crit = rho_crit
        self.tanh_gamma = gamma

    def compute_coherence(self, rho_vis: np.ndarray, rho_0=None) -> np.ndarray:
        C_vis = np.tanh(self.tanh_gamma * np.log(rho_vis / self.rho_crit + 1.0))
        return np.clip(C_vis, 0.0, 0.99999)


def load_virial_params():
    """Load virial scaling parameters from Session #42."""
    results_path = Path(__file__).parent / 'session42_virial_predictor_results.json'

    with open(results_path, 'r') as f:
        data = json.load(f)

    return data['virial_params']


def predict_rho_crit_from_v_max(v_max: float, virial_params: dict) -> float:
    """Predict ρ_crit from v_max using Session #42 virial scaling."""
    A = virial_params['A']
    B = virial_params['B']
    return A * v_max**B


def test_combined_predictor(coherence_type: str = 'tanh',
                             gamma_fixed: float = 1.0,
                             virial_params: dict = None,
                             limit=None):
    """
    Test combined predictor with predicted ρ_crit and fixed γ.

    Args:
        coherence_type: 'tanh' or 'exponential'
        gamma_fixed: Fixed γ value for all galaxies
        virial_params: Virial scaling parameters for ρ_crit prediction
        limit: Limit number of galaxies (for testing)
    """
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies(limit=limit)

    results = []

    for galaxy in galaxies:
        # Compute v_max
        v_max = np.max(galaxy.v_obs)

        # Predict ρ_crit from v_max
        rho_crit_pred = predict_rho_crit_from_v_max(v_max, virial_params)

        # Create predictor with predicted ρ_crit and fixed γ
        if coherence_type == 'tanh':
            predictor = TanhCoherencePredictor(rho_crit=rho_crit_pred, gamma=gamma_fixed)
        else:  # exponential
            predictor = RefinedCoherencePredictor(rho_crit=rho_crit_pred, gamma=gamma_fixed)

        # Fit only α (both ρ_crit and γ are predicted/fixed)
        fit_result = predictor.fit_alpha(galaxy)

        results.append({
            'name': galaxy.name,
            'v_max': float(v_max),
            'rho_crit_predicted': float(rho_crit_pred),
            'gamma_fixed': float(gamma_fixed),
            'alpha': float(fit_result['alpha_best']),
            'chi2_red': float(fit_result['chi2_red'])
        })

    return results


def run_combined_predictor_analysis():
    """Test multiple combined predictor configurations."""

    print("\n" + "="*80)
    print("SESSION #43 TRACK C: COMBINED PREDICTOR ANALYSIS")
    print("="*80)

    # Load virial params
    virial_params = load_virial_params()
    print(f"\nVirial scaling: ρ_crit = {virial_params['A']:.2f} × v_max^{virial_params['B']:.2f}")

    # Test configurations
    configs = [
        ('tanh', 0.30, 'Tanh + γ=0.30'),
        ('tanh', 0.50, 'Tanh + γ=0.50'),
        ('tanh', 1.00, 'Tanh + γ=1.00'),
        ('tanh', 2.00, 'Tanh + γ=2.00'),
        ('tanh', 5.00, 'Tanh + γ=5.00'),
        ('exponential', 0.30, 'Exponential + γ=0.30 (Session #42 virial)'),
    ]

    all_results = {}

    print(f"\nTesting {len(configs)} configurations...")
    print("-"*80)

    for ctype, gamma, label in configs:
        results = test_combined_predictor(
            coherence_type=ctype,
            gamma_fixed=gamma,
            virial_params=virial_params
        )

        chi2_values = np.array([r['chi2_red'] for r in results])
        success_rate = 100.0 * np.sum(chi2_values < 5.0) / len(results)
        median_chi2 = np.median(chi2_values)

        all_results[label] = {
            'coherence_type': ctype,
            'gamma': gamma,
            'success_rate': float(success_rate),
            'median_chi2': float(median_chi2),
            'results': results
        }

        print(f"  {label:40s}: {success_rate:5.1f}% success, median χ² = {median_chi2:.2f}")

    # Summary
    print("\n" + "="*80)
    print("COMBINED PREDICTOR COMPARISON")
    print("="*80)

    print("\n### Success Rates (ZERO galaxy-specific parameters):\n")

    # Sort by success rate
    sorted_configs = sorted(all_results.items(), key=lambda x: x[1]['success_rate'], reverse=True)

    for label, data in sorted_configs:
        print(f"  {label:45s}: {data['success_rate']:5.1f}%  (median χ² = {data['median_chi2']:.2f})")

    # Best configuration
    best_label, best_data = sorted_configs[0]
    print(f"\n### BEST COMBINED PREDICTOR:")
    print(f"  {best_label}")
    print(f"  Success rate: {best_data['success_rate']:.1f}%")
    print(f"  Median χ²: {best_data['median_chi2']:.2f}")
    print(f"  Parameters: ρ_crit = {virial_params['A']:.2f} × v_max^{virial_params['B']:.2f}, γ = {best_data['gamma']}")

    # Compare to baselines
    print(f"\n### Comparison to Baselines:")
    print(f"  Session #42 tanh (optimized γ):     64.6%")
    print(f"  Session #42 virial (optimized γ):   44.6%")
    print(f"  Session #43 combined (fixed γ):     {best_data['success_rate']:.1f}%")

    # Save results
    output_path = Path(__file__).parent / 'session43_combined_predictor_results.json'

    save_data = {
        'session': 43,
        'track': 'C - Combined Predictor',
        'date': datetime.now().isoformat(),
        'virial_params': virial_params,
        'configurations': {
            k: {
                'coherence_type': v['coherence_type'],
                'gamma': v['gamma'],
                'success_rate': v['success_rate'],
                'median_chi2': v['median_chi2'],
                'n_galaxies': len(v['results'])
            }
            for k, v in all_results.items()
        },
        'best_config': {
            'label': best_label,
            'success_rate': best_data['success_rate'],
            'median_chi2': best_data['median_chi2'],
            'gamma': best_data['gamma']
        }
    }

    with open(output_path, 'w') as f:
        json.dump(save_data, f, indent=2)

    print(f"\n\nResults saved to: {output_path}")

    return save_data


if __name__ == '__main__':
    results = run_combined_predictor_analysis()

    print("\n" + "="*80)
    print("SESSION #43 TRACK C COMPLETE")
    print("="*80)
