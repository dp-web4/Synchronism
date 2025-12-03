#!/usr/bin/env python3
"""
Compare Empirical vs First-Principles Derived Models on SPARC

The manuscript versions use different parameter sets:
- arxiv-v1.tex (Nov 25): A = 0.25, B = 1.62 (empirically fitted)
- draft_v1.md (Dec 1): A = 0.028, B = 0.5 (first-principles derived)

This script tests BOTH on SPARC rotation curves to answer:
"How would the derived model perform vs the empirical one?"

Author: CBP Autonomous Synchronism Research
Date: 2025-12-02
"""

import numpy as np
from pathlib import Path
import sys
import json
from datetime import datetime

# Import utilities
sys.path.append(str(Path(__file__).parent))
from synchronism_real_sparc_validation import RealSPARCLoader
from synchronism_sparc_validation import SPARCGalaxy, SynchronismPredictor


class TanhCoherenceWithVirial(SynchronismPredictor):
    """
    Tanh coherence with virial-based rho_crit prediction.

    C = tanh(gamma * log(rho_vis / rho_crit + 1))
    rho_crit = A * v_max^B
    """

    def __init__(self, A: float, B: float, gamma: float = 2.0, beta: float = 0.30):
        """
        Initialize predictor.

        Args:
            A: Normalization for rho_crit = A * v_max^B
            B: Exponent for velocity scaling
            gamma: Decoherence exponent (fixed at 2.0 in both models)
            beta: Dark matter scaling exponent
        """
        super().__init__(gamma=gamma, beta=beta)
        self.A = A
        self.B = B
        self.tanh_gamma = gamma

    def predict_rho_crit(self, v_max: float) -> float:
        """Predict critical density from v_max."""
        return self.A * (v_max ** self.B)

    def compute_coherence(self, rho_vis: np.ndarray, rho_0=None) -> np.ndarray:
        """Compute coherence using tanh formula with predicted rho_crit."""
        # C = tanh(gamma * log(rho/rho_crit + 1))
        arg = self.tanh_gamma * np.log(rho_vis / self.rho_crit + 1.0)
        C = np.tanh(arg)
        return np.clip(C, 0.0, 0.99999)

    def fit_alpha_with_virial(self, galaxy: SPARCGalaxy,
                               alpha_range=np.logspace(-2, 2, 100)) -> dict:
        """
        Fit alpha parameter for galaxy using virial-predicted rho_crit.

        Returns dict with best alpha and chi2.
        """
        # Compute v_max and predict rho_crit
        v_max = np.max(galaxy.v_obs)
        self.rho_crit = self.predict_rho_crit(v_max)

        best_chi2 = np.inf
        best_alpha = 1.0

        for alpha in alpha_range:
            # Get visible matter density
            rho_vis = galaxy.total_baryonic_density()

            # Compute coherence with virial-predicted rho_crit
            C_vis = self.compute_coherence(rho_vis)

            # Compute dark matter density
            rho_DM = self.compute_dark_matter_density(rho_vis, C_vis, alpha)

            # Total density and enclosed mass
            rho_total = rho_vis + rho_DM
            M_enc = self.compute_enclosed_mass(galaxy.radius, rho_total)

            # Rotation velocity
            v_pred = np.sqrt(self.G * M_enc / np.maximum(galaxy.radius, 1e-6))

            # Chi-squared
            residuals = (galaxy.v_obs - v_pred) ** 2
            weights = 1.0 / (galaxy.v_err ** 2 + 1e-6)
            chi2 = np.sum(residuals * weights)
            dof = len(galaxy.v_obs) - 1  # 1 free parameter (alpha)
            chi2_red = chi2 / max(dof, 1)

            if chi2_red < best_chi2:
                best_chi2 = chi2_red
                best_alpha = alpha

        return {
            'alpha_best': float(best_alpha),
            'chi2_red': float(best_chi2),
            'n_points': len(galaxy.v_obs),
            'v_max': float(v_max),
            'rho_crit': float(self.rho_crit)
        }


def run_comparison():
    """Compare empirical vs derived models on SPARC."""

    print("\n" + "="*80)
    print("COMPARING EMPIRICAL VS DERIVED MODELS ON SPARC")
    print("="*80)

    # Model configurations
    models = {
        'Empirical (arxiv-v1)': {
            'A': 0.25,
            'B': 1.62,
            'gamma': 2.0,
            'source': 'Session #42-43 fitted parameters'
        },
        'Derived (draft_v1)': {
            'A': 0.028,
            'B': 0.5,
            'gamma': 2.0,
            'source': 'First-principles: A from Jeans, B from virial+TF'
        }
    }

    # Load SPARC galaxies
    loader = RealSPARCLoader()
    galaxies = loader.load_all_galaxies()

    if not galaxies:
        print("ERROR: No SPARC galaxies loaded!")
        return None

    print(f"\nLoaded {len(galaxies)} SPARC galaxies")
    print("-" * 80)

    results = {model_name: [] for model_name in models}

    # Test each model
    for model_name, params in models.items():
        print(f"\nTesting: {model_name}")
        print(f"  A = {params['A']}, B = {params['B']}, γ = {params['gamma']}")
        print(f"  Source: {params['source']}")

        predictor = TanhCoherenceWithVirial(
            A=params['A'],
            B=params['B'],
            gamma=params['gamma']
        )

        for i, galaxy in enumerate(galaxies):
            try:
                fit = predictor.fit_alpha_with_virial(galaxy)

                results[model_name].append({
                    'name': galaxy.name,
                    'v_max': fit['v_max'],
                    'alpha': fit['alpha_best'],
                    'chi2_red': fit['chi2_red'],
                    'n_points': fit['n_points'],
                    'rho_crit': fit['rho_crit']
                })

                if (i + 1) % 25 == 0:
                    print(f"    Processed {i+1}/{len(galaxies)}...")

            except Exception as e:
                print(f"    Error on {galaxy.name}: {e}")

    # Compute summary statistics
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)

    summary = {}

    for model_name in models:
        model_results = results[model_name]
        chi2_values = np.array([r['chi2_red'] for r in model_results])
        v_max_values = np.array([r['v_max'] for r in model_results])

        # Overall success rate (chi2 < 5)
        success = np.sum(chi2_values < 5.0)
        success_rate = 100.0 * success / len(model_results)

        # By velocity class
        dwarf_mask = v_max_values < 50
        intermediate_mask = (v_max_values >= 50) & (v_max_values < 100)
        massive_mask = v_max_values >= 100

        dwarf_success = 100.0 * np.sum(chi2_values[dwarf_mask] < 5.0) / max(np.sum(dwarf_mask), 1)
        intermediate_success = 100.0 * np.sum(chi2_values[intermediate_mask] < 5.0) / max(np.sum(intermediate_mask), 1)
        massive_success = 100.0 * np.sum(chi2_values[massive_mask] < 5.0) / max(np.sum(massive_mask), 1)

        summary[model_name] = {
            'n_galaxies': len(model_results),
            'success_rate': float(success_rate),
            'median_chi2': float(np.median(chi2_values)),
            'mean_chi2': float(np.mean(chi2_values)),
            'dwarf_success': float(dwarf_success),
            'dwarf_n': int(np.sum(dwarf_mask)),
            'intermediate_success': float(intermediate_success),
            'intermediate_n': int(np.sum(intermediate_mask)),
            'massive_success': float(massive_success),
            'massive_n': int(np.sum(massive_mask))
        }

        print(f"\n{model_name}:")
        print(f"  Overall: {success_rate:.1f}% success ({success}/{len(model_results)})")
        print(f"  Median χ²: {np.median(chi2_values):.2f}")
        print(f"  ")
        print(f"  By velocity class:")
        print(f"    Dwarfs (v < 50):     {dwarf_success:.1f}% ({np.sum(dwarf_mask)} galaxies)")
        print(f"    Intermediate (50-100): {intermediate_success:.1f}% ({np.sum(intermediate_mask)} galaxies)")
        print(f"    Massive (v > 100):   {massive_success:.1f}% ({np.sum(massive_mask)} galaxies)")

    # Direct comparison
    print("\n" + "="*80)
    print("HEAD-TO-HEAD COMPARISON")
    print("="*80)

    empirical_chi2 = np.array([r['chi2_red'] for r in results['Empirical (arxiv-v1)']])
    derived_chi2 = np.array([r['chi2_red'] for r in results['Derived (draft_v1)']])

    empirical_wins = np.sum(empirical_chi2 < derived_chi2)
    derived_wins = np.sum(derived_chi2 < empirical_chi2)
    ties = np.sum(empirical_chi2 == derived_chi2)

    print(f"\nPer-galaxy comparison:")
    print(f"  Empirical better: {empirical_wins} galaxies")
    print(f"  Derived better:   {derived_wins} galaxies")
    print(f"  Ties:             {ties} galaxies")

    chi2_ratio = derived_chi2 / (empirical_chi2 + 1e-10)
    print(f"\nMedian χ² ratio (derived/empirical): {np.median(chi2_ratio):.2f}")
    print(f"  <1 means derived is better, >1 means empirical is better")

    # Save results
    output = {
        'date': datetime.now().isoformat(),
        'n_galaxies': len(galaxies),
        'models': {
            model_name: {
                'params': models[model_name],
                'summary': summary[model_name]
            }
            for model_name in models
        },
        'comparison': {
            'empirical_wins': int(empirical_wins),
            'derived_wins': int(derived_wins),
            'ties': int(ties),
            'median_chi2_ratio': float(np.median(chi2_ratio))
        },
        'detailed_results': results
    }

    output_path = Path(__file__).parent / 'empirical_vs_derived_results.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to: {output_path}")

    return output


if __name__ == '__main__':
    results = run_comparison()

    print("\n" + "="*80)
    print("COMPARISON COMPLETE")
    print("="*80)
    print("""
CONCLUSION:
═══════════════════════════════════════════════════════════════════════════════

This comparison answers the question:
"Why wasn't the derived model (A=0.028, B=0.5) tested on SPARC rotation curves?"

The results above show the performance difference between:
- Empirical parameters (fitted to maximize SPARC success)
- Derived parameters (from first principles)

If derived model performs significantly worse, it suggests the theoretical
derivation has gaps that the empirical fitting compensates for.

If derived model performs comparably, it validates the theoretical approach.

═══════════════════════════════════════════════════════════════════════════════
""")
