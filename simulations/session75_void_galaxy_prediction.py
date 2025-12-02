"""
Session #75 Track B: Void Galaxy Dynamics Prediction

From Session #74: Binary pulsars are NOT discriminating tests for Synchronism
because C ~ 1 in all high-density regions.

Discriminating tests need LOW-DENSITY environments where C << 1.
Void galaxies are ideal: galaxies in cosmic voids have lower external density.

Key prediction: Void galaxies should have DIFFERENT Tully-Fisher relation
than cluster galaxies if Synchronism is correct.

Physics:
- C_formation depends on environment during galaxy formation
- Void galaxies form in lower-density environments → lower C_formation
- Lower C → higher G_eff = G/C → different dynamics
"""

import numpy as np
import json
import os

# Physical constants
G = 6.67430e-11         # m³/kg/s²
C_LIGHT = 299792458     # m/s
M_SUN = 1.989e30        # kg
PC = 3.086e16           # parsec in m
KPC = 1000 * PC


class CoherenceModel:
    """
    Coherence model from Sessions #72-74.

    C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))
    """

    def __init__(self, gamma=2.0, rho_crit=1e-22):
        """
        gamma: steepness parameter (derived Session #64)
        rho_crit: critical density kg/m³
        """
        self.gamma = gamma
        self.rho_crit = rho_crit

    def coherence(self, rho):
        """Compute coherence from density."""
        x = np.log(rho / self.rho_crit + 1)
        return np.tanh(self.gamma * x)

    def effective_G(self, rho):
        """Effective gravitational constant G_eff = G/C."""
        C = self.coherence(rho)
        return G / np.maximum(C, 0.01)  # Avoid division by zero


class GalaxyEnvironment:
    """
    Model different galaxy environments.

    Note: These are EXTERNAL environment densities, not internal galaxy density.
    The coherence that matters is the TOTAL density (galaxy + environment)
    at the location of interest.

    For formation coherence, use characteristic density during formation.
    """

    # Typical TOTAL densities at galaxy disk location in kg/m³
    # Galaxy disk ρ ~ 10⁻²¹ to 10⁻²³ kg/m³
    # Environment adds to this

    CLUSTER_CENTER = 1e-21      # Dense cluster - high total density
    CLUSTER_OUTSKIRTS = 5e-22   # Moderate cluster
    FIELD = 1e-22               # Field galaxy
    VOID_CENTER = 1e-23         # Void galaxy - lower formation density

    def __init__(self, environment_type='field'):
        self.type = environment_type
        self.densities = {
            'cluster_center': self.CLUSTER_CENTER,
            'cluster_outskirts': self.CLUSTER_OUTSKIRTS,
            'field': self.FIELD,
            'void': self.VOID_CENTER
        }
        self.rho_env = self.densities.get(environment_type, self.FIELD)

    def formation_coherence(self, model: CoherenceModel):
        """
        Coherence at galaxy formation.

        Galaxies form with coherence determined by environment.
        This "freezes in" at formation and affects subsequent dynamics.
        """
        return model.coherence(self.rho_env)


class TullyFisherRelation:
    """
    Tully-Fisher relation: L ∝ v_max^α or M_baryonic ∝ v_max^β

    Standard: L ∝ v_max^4 (baryonic TF: M_b ∝ v_max^4)

    In Synchronism with G_eff = G/C:
    v² = G_eff M / r = G M / (C r)

    So v² ∝ 1/C for fixed M, r

    If C differs between environments, TF relation shifts!
    """

    def __init__(self, model: CoherenceModel):
        self.model = model

    def standard_TF(self, M_baryonic, r_disk):
        """
        Standard Tully-Fisher: v_max = sqrt(G M / r)
        """
        return np.sqrt(G * M_baryonic / r_disk)

    def synchronism_TF(self, M_baryonic, r_disk, C_formation):
        """
        Synchronism Tully-Fisher: v_max = sqrt(G M / (C r))
        """
        G_eff = G / C_formation
        return np.sqrt(G_eff * M_baryonic / r_disk)

    def predict_velocity_ratio(self, C_void, C_cluster):
        """
        Ratio of void to cluster galaxy velocities at same M, r.

        v_void / v_cluster = sqrt(C_cluster / C_void)

        If C_void < C_cluster, void galaxies have HIGHER v_max!
        """
        return np.sqrt(C_cluster / C_void)


class VoidGalaxyPrediction:
    """
    Quantitative predictions for void vs cluster galaxies.
    """

    def __init__(self):
        self.model = CoherenceModel()
        self.TF = TullyFisherRelation(self.model)

    def compute_coherence_by_environment(self):
        """
        Compute formation coherence for each environment.
        """
        environments = ['cluster_center', 'cluster_outskirts', 'field', 'void']
        results = {}

        for env in environments:
            galaxy_env = GalaxyEnvironment(env)
            C = galaxy_env.formation_coherence(self.model)
            results[env] = {
                'rho_env': float(galaxy_env.rho_env),
                'C_formation': float(C),
                'G_eff_ratio': float(1/C)  # G_eff/G
            }

        return results

    def predict_TF_shift(self):
        """
        Predict shift in Tully-Fisher relation between environments.
        """
        # Coherences
        C_cluster = self.model.coherence(GalaxyEnvironment.CLUSTER_CENTER)
        C_void = self.model.coherence(GalaxyEnvironment.VOID_CENTER)

        # Velocity ratio
        v_ratio = self.TF.predict_velocity_ratio(C_void, C_cluster)

        # In magnitudes (TF is log-log)
        # If v_void = k × v_cluster, then log(v_void) = log(k) + log(v_cluster)
        # TF offset in log-v = log(v_ratio) = 0.5 × log(C_cluster/C_void)

        delta_log_v = np.log10(v_ratio)

        # Standard TF: M ∝ v^4, so log(M) = 4 log(v) + const
        # Shift in log(M) at fixed v = -4 × delta_log_v

        delta_log_M = -4 * delta_log_v

        return {
            'C_cluster': float(C_cluster),
            'C_void': float(C_void),
            'v_ratio_void_to_cluster': float(v_ratio),
            'delta_log_v': float(delta_log_v),
            'delta_log_M_at_fixed_v': float(delta_log_M),
            'interpretation': f'''
Void galaxies have C_formation = {C_void:.4f} vs cluster {C_cluster:.4f}.

At fixed baryonic mass:
- Void galaxy v_max is {(v_ratio-1)*100:.1f}% HIGHER than cluster galaxy

At fixed v_max:
- Void galaxy appears {10**delta_log_M:.2f}x more massive in baryonic TF

This is a TESTABLE PREDICTION:
Void galaxies should lie ABOVE the standard baryonic Tully-Fisher relation
(higher v at fixed M, or appearing more massive at fixed v).
            '''
        }

    def rotation_curve_prediction(self):
        """
        Predict rotation curve differences.
        """
        # Test galaxy: M = 10^10 M_sun, r_disk = 3 kpc
        M = 1e10 * M_SUN
        r_disk = 3 * KPC

        # Radii from 1 to 30 kpc
        r = np.linspace(1, 30, 100) * KPC

        # Standard rotation curve (Keplerian at large r)
        v_std = np.sqrt(G * M / r)

        # Synchronism: G_eff = G/C(r)
        # Need density profile to compute C(r)
        # Use exponential disk: ρ(r) = ρ_0 exp(-r/r_disk)
        rho_0_cluster = 1e-22  # cluster galaxy
        rho_0_void = 1e-24     # void galaxy

        rho_cluster = rho_0_cluster * np.exp(-r / r_disk)
        rho_void = rho_0_void * np.exp(-r / r_disk)

        C_cluster = self.model.coherence(rho_cluster)
        C_void = self.model.coherence(rho_void)

        v_cluster = np.sqrt(G * M / (C_cluster * r))
        v_void = np.sqrt(G * M / (C_void * r))

        return {
            'r_kpc': (r / KPC).tolist(),
            'v_standard_km_s': (v_std / 1000).tolist(),
            'v_cluster_km_s': (v_cluster / 1000).tolist(),
            'v_void_km_s': (v_void / 1000).tolist(),
            'interpretation': '''
Rotation curve predictions:

Standard: v(r) = sqrt(GM/r) - falls as r^(-1/2)

Synchronism cluster galaxy:
- Inner regions: C ~ 1, v ~ v_standard
- Outer regions: C < 1, v > v_standard (flat curve)

Synchronism void galaxy:
- C is LOWER throughout (formed in low-density)
- v is HIGHER throughout
- Even FLATTER rotation curve

Observable difference:
- Void galaxies: more extended flat regions
- Higher v_max at fixed luminosity
- Different dark matter halo inference
            '''
        }


def compute_predictions():
    """
    Compute all void galaxy predictions.
    """
    results = {
        'session': 75,
        'track': 'B',
        'title': 'Void Galaxy Dynamics Prediction'
    }

    predictor = VoidGalaxyPrediction()

    # Environment coherences
    env_results = predictor.compute_coherence_by_environment()
    results['environment_coherence'] = env_results

    # TF shift prediction
    tf_results = predictor.predict_TF_shift()
    results['tully_fisher_shift'] = tf_results

    # Rotation curve prediction
    rc_results = predictor.rotation_curve_prediction()
    results['rotation_curves'] = rc_results

    # Summary
    results['summary'] = '''
VOID GALAXY PREDICTIONS - DISCRIMINATING TEST FOR SYNCHRONISM

PREDICTION 1: Tully-Fisher Offset
- Void galaxies have {:.1f}% higher v_max at fixed baryonic mass
- Or appear {:.1f}x more massive at fixed v_max
- Testable with SDSS void galaxy sample + HI surveys

PREDICTION 2: Rotation Curve Shape
- Void galaxies have more extended flat regions
- Higher asymptotic velocity
- Different inferred dark matter halo

PREDICTION 3: Dark Matter Halo Inference
- Standard ΛCDM: halo mass independent of environment
- Synchronism: apparent halo is coherence effect
- Void galaxies should have "larger" apparent halos

FALSIFICATION CRITERIA:
- If void and cluster galaxies show IDENTICAL TF relation: Synchronism falsified
- If void galaxies have LOWER v_max: Synchronism falsified
- Prediction requires >{:.0f}% velocity difference to be detectable

CURRENT OBSERVATIONAL STATUS:
- Some studies suggest void galaxies may have unusual properties
- Need systematic comparison with explicit coherence model
- SDSS + ALFALFA data could test this
    '''.format(
        (tf_results['v_ratio_void_to_cluster'] - 1) * 100,
        10**tf_results['delta_log_M_at_fixed_v'],
        (tf_results['v_ratio_void_to_cluster'] - 1) * 100
    )

    return results


if __name__ == '__main__':
    print("="*60)
    print("Session #75 Track B: Void Galaxy Predictions")
    print("="*60)

    results = compute_predictions()

    print("\n--- Environment Coherences ---")
    for env, data in results['environment_coherence'].items():
        print(f"{env:20s}: C = {data['C_formation']:.4f}, G_eff/G = {data['G_eff_ratio']:.2f}")

    print("\n--- Tully-Fisher Shift ---")
    tf = results['tully_fisher_shift']
    print(f"C_cluster = {tf['C_cluster']:.4f}")
    print(f"C_void = {tf['C_void']:.4f}")
    print(f"v_ratio (void/cluster) = {tf['v_ratio_void_to_cluster']:.4f}")
    print(f"At fixed M: void v_max is {(tf['v_ratio_void_to_cluster']-1)*100:.1f}% higher")

    print("\n--- Summary ---")
    print(results['summary'][:500] + "...")

    # Save
    output_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(output_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    output_file = os.path.join(results_dir, 'session75_void_galaxy_prediction.json')

    # Clean rotation curve arrays for JSON
    rc = results['rotation_curves']
    rc['r_kpc'] = [float(x) for x in rc['r_kpc'][:10]] + ['...']  # Truncate for readability

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults saved to: {output_file}")
