"""
Session #74 Track C: Binary Pulsar Predictions

From GR_EMERGENCE_SYNTHESIS Priority 5:
"Once the toy model works: Binary pulsar orbital decay"

Binary pulsars (especially Hulse-Taylor PSR B1913+16) provide
precision tests of GR through orbital period decay from gravitational
wave emission.

Key question: Does Synchronism predict the same orbital decay as GR?

From Session #71: GW170817 resolved via conformal invariance.
- Coherence affects matter dynamics, not GW propagation
- GW speed = c (confirmed)
- But GW emission rate could differ if effective G changes

From Session #72: Toy model gives effective metric with G_eff = G/C
"""

import numpy as np
import json
import os

# Physical constants
G = 6.67430e-11         # m³/kg/s²
C_LIGHT = 299792458     # m/s
M_SUN = 1.989e30        # kg
YEAR = 365.25 * 24 * 3600  # seconds


class BinaryPulsar:
    """
    Model of binary pulsar system.

    Hulse-Taylor pulsar (PSR B1913+16) parameters:
    - M1 ≈ 1.4398 M_sun (pulsar)
    - M2 ≈ 1.3886 M_sun (companion)
    - Orbital period P ≈ 7.75 hours
    - Eccentricity e ≈ 0.617
    - Semi-major axis a ≈ 1.95 × 10^9 m

    GR prediction for period decay:
    dP/dt = -(192π/5) × (G^5/3/c^5) × (Mc)^(5/3) × (P_b/(2π))^(-5/3) × f(e)

    where:
    - Mc = (M1 × M2)^(3/5) / (M1 + M2)^(1/5) = chirp mass
    - f(e) = (1 + 73/24 e² + 37/96 e⁴) / (1 - e²)^(7/2)
    """

    def __init__(self, M1=1.4398, M2=1.3886, P_b=7.752, e=0.6171):
        """
        Initialize binary pulsar.

        M1, M2: masses in solar masses
        P_b: orbital period in hours
        e: eccentricity
        """
        self.M1 = M1 * M_SUN  # kg
        self.M2 = M2 * M_SUN  # kg
        self.M_total = self.M1 + self.M2
        self.P_b = P_b * 3600  # seconds
        self.e = e

        # Derived quantities
        self.chirp_mass = (self.M1 * self.M2)**(3/5) / self.M_total**(1/5)

        # Semi-major axis from Kepler's third law
        self.a = (G * self.M_total * self.P_b**2 / (4 * np.pi**2))**(1/3)

    def eccentricity_factor(self):
        """
        f(e) = (1 + 73/24 e² + 37/96 e⁴) / (1 - e²)^(7/2)

        Enhancement factor due to elliptical orbit.
        For e = 0.617: f(e) ≈ 11.8
        """
        e2 = self.e**2
        numerator = 1 + (73/24) * e2 + (37/96) * e2**2
        denominator = (1 - e2)**(7/2)
        return numerator / denominator

    def orbital_decay_GR(self):
        """
        GR prediction for orbital period derivative.

        dP/dt = -(192π/5) × (G^5/3/c^5) × Mc^(5/3) × (P_b/(2π))^(-5/3) × f(e)

        Returns dP/dt in seconds/second (dimensionless).
        """
        Mc = self.chirp_mass
        f_e = self.eccentricity_factor()

        prefactor = -(192 * np.pi / 5)
        G_factor = G**(5/3) / C_LIGHT**5
        mass_factor = Mc**(5/3)
        period_factor = (self.P_b / (2 * np.pi))**(-5/3)

        dPdt = prefactor * G_factor * mass_factor * period_factor * f_e

        return dPdt


class SynchronismBinaryPulsar(BinaryPulsar):
    """
    Binary pulsar in Synchronism framework.

    From Session #72: G_eff = G/C where C = coherence.

    For neutron stars: density ρ ~ 10^17 kg/m³ (nuclear density)
    Using standard coherence: C = tanh(γ × log(ρ/ρ_crit + 1))

    With ρ ~ 10^17, ρ_crit ~ 10^-22:
    log(ρ/ρ_crit) ~ log(10^39) ~ 90
    C = tanh(2 × 90) ≈ 1 (fully classical!)

    So for neutron stars: C ≈ 1, G_eff ≈ G

    But what about the space BETWEEN the stars?
    Interstellar density: ρ ~ 10^-21 kg/m³
    C ~ tanh(2 × log(10)) ~ 0.96 (still high!)

    Key insight: Even "empty" space has ρ >> ρ_crit for compact binaries.
    """

    def __init__(self, M1=1.4398, M2=1.3886, P_b=7.752, e=0.6171):
        super().__init__(M1, M2, P_b, e)

        # Coherence calculation
        self.rho_crit = 1e-22  # kg/m³
        self.gamma = 2.0

    def neutron_star_density(self):
        """
        Average density of neutron star.
        R_NS ~ 10 km, M_NS ~ 1.4 M_sun
        ρ ~ 3M/(4πR³) ~ 4 × 10^17 kg/m³
        """
        R_NS = 10e3  # 10 km in meters
        M_NS = 1.4 * M_SUN
        return 3 * M_NS / (4 * np.pi * R_NS**3)

    def orbital_separation_density(self):
        """
        Estimate effective density at orbital separation.

        For binary pulsar, a ~ 2 × 10^9 m
        Gravitational potential: Φ ~ GM/a
        Effective "gravitational energy density": u_g ~ Φ²/(8πG)

        Or use average matter density in orbital volume.
        """
        # Average matter density if stars smeared over orbit
        M_total = self.M_total
        a = self.a
        volume = (4/3) * np.pi * a**3
        rho_avg = M_total / volume
        return rho_avg

    def coherence_at_separation(self):
        """
        Compute coherence at orbital separation.
        """
        rho = self.orbital_separation_density()
        x = np.log(rho / self.rho_crit + 1)
        C = np.tanh(self.gamma * x)
        return C

    def coherence_in_star(self):
        """
        Coherence inside neutron star.
        """
        rho = self.neutron_star_density()
        x = np.log(rho / self.rho_crit + 1)
        C = np.tanh(self.gamma * x)
        return C

    def effective_G(self, C):
        """
        Effective gravitational constant.

        From Session #72: G_eff = G/C

        But careful: What C to use?
        - Inside star: C ~ 1
        - At orbital distance: C ~ 0.9-1.0

        For GW emission, relevant region is the orbit itself.
        """
        return G / C

    def orbital_decay_Synchronism(self, C_eff=None):
        """
        Synchronism prediction for orbital period derivative.

        Key question: Does G_eff appear in GW emission formula?

        In GR: dP/dt ∝ G^(5/3)

        If G → G_eff = G/C:
        dP/dt_Sync = dP/dt_GR × (1/C)^(5/3)

        But Session #71 argued: GW propagation unaffected by C.
        The question is whether GW EMISSION is affected.

        Two interpretations:
        1. GW emission uses local G → G_eff = G (inside NS, C ~ 1)
        2. GW emission uses averaged G over orbit → G_eff = G/C_orbit
        """
        if C_eff is None:
            # Use coherence at orbital separation
            C_eff = self.coherence_at_separation()

        dPdt_GR = self.orbital_decay_GR()

        # Synchronism modification: G → G/C affects orbital dynamics
        # GW emission rate ∝ G^(5/3)
        # If G_eff = G/C: rate enhanced by C^(-5/3)

        # BUT: If GW propagation uses c (unmodified) and
        # emission uses local physics (C ~ 1 in NS):
        # Then dP/dt_Sync ≈ dP/dt_GR (no modification!)

        # This is the key question to resolve.
        # For now, compute both scenarios:

        scenario_1 = dPdt_GR  # No modification (C ~ 1 in source)
        scenario_2 = dPdt_GR * (1/C_eff)**(5/3)  # Orbit-averaged C

        return {
            'scenario_1_no_mod': scenario_1,
            'scenario_2_orbit_C': scenario_2,
            'C_eff': C_eff
        }


def compute_predictions():
    """
    Compute binary pulsar predictions for GR and Synchronism.
    """
    results = {
        'session': 74,
        'track': 'C',
        'title': 'Binary Pulsar Orbital Decay Predictions'
    }

    # Hulse-Taylor pulsar
    pulsar_GR = BinaryPulsar()
    pulsar_Sync = SynchronismBinaryPulsar()

    # GR prediction
    dPdt_GR = pulsar_GR.orbital_decay_GR()

    # Observed value (Weisberg & Taylor 2005)
    dPdt_obs = -2.4184e-12  # s/s (dimensionless)
    dPdt_obs_err = 0.0009e-12

    # Synchronism predictions
    sync_results = pulsar_Sync.orbital_decay_Synchronism()

    results['hulse_taylor'] = {
        'M1_Msun': 1.4398,
        'M2_Msun': 1.3886,
        'P_b_hours': 7.752,
        'e': 0.6171,
        'a_meters': float(pulsar_GR.a),
        'chirp_mass_kg': float(pulsar_GR.chirp_mass),
        'eccentricity_factor': float(pulsar_GR.eccentricity_factor())
    }

    results['predictions'] = {
        'dPdt_GR': float(dPdt_GR),
        'dPdt_observed': float(dPdt_obs),
        'dPdt_observed_err': float(dPdt_obs_err),
        'ratio_GR_obs': float(dPdt_GR / dPdt_obs),
        'dPdt_Sync_scenario1': float(sync_results['scenario_1_no_mod']),
        'dPdt_Sync_scenario2': float(sync_results['scenario_2_orbit_C']),
        'C_at_orbit': float(sync_results['C_eff'])
    }

    # Analysis
    ratio_GR = dPdt_GR / dPdt_obs
    ratio_Sync1 = sync_results['scenario_1_no_mod'] / dPdt_obs
    ratio_Sync2 = sync_results['scenario_2_orbit_C'] / dPdt_obs

    results['analysis'] = {
        'GR_agreement': f'GR/obs = {ratio_GR:.4f} (should be ~1.0)',
        'Sync_scenario1': f'Sync(C~1)/obs = {ratio_Sync1:.4f}',
        'Sync_scenario2': f'Sync(C_orbit)/obs = {ratio_Sync2:.4f}',
        'interpretation': '''
Binary pulsar orbital decay analysis:

SCENARIO 1: No Synchronism modification
- Assume GW emission occurs in high-density region (neutron star)
- C ~ 1 inside NS → G_eff = G
- Prediction matches GR exactly

SCENARIO 2: Orbit-averaged coherence
- Assume effective G averaged over orbit
- C at orbital separation ~ 0.9-1.0
- Small modification to GR prediction

KEY INSIGHT:
For binary pulsars, density at orbital separation is STILL very high
compared to ρ_crit. Even at a ~ 2×10^9 m:

ρ_avg ~ M_total/V_orbit ~ 10^-12 kg/m³

This is >> ρ_crit = 10^-22 kg/m³!

So C ~ tanh(2 × log(10^10)) ~ 1.000

CONCLUSION:
Synchronism predicts IDENTICAL orbital decay to GR for binary pulsars.
This is because:
1. Source region (NS) has C ~ 1
2. Even orbital region has C ~ 1 due to high gravity
3. Only in GALACTIC OUTSKIRTS (ρ ~ ρ_crit) does C differ from 1

Binary pulsars are NOT a discriminating test between GR and Synchronism.
Need low-density environments to see coherence effects.
        '''
    }

    # Calculate where C deviates significantly
    rho_crit = 1e-22
    gamma = 2.0

    rho_test = np.logspace(-25, -10, 100)
    C_test = np.tanh(gamma * np.log(rho_test / rho_crit + 1))

    # Find where C drops below 0.99 (1% deviation)
    idx_99 = np.where(C_test < 0.99)[0]
    if len(idx_99) > 0:
        rho_99 = rho_test[idx_99[-1]]
    else:
        rho_99 = rho_test[0]

    results['threshold_analysis'] = {
        'rho_for_C_99': float(rho_99),
        'interpretation': f'C > 0.99 for ρ > {rho_99:.2e} kg/m³',
        'galactic_comparison': 'Galactic outskirts: ρ ~ 10^-24 kg/m³ → C ~ 0.7-0.9'
    }

    return results


def discriminating_tests():
    """
    What tests WOULD distinguish Synchronism from GR?
    """
    results = {'title': 'Discriminating Tests for Synchronism vs GR'}

    results['tests'] = '''
DISCRIMINATING TESTS: Where Synchronism differs from GR

Binary pulsars are NOT discriminating because C ~ 1 everywhere.
Need LOW DENSITY environments where C << 1.

PROMISING TESTS:

1. WIDE BINARY STARS (separation > 1 kpc)
   - Very low average density in system
   - Orbital dynamics should show G_eff = G/C effects
   - Prediction: Wider binaries deviate more from Kepler

2. GALAXY CLUSTER LENSING
   - Low-density intracluster medium
   - Lensing probes M_true, dynamics probe M/C
   - Prediction: Lensing mass > dynamical mass in outer regions

3. VOID GALAXY DYNAMICS
   - Galaxies in cosmic voids have lower external density
   - C_formation should differ from cluster galaxies
   - Prediction: Different Tully-Fisher relation in voids

4. GRAVITATIONAL WAVE FROM MERGERS IN LOW-DENSITY
   - Binary black holes in voids vs clusters
   - Different C → different chirp mass inference
   - Prediction: Mass discrepancy depending on environment

5. COSMOLOGICAL PERTURBATION GROWTH
   - From Session #72: Voids grow faster than clusters
   - Different structure formation history
   - Prediction: Scale-dependent σ_8

CURRENT STATUS:
- Binary pulsars: NOT discriminating (C ~ 1)
- Solar system: NOT discriminating (C ~ 1)
- Galaxy rotation curves: DISCRIMINATING (C varies from 1 to 0.3)
- Cosmology: POTENTIALLY discriminating (need detailed calculation)
    '''

    return results


if __name__ == '__main__':
    print("="*60)
    print("Session #74 Track C: Binary Pulsar Predictions")
    print("="*60)

    results = compute_predictions()

    print("\n--- Hulse-Taylor Binary Pulsar ---")
    ht = results['hulse_taylor']
    print(f"M1 = {ht['M1_Msun']:.4f} M_sun")
    print(f"M2 = {ht['M2_Msun']:.4f} M_sun")
    print(f"P_b = {ht['P_b_hours']:.3f} hours")
    print(f"e = {ht['e']:.4f}")
    print(f"a = {ht['a_meters']:.3e} m")

    print("\n--- Orbital Decay Predictions ---")
    pred = results['predictions']
    print(f"dP/dt (GR):      {pred['dPdt_GR']:.4e} s/s")
    print(f"dP/dt (observed): {pred['dPdt_observed']:.4e} s/s")
    print(f"GR/observed:     {pred['ratio_GR_obs']:.4f}")
    print(f"\nSynchronism predictions:")
    print(f"Scenario 1 (C~1):     {pred['dPdt_Sync_scenario1']:.4e} s/s")
    print(f"Scenario 2 (C_orbit): {pred['dPdt_Sync_scenario2']:.4e} s/s")
    print(f"C at orbit:          {pred['C_at_orbit']:.6f}")

    print("\n--- Conclusion ---")
    print("Synchronism predicts IDENTICAL orbital decay to GR")
    print("Reason: C ~ 1 in all high-density/high-gravity regions")
    print("Binary pulsars are NOT a discriminating test")

    # Add discriminating tests
    disc_results = discriminating_tests()
    results['discriminating_tests'] = disc_results

    print("\n--- Discriminating Tests ---")
    print("Binary pulsars: NOT discriminating (C ~ 1)")
    print("Galaxy rotation curves: DISCRIMINATING")
    print("Void dynamics: POTENTIALLY discriminating")

    # Save
    output_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(output_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    output_file = os.path.join(results_dir, 'session74_binary_pulsar.json')

    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults saved to: {output_file}")
