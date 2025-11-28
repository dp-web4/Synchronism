#!/usr/bin/env python3
"""
Session #59: Gravitational Wave Coherence Analysis

Applies the Synchronism coherence framework to gravitational wave physics.
Analyzes GW170817 and makes predictions for future observations.

Key equations:
    C = tanh(γ × log(ρ/ρ_crit + 1))   [same as dark matter]
    c_g/c = 1 + α × (1 - C_avg)       [GW speed modification]

Parameters from dark matter validation (Sessions #52-58):
    γ = 2.0, A = 0.028 M_sun/pc³, B = 0.5
"""

import numpy as np
from typing import Dict, Tuple, List
from dataclasses import dataclass
import json

# =============================================================================
# CONSTANTS
# =============================================================================

c = 299792.458  # km/s - speed of light
G = 4.30e-6     # kpc³/(M_sun × Gyr²)
Mpc_to_kpc = 1000
Mpc_to_m = 3.086e22

# Dark matter coherence parameters (from arXiv paper)
GAMMA = 2.0
A = 0.028       # M_sun/pc³
B = 0.5

# =============================================================================
# COHERENCE FUNCTIONS (from dark matter model)
# =============================================================================

def coherence(rho: float, rho_crit: float, gamma: float = GAMMA) -> float:
    """
    Calculate coherence from density ratio.
    Same function validated on 195 astrophysical systems.

    Parameters
    ----------
    rho : float
        Local density (M_sun/pc³)
    rho_crit : float
        Critical density (M_sun/pc³)

    Returns
    -------
    C : float
        Coherence (0-1)
    """
    ratio = rho / rho_crit
    return np.tanh(gamma * np.log(ratio + 1))


def critical_density(V: float, A: float = A, B: float = B) -> float:
    """
    Calculate critical density for coherence.

    Parameters
    ----------
    V : float
        Characteristic velocity (km/s)

    Returns
    -------
    rho_crit : float
        Critical density (M_sun/pc³)
    """
    return A * V**B


# =============================================================================
# GW170817 ANALYSIS
# =============================================================================

@dataclass
class GWEvent:
    """Gravitational wave event parameters."""
    name: str
    distance_Mpc: float      # Luminosity distance
    distance_err: float      # Distance uncertainty
    gw_arrival_gps: float    # GPS time of GW arrival
    em_arrival_gps: float    # GPS time of EM arrival (if multi-messenger)
    em_delay_s: float        # Delay between GW and EM (seconds)
    host_mass_Msun: float    # Host galaxy stellar mass
    host_fdm: float          # Host galaxy DM fraction (estimated)
    ra_deg: float            # Right ascension
    dec_deg: float           # Declination


# GW170817 - First multi-messenger event
GW170817 = GWEvent(
    name="GW170817",
    distance_Mpc=40.0,           # 40 ± 8 Mpc
    distance_err=8.0,
    gw_arrival_gps=1187008882.4, # GPS time
    em_arrival_gps=1187008884.3, # GRB 170817A
    em_delay_s=1.74,             # 1.74 ± 0.05 seconds
    host_mass_Msun=3e10,         # NGC 4993 stellar mass
    host_fdm=0.3,                # Estimated from morphology
    ra_deg=197.450,
    dec_deg=-23.381
)


def estimate_los_coherence(event: GWEvent) -> Dict:
    """
    Estimate average coherence along line-of-sight.

    Uses simplified model:
    - Intergalactic medium: ρ ~ 10^-7 M_sun/pc³, C ~ 0
    - Galaxy intersections: Add high-C patches
    - Host galaxy: Final crossing

    Parameters
    ----------
    event : GWEvent
        The gravitational wave event

    Returns
    -------
    result : dict
        Line-of-sight coherence analysis
    """
    # Typical intergalactic density
    rho_igm = 1e-7  # M_sun/pc³ (very low)

    # Characteristic velocity for IGM (from Hubble flow)
    V_igm = 70 * event.distance_Mpc  # km/s (rough)
    V_igm = max(V_igm, 100)  # Minimum 100 km/s

    # Critical density for IGM
    rho_crit_igm = critical_density(V_igm)

    # IGM coherence (should be ~0)
    C_igm = coherence(rho_igm, rho_crit_igm)

    # Host galaxy contribution
    # Assume GW travels through ~10 kpc of host
    host_path_kpc = 10.0

    # Host density estimate
    rho_host = event.host_mass_Msun / (4/3 * np.pi * (10 * 1000)**3)  # M_sun/pc³
    V_host = 200  # km/s typical
    rho_crit_host = critical_density(V_host)
    C_host = coherence(rho_host, rho_crit_host)

    # Path lengths
    total_path_kpc = event.distance_Mpc * Mpc_to_kpc
    igm_fraction = (total_path_kpc - host_path_kpc) / total_path_kpc
    host_fraction = host_path_kpc / total_path_kpc

    # Weighted average coherence
    C_avg = igm_fraction * C_igm + host_fraction * C_host

    return {
        'event': event.name,
        'distance_Mpc': event.distance_Mpc,
        'total_path_kpc': total_path_kpc,
        'rho_igm': rho_igm,
        'rho_crit_igm': rho_crit_igm,
        'C_igm': float(C_igm),
        'rho_host': rho_host,
        'rho_crit_host': rho_crit_host,
        'C_host': float(C_host),
        'igm_fraction': igm_fraction,
        'host_fraction': host_fraction,
        'C_avg': float(C_avg),
        '1_minus_C_avg': float(1 - C_avg)
    }


def gw_speed_modification(C_avg: float, alpha: float) -> Dict:
    """
    Calculate GW speed modification from coherence.

    c_g/c = 1 + α × (1 - C_avg)

    Parameters
    ----------
    C_avg : float
        Average coherence along line-of-sight
    alpha : float
        Coupling parameter (<<1)

    Returns
    -------
    result : dict
        Speed modification analysis
    """
    delta_c_over_c = alpha * (1 - C_avg)
    c_g = c * (1 + delta_c_over_c)

    return {
        'alpha': alpha,
        'C_avg': C_avg,
        '1_minus_C': 1 - C_avg,
        'delta_c_over_c': delta_c_over_c,
        'c_g_km_s': c_g,
        'c_g_minus_c_km_s': c_g - c
    }


def gw170817_constraint() -> Dict:
    """
    Derive constraint on α from GW170817.

    GW170817 + GRB 170817A gave:
    |c_g - c|/c < 3 × 10^-15 (Abbott et al. 2017)

    Returns
    -------
    constraint : dict
        Constraint analysis
    """
    # Observed constraint
    delta_c_limit = 3e-15

    # Calculate line-of-sight coherence for GW170817
    los = estimate_los_coherence(GW170817)

    # Constraint on alpha
    one_minus_C = los['1_minus_C_avg']
    alpha_max = delta_c_limit / one_minus_C

    return {
        'event': 'GW170817',
        'observed_delta_c_limit': delta_c_limit,
        'C_avg': los['C_avg'],
        '1_minus_C_avg': one_minus_C,
        'alpha_max': alpha_max,
        'constraint': f'α < {alpha_max:.2e}'
    }


# =============================================================================
# PREDICTIONS FOR FUTURE OBSERVATIONS
# =============================================================================

def predict_correlation(events: List[GWEvent], alpha: float) -> Dict:
    """
    Predict correlation between GW delay and DM column.

    Parameters
    ----------
    events : list
        List of GW events
    alpha : float
        Coupling parameter

    Returns
    -------
    prediction : dict
        Predicted correlation
    """
    delays = []
    dm_columns = []

    for event in events:
        los = estimate_los_coherence(event)

        # Predicted time delay relative to EM
        # Δt/D × c = α × (1 - C_avg)
        delta_t_over_D = alpha * (1 - los['C_avg']) / c  # s/km

        # Convert to s/Mpc
        delta_t_per_Mpc = delta_t_over_D * Mpc_to_kpc * 1000  # s/Mpc

        delays.append(delta_t_per_Mpc * event.distance_Mpc)
        dm_columns.append(1 - los['C_avg'])

    # Correlation coefficient
    if len(events) > 1:
        correlation = np.corrcoef(delays, dm_columns)[0, 1]
    else:
        correlation = None

    return {
        'n_events': len(events),
        'alpha': alpha,
        'delays_s': delays,
        'dm_columns': dm_columns,
        'correlation': correlation
    }


def ringdown_prediction(event: GWEvent, delta: float = 1e-4) -> Dict:
    """
    Predict ringdown frequency modification in DM-rich environments.

    f_ringdown = f_GR × (1 + δ × f_DM,host)

    Parameters
    ----------
    event : GWEvent
        The gravitational wave event
    delta : float
        Coupling parameter for ringdown

    Returns
    -------
    prediction : dict
        Ringdown modification
    """
    # Typical ringdown frequency for BNS merger remnant
    # (If it forms a black hole before collapse)
    f_ringdown_GR = 2000  # Hz (rough)

    # Modification from host DM
    f_modification = delta * event.host_fdm
    f_ringdown_sync = f_ringdown_GR * (1 + f_modification)

    return {
        'event': event.name,
        'host_fdm': event.host_fdm,
        'delta': delta,
        'f_ringdown_GR_Hz': f_ringdown_GR,
        'f_modification': f_modification,
        'f_ringdown_sync_Hz': f_ringdown_sync,
        'delta_f_Hz': f_ringdown_sync - f_ringdown_GR
    }


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    """Run GW coherence analysis for Session #59."""

    print("=" * 70)
    print("SESSION #59: GRAVITATIONAL WAVE COHERENCE ANALYSIS")
    print("=" * 70)

    # 1. GW170817 Line-of-Sight Analysis
    print("\n1. GW170817 LINE-OF-SIGHT COHERENCE")
    print("-" * 50)

    los = estimate_los_coherence(GW170817)

    print(f"   Event: {los['event']}")
    print(f"   Distance: {los['distance_Mpc']:.1f} Mpc")
    print(f"   Total path: {los['total_path_kpc']:.0e} kpc")
    print(f"\n   IGM properties:")
    print(f"     ρ_IGM: {los['rho_igm']:.1e} M_sun/pc³")
    print(f"     ρ_crit: {los['rho_crit_igm']:.1e} M_sun/pc³")
    print(f"     C_IGM: {los['C_igm']:.4f}")
    print(f"\n   Host galaxy (NGC 4993):")
    print(f"     ρ_host: {los['rho_host']:.2e} M_sun/pc³")
    print(f"     ρ_crit: {los['rho_crit_host']:.4f} M_sun/pc³")
    print(f"     C_host: {los['C_host']:.4f}")
    print(f"\n   Weighted average:")
    print(f"     C_avg: {los['C_avg']:.4f}")
    print(f"     1 - C_avg: {los['1_minus_C_avg']:.4f}")

    # 2. Constraint from GW170817
    print("\n\n2. CONSTRAINT ON α FROM GW170817")
    print("-" * 50)

    constraint = gw170817_constraint()

    print(f"   Observed limit: |Δc/c| < {constraint['observed_delta_c_limit']:.0e}")
    print(f"   Line-of-sight 1-C: {constraint['1_minus_C_avg']:.4f}")
    print(f"\n   CONSTRAINT: {constraint['constraint']}")

    # 3. Speed modification at limit
    print("\n\n3. SPEED MODIFICATION AT CONSTRAINT LIMIT")
    print("-" * 50)

    alpha_limit = constraint['alpha_max']
    speed = gw_speed_modification(los['C_avg'], alpha_limit)

    print(f"   α = {speed['alpha']:.2e}")
    print(f"   Δc/c = {speed['delta_c_over_c']:.2e}")
    print(f"   c_g = {speed['c_g_km_s']:.6f} km/s")
    print(f"   c_g - c = {speed['c_g_minus_c_km_s']:.6e} km/s")

    # 4. Ringdown prediction
    print("\n\n4. RINGDOWN FREQUENCY PREDICTION")
    print("-" * 50)

    ringdown = ringdown_prediction(GW170817, delta=1e-4)

    print(f"   Host f_DM: {ringdown['host_fdm']:.2f}")
    print(f"   δ parameter: {ringdown['delta']:.0e}")
    print(f"\n   GR prediction: {ringdown['f_ringdown_GR_Hz']:.0f} Hz")
    print(f"   Synchronism prediction: {ringdown['f_ringdown_sync_Hz']:.1f} Hz")
    print(f"   Δf: {ringdown['delta_f_Hz']:.2f} Hz")

    # 5. Key predictions summary
    print("\n\n5. TESTABLE PREDICTIONS SUMMARY")
    print("-" * 50)

    predictions = {
        'GW_speed_vs_DM_column': {
            'prediction': 'c_g/c = 1 + α(1-C_avg)',
            'constraint': f'α < {alpha_limit:.2e}',
            'test': 'Correlation of Δt with DM column across multi-messenger events',
            'required_events': '20-50 for 3σ detection',
            'falsification': 'No correlation to 10^-16 level'
        },
        'ringdown_anomaly': {
            'prediction': 'f_ring = f_GR × (1 + δ × f_DM)',
            'estimate': 'δ ~ 10^-4 to 10^-5',
            'test': 'Compare ringdown in high-DM vs low-DM host galaxies',
            'required_events': '~100 with identified hosts',
            'falsification': 'Ringdown exactly matches GR independent of host'
        },
        'stochastic_background': {
            'prediction': 'SGWB anisotropic following DM distribution',
            'test': 'Correlate SGWB with large-scale structure',
            'instrument': 'LISA, Einstein Telescope',
            'falsification': 'SGWB isotropic to arcminute scales'
        }
    }

    for name, pred in predictions.items():
        print(f"\n   {name}:")
        for key, value in pred.items():
            print(f"     {key}: {value}")

    # 6. Connection to dark matter paper
    print("\n\n6. CONNECTION TO DARK MATTER VALIDATION")
    print("-" * 50)

    print("   Same coherence function C = tanh(γ log(ρ/ρ_crit + 1))")
    print("   Same parameters: γ=2.0, A=0.028, B=0.5")
    print(f"\n   Dark matter: f_DM = 1 - C (validated on 195 systems)")
    print(f"   GW speed:    c_g/c = 1 + α(1-C)")
    print("\n   UNIFIED PREDICTION: Regions with high f_DM show GW effects")

    # Save results
    results = {
        'los_analysis': los,
        'constraint': constraint,
        'speed_at_limit': speed,
        'ringdown': ringdown,
        'predictions': predictions
    }

    # Convert numpy types for JSON
    def convert_types(obj):
        if isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [convert_types(v) for v in obj]
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open('session59_gw_results.json', 'w') as f:
        json.dump(convert_types(results), f, indent=2)

    print("\n\n" + "=" * 70)
    print("ANALYSIS COMPLETE - Results saved to session59_gw_results.json")
    print("=" * 70)

    return results


if __name__ == "__main__":
    results = main()
