#!/usr/bin/env python3
"""
Synchronism: Code Supplementary Material for arXiv Submission
=============================================================

This module provides the core functions and validation routines used in:
"Synchronism: Dark Matter Phenomenology from Quantum Coherence in Galactic Systems"

Contents:
---------
1. Core Synchronism Functions (coherence, critical density)
2. Galaxy Validation Routines (DM fraction prediction)
3. Star Cluster Analysis
4. Galaxy Cluster Analysis with ICM Correction
5. Parameter Sensitivity Analysis

Usage:
------
    from synchronism_validation_code import SynchronismModel

    model = SynchronismModel()
    f_DM = model.predict_dm_fraction(rho=0.01, V=200)

Parameters:
-----------
- A = 0.028 M_sun/pc^3 (critical density coefficient)
- B = 0.5 (velocity exponent, derived from R ∝ V^0.75)
- gamma = 2.0 (decoherence exponent, derived from Γ ∝ ΔE^2)

Author: Synchronism Research Team
Date: November 2025
License: MIT
"""

import numpy as np
from typing import Tuple, Dict, List, Optional
from dataclasses import dataclass

# ==============================================================================
# PHYSICAL CONSTANTS
# ==============================================================================

G = 4.30e-6  # kpc³/(M_sun × Gyr²) - Gravitational constant
k_B = 1.38e-23  # J/K - Boltzmann constant
m_p = 1.67e-27  # kg - Proton mass
m_e = 9.11e-31  # kg - Electron mass
e = 1.6e-19  # C - Elementary charge
epsilon_0 = 8.85e-12  # F/m - Vacuum permittivity
M_sun = 1.989e30  # kg - Solar mass

# ==============================================================================
# SYNCHRONISM MODEL CLASS
# ==============================================================================

@dataclass
class SynchronismParams:
    """Parameters for the Synchronism model."""
    A: float = 0.028  # M_sun/pc³ - Critical density coefficient
    B: float = 0.5    # Velocity exponent
    gamma: float = 2.0  # Decoherence exponent


class SynchronismModel:
    """
    Core Synchronism model for predicting dark matter fractions.

    The coherence function is:
        C = tanh(γ × log(ρ/ρ_crit + 1))

    Where:
        ρ_crit = A × V^B

    The dark matter fraction is:
        f_DM = 1 - C
    """

    def __init__(self, params: Optional[SynchronismParams] = None):
        """Initialize with default or custom parameters."""
        self.params = params or SynchronismParams()

    def critical_density(self, V: float) -> float:
        """
        Calculate critical density for coherence transition.

        Parameters
        ----------
        V : float
            Velocity dispersion or circular velocity (km/s)

        Returns
        -------
        rho_crit : float
            Critical density (M_sun/pc³)
        """
        return self.params.A * V**self.params.B

    def coherence(self, rho: float, rho_crit: float) -> float:
        """
        Calculate coherence from density ratio.

        Parameters
        ----------
        rho : float
            Local baryon density (M_sun/pc³)
        rho_crit : float
            Critical density (M_sun/pc³)

        Returns
        -------
        C : float
            Coherence (0 to 1)
        """
        ratio = rho / rho_crit
        return np.tanh(self.params.gamma * np.log(ratio + 1))

    def predict_dm_fraction(self, rho: float, V: float) -> Tuple[float, float, float]:
        """
        Predict dark matter fraction for a given density and velocity.

        Parameters
        ----------
        rho : float
            Mean baryon density (M_sun/pc³)
        V : float
            Velocity dispersion or circular velocity (km/s)

        Returns
        -------
        f_DM : float
            Predicted dark matter fraction
        C : float
            Coherence
        rho_crit : float
            Critical density
        """
        rho_crit = self.critical_density(V)
        C = self.coherence(rho, rho_crit)
        f_DM = 1 - C
        return f_DM, C, rho_crit

    def predict_from_galaxy_params(
        self,
        M_star: float,
        R_eff: float,
        V: float
    ) -> Tuple[float, float, Dict]:
        """
        Predict DM fraction from observable galaxy parameters.

        Parameters
        ----------
        M_star : float
            Stellar mass (M_sun)
        R_eff : float
            Effective radius (kpc)
        V : float
            Velocity dispersion (km/s)

        Returns
        -------
        f_DM : float
            Predicted dark matter fraction
        C : float
            Coherence
        info : dict
            Additional calculated quantities
        """
        # Convert to mean density
        R_pc = R_eff * 1000  # kpc to pc
        volume = (4/3) * np.pi * R_pc**3
        rho_mean = M_star / volume

        # Calculate prediction
        f_DM, C, rho_crit = self.predict_dm_fraction(rho_mean, V)

        info = {
            'rho_mean': rho_mean,
            'rho_crit': rho_crit,
            'rho_ratio': rho_mean / rho_crit,
            'regime': self._classify_regime(rho_mean / rho_crit)
        }

        return f_DM, C, info

    @staticmethod
    def _classify_regime(rho_ratio: float) -> str:
        """Classify density regime."""
        if rho_ratio > 100:
            return "high_density"
        elif rho_ratio > 0.1:
            return "transition"
        else:
            return "low_density"


# ==============================================================================
# STAR CLUSTER VALIDATION
# ==============================================================================

def validate_star_cluster(
    M: float,
    R_pc: float,
    sigma_v: float,
    model: Optional[SynchronismModel] = None
) -> Dict:
    """
    Validate Synchronism prediction for a star cluster.

    Star clusters are expected to have f_DM ≈ 0 (no dark matter)
    because they are in the high-density regime.

    Parameters
    ----------
    M : float
        Total mass (M_sun)
    R_pc : float
        Half-light radius (pc)
    sigma_v : float
        Velocity dispersion (km/s)
    model : SynchronismModel, optional
        Model instance (uses default if not provided)

    Returns
    -------
    result : dict
        Validation results including prediction and regime
    """
    model = model or SynchronismModel()

    # Calculate mean density
    volume = (4/3) * np.pi * R_pc**3
    rho_mean = M / volume

    # Get prediction
    f_DM, C, rho_crit = model.predict_dm_fraction(rho_mean, sigma_v)

    # Expected result: f_DM ≈ 0 for all star clusters
    success = f_DM < 0.05

    return {
        'M': M,
        'R_pc': R_pc,
        'sigma_v': sigma_v,
        'rho_mean': rho_mean,
        'rho_crit': rho_crit,
        'rho_ratio': rho_mean / rho_crit,
        'C': C,
        'f_DM_pred': f_DM,
        'f_DM_obs': 0.0,  # Star clusters have no DM
        'success': success,
        'regime': 'high_density' if rho_mean / rho_crit > 10 else 'transition'
    }


# ==============================================================================
# GALAXY CLUSTER VALIDATION WITH ICM CORRECTION
# ==============================================================================

def icm_coherence(
    T: float,
    n_e: float,
    B_uG: float,
    R_kpc: float
) -> Tuple[float, Dict]:
    """
    Calculate ICM coherence from plasma physics.

    The ICM maintains partial coherence due to:
    1. Magnetic confinement (Larmor radius << cluster size)
    2. Plasma collective effects (N_D >> 1)
    3. Thermal coupling (sound crossing < Hubble time)

    Parameters
    ----------
    T : float
        ICM temperature (K)
    n_e : float
        Electron density (cm^-3)
    B_uG : float
        Magnetic field (microGauss)
    R_kpc : float
        Cluster radius (kpc)

    Returns
    -------
    C_icm : float
        ICM coherence (typically ~0.97)
    mechanisms : dict
        Individual mechanism contributions
    """
    # Convert units
    n_e_m3 = n_e * 1e6
    B_T = B_uG * 1e-10
    R_m = R_kpc * 3.086e19

    mechanisms = {}

    # 1. Magnetic coherence (Larmor radius)
    L_larmor = np.sqrt(2 * k_B * T * m_e) / (e * B_T) if B_T > 0 else 1e6
    L_larmor_kpc = L_larmor / 3.086e19
    C_magnetic = np.tanh(R_kpc / L_larmor_kpc / 1e6)
    mechanisms['magnetic'] = C_magnetic

    # 2. Plasma parameter (Debye sphere)
    lambda_D = np.sqrt(epsilon_0 * k_B * T / (n_e_m3 * e**2)) if n_e_m3 > 0 else 0
    N_D = (4/3) * np.pi * n_e_m3 * lambda_D**3 if lambda_D > 0 else 0
    C_plasma = np.tanh(np.log10(N_D + 1) / 10) if N_D > 0 else 0
    mechanisms['plasma'] = C_plasma

    # 3. Thermal coupling (sound crossing)
    c_s = np.sqrt(5/3 * k_B * T / m_p)
    c_s_kpc_Gyr = c_s * 3.156e16 / 3.086e19
    tau_sound = R_kpc / c_s_kpc_Gyr if c_s_kpc_Gyr > 0 else 1e10
    C_thermal = np.tanh(10 / tau_sound) if tau_sound > 0 else 0
    mechanisms['thermal'] = C_thermal

    # Combined coherence (geometric mean)
    C_icm = (C_magnetic * C_plasma * C_thermal)**(1/3)

    return C_icm, mechanisms


def validate_galaxy_cluster(
    M: float,
    R_200: float,
    sigma_v: float,
    f_icm: float,
    T_icm: float,
    n_e: float,
    B_uG: float,
    f_DM_obs: float,
    model: Optional[SynchronismModel] = None
) -> Dict:
    """
    Validate Synchronism prediction for a galaxy cluster with ICM correction.

    Parameters
    ----------
    M : float
        Total mass (M_sun)
    R_200 : float
        Virial radius (kpc)
    sigma_v : float
        Velocity dispersion (km/s)
    f_icm : float
        ICM mass fraction (typically 0.10-0.15)
    T_icm : float
        ICM temperature (K)
    n_e : float
        Electron density (cm^-3)
    B_uG : float
        Magnetic field (microGauss)
    f_DM_obs : float
        Observed dark matter fraction
    model : SynchronismModel, optional
        Model instance

    Returns
    -------
    result : dict
        Validation results with and without ICM correction
    """
    model = model or SynchronismModel()

    # Mean cluster density
    R_pc = R_200 * 1000
    volume = (4/3) * np.pi * R_pc**3
    rho_mean = M / volume

    # Base prediction (without ICM)
    f_DM_orig, C_base, rho_crit = model.predict_dm_fraction(rho_mean, sigma_v)

    # ICM coherence correction
    C_icm, mechanisms = icm_coherence(T_icm, n_e, B_uG, R_200)
    C_effective = f_icm * C_icm
    f_DM_corrected = 1.0 - C_effective

    return {
        'M': M,
        'R_200': R_200,
        'rho_mean': rho_mean,
        'rho_crit': rho_crit,
        'f_DM_obs': f_DM_obs,
        'f_DM_orig': f_DM_orig,
        'f_DM_corrected': f_DM_corrected,
        'C_icm': C_icm,
        'C_effective': C_effective,
        'error_orig': abs(f_DM_orig - f_DM_obs),
        'error_corrected': abs(f_DM_corrected - f_DM_obs),
        'mechanisms': mechanisms
    }


# ==============================================================================
# PARAMETER SENSITIVITY ANALYSIS
# ==============================================================================

def parameter_sensitivity(
    rho: float,
    V: float,
    perturbation: float = 0.20
) -> Dict:
    """
    Analyze sensitivity of prediction to parameter changes.

    Parameters
    ----------
    rho : float
        Baryon density (M_sun/pc³)
    V : float
        Velocity (km/s)
    perturbation : float
        Fractional perturbation (default 0.20 = ±20%)

    Returns
    -------
    sensitivity : dict
        Sensitivity to each parameter (A, B, gamma)
    """
    base_model = SynchronismModel()
    base_f_DM, _, _ = base_model.predict_dm_fraction(rho, V)

    sensitivity = {}

    # Vary A
    A_values = [0.028 * (1 - perturbation), 0.028, 0.028 * (1 + perturbation)]
    A_predictions = []
    for A in A_values:
        model = SynchronismModel(SynchronismParams(A=A))
        f_DM, _, _ = model.predict_dm_fraction(rho, V)
        A_predictions.append(f_DM)
    sensitivity['A'] = {
        'values': A_values,
        'predictions': A_predictions,
        'range': max(A_predictions) - min(A_predictions)
    }

    # Vary B
    B_values = [0.5 * (1 - perturbation), 0.5, 0.5 * (1 + perturbation)]
    B_predictions = []
    for B in B_values:
        model = SynchronismModel(SynchronismParams(B=B))
        f_DM, _, _ = model.predict_dm_fraction(rho, V)
        B_predictions.append(f_DM)
    sensitivity['B'] = {
        'values': B_values,
        'predictions': B_predictions,
        'range': max(B_predictions) - min(B_predictions)
    }

    # Vary gamma
    gamma_values = [2.0 * (1 - perturbation), 2.0, 2.0 * (1 + perturbation)]
    gamma_predictions = []
    for gamma in gamma_values:
        model = SynchronismModel(SynchronismParams(gamma=gamma))
        f_DM, _, _ = model.predict_dm_fraction(rho, V)
        gamma_predictions.append(f_DM)
    sensitivity['gamma'] = {
        'values': gamma_values,
        'predictions': gamma_predictions,
        'range': max(gamma_predictions) - min(gamma_predictions)
    }

    return sensitivity


# ==============================================================================
# EXAMPLE VALIDATION
# ==============================================================================

if __name__ == "__main__":
    print("="*70)
    print("SYNCHRONISM VALIDATION CODE - DEMONSTRATION")
    print("="*70)

    model = SynchronismModel()

    # Test 1: Dwarf galaxy (DM-dominated)
    print("\n1. DWARF GALAXY (DM-dominated)")
    f_DM, C, info = model.predict_from_galaxy_params(
        M_star=1e8,  # M_sun
        R_eff=1.0,   # kpc
        V=40         # km/s
    )
    print(f"   f_DM = {f_DM:.4f}, C = {C:.4f}")
    print(f"   ρ/ρ_crit = {info['rho_ratio']:.2e}")
    print(f"   Regime: {info['regime']}")

    # Test 2: Globular cluster (baryon-dominated)
    print("\n2. GLOBULAR CLUSTER (baryon-dominated)")
    result = validate_star_cluster(
        M=1e6,       # M_sun
        R_pc=5,      # pc
        sigma_v=10   # km/s
    )
    print(f"   f_DM = {result['f_DM_pred']:.4f}")
    print(f"   ρ/ρ_crit = {result['rho_ratio']:.2e}")
    print(f"   Success: {result['success']}")

    # Test 3: Galaxy cluster with ICM
    print("\n3. GALAXY CLUSTER (with ICM correction)")
    result = validate_galaxy_cluster(
        M=1e15,           # M_sun
        R_200=2000,       # kpc
        sigma_v=1000,     # km/s
        f_icm=0.15,
        T_icm=8.5e7,      # K
        n_e=3e-3,         # cm^-3
        B_uG=5,           # μG
        f_DM_obs=0.87
    )
    print(f"   f_DM_orig = {result['f_DM_orig']:.4f}")
    print(f"   f_DM_corrected = {result['f_DM_corrected']:.4f}")
    print(f"   f_DM_obs = {result['f_DM_obs']:.2f}")
    print(f"   Error reduction: {(1 - result['error_corrected']/result['error_orig'])*100:.1f}%")

    # Test 4: Parameter sensitivity
    print("\n4. PARAMETER SENSITIVITY (transition regime)")
    sens = parameter_sensitivity(rho=0.5, V=350)
    print(f"   Δf_DM(A): {sens['A']['range']:.4f}")
    print(f"   Δf_DM(B): {sens['B']['range']:.4f}")
    print(f"   Δf_DM(γ): {sens['gamma']['range']:.4f}")
    print(f"   Most sensitive: {'B' if sens['B']['range'] > max(sens['A']['range'], sens['gamma']['range']) else 'A or gamma'}")

    print("\n" + "="*70)
    print("DEMONSTRATION COMPLETE")
    print("="*70)
