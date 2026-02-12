#!/usr/bin/env python3
"""
MOND Offset Predictor — Standalone Tool
========================================

Predicts the Radial Acceleration Relation (RAR) offset for disk galaxies
using 3 observable inputs: V_flat, luminosity, and gas fraction.

The offset measures the deviation from the standard MOND prediction
(McGaugh 2016 interpolation function with assumed M/L = 0.5 for disks):

    offset = log10(g_obs) - log10(g_bar * nu(g_bar/a0))

This is primarily an M/L correction: galaxies with higher true M/L
have positive offsets, gas-rich galaxies have offsets closer to zero.

Models (from Sessions #483-587 of the Synchronism research program):
    3-var: offset = -3.238 + 1.739*logV - 0.450*logL - 0.374*f_gas
    6-var: offset = -3.379 + 1.897*logV - 0.548*logL - 0.218*c_V
                    - 0.451*f_gas + 0.147*logV*c_V + 0.181*logL*f_gas

Performance on SPARC (175 galaxies, Lelli+ 2016):
    3-var: LOO R^2 = 0.854, RMS = 0.060 dex (4 free parameters)
    6-var: LOO R^2 = 0.885, RMS = 0.053 dex (7 free parameters)

Physical basis (Session #587):
    - V-L ratio = 3.96 at mean galaxy (MOND predicts 4.0)
    - All coefficients derivable from MOND (V^4 = G*M_bar*a0) + M/L physics
    - Implied mean stellar M/L = 0.35 (literature: 0.44-0.50 at 3.6 um)
    - Gas fraction suppresses offset because gas M/L is precisely known

Usage:
    from mond_offset_predictor import predict_offset, predict_corrected_rar

    # Single galaxy
    result = predict_offset(vflat=120.0, luminosity=1e10, f_gas=0.15)
    print(result['offset'], result['uncertainty'])

    # Corrected RAR for a rotation curve
    g_obs = v_obs**2 / radius  # observed acceleration
    g_bar = (v_disk**2 + v_gas**2) / radius  # baryonic acceleration
    g_corrected = predict_corrected_rar(g_obs, g_bar, vflat=120.0,
                                         luminosity=1e10, f_gas=0.15)

Reference:
    Synchronism Research Program, Sessions #376-587
    SPARC database: Lelli, McGaugh & Schombert (2016, AJ 152, 157)
    MOND interpolation: McGaugh, Lelli & Schombert (2016, PRL 117, 201101)
"""

import numpy as np

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

A0_MOND = 1.2e-10    # MOND acceleration scale [m/s^2]
KPC_TO_M = 3.086e19  # kpc to meters
KMS_TO_MS = 1e3       # km/s to m/s

# ============================================================================
# MODEL COEFFICIENTS (fitted on SPARC, Sessions #483-585)
# ============================================================================

# 3-variable model (Session #585): offset = b0 + b1*logV + b2*logL + b3*f_gas
COEFF_3VAR = {
    'intercept': -3.238,
    'logV': 1.739,
    'logL': -0.450,
    'f_gas': -0.374,
}

# 6-variable model (Sessions #483-588): includes RC shape and interactions
# Coefficients fitted on the same 135-galaxy sample as the 3-var model
COEFF_6VAR = {
    'intercept': -3.097,
    'logV': 1.647,
    'logL': -0.486,
    'c_V': -0.332,
    'f_gas': -0.308,
    'logV_x_cV': 0.201,
    'logL_x_fgas': 0.163,
}

# LOO performance (calibration uncertainties)
LOO_RMS_3VAR = 0.060   # dex
LOO_RMS_6VAR = 0.053   # dex

# SPARC training sample statistics (for z-score outlier detection)
# NOTE: logL is in SPARC units = log10(luminosity / 10^9 L_sun)
# NOT in log10(L_sun). The previous value of 9.259 was a bug.
SPARC_STATS = {
    'logV_mean': 2.006, 'logV_std': 0.297,
    'logL_mean': 0.259, 'logL_std': 1.091,   # SPARC units: 10^9 L_sun
    'f_gas_mean': 0.184, 'f_gas_std': 0.196,
    'n_galaxies': 135,
}


# ============================================================================
# MOND INTERPOLATION
# ============================================================================

def nu_mcgaugh(x):
    """McGaugh (2016) MOND interpolation function.

    nu(x) = 1 / (1 - exp(-sqrt(x)))
    where x = g_bar / a0.

    Args:
        x: dimensionless acceleration g_bar / a0 (scalar or array)

    Returns:
        nu(x): the MOND boost factor g_obs/g_bar
    """
    x = np.asarray(x, dtype=float)
    x = np.clip(x, 1e-10, None)
    return 1.0 / (1.0 - np.exp(-np.sqrt(x)))


# ============================================================================
# PREDICTION FUNCTIONS
# ============================================================================

def predict_offset(vflat, luminosity, f_gas, c_V=None, model='3var'):
    """Predict the RAR offset for a galaxy.

    The offset is defined as:
        offset = log10(g_obs) - log10(g_bar * nu(g_bar/a0))
    and primarily reflects the M/L correction needed.

    Args:
        vflat: flat rotation velocity [km/s]
        luminosity: total luminosity at 3.6 um [SPARC units, 10^9 L_sun]
        f_gas: gas fraction = M_gas / M_baryonic (0 to 1)
        c_V: rotation curve concentration = <V_inner> / <V_outer>
             (required for model='6var', ignored for '3var')
        model: '3var' (default) or '6var'

    Returns:
        dict with keys:
            'offset': predicted offset [dex]
            'uncertainty': LOO RMS [dex] (calibration uncertainty)
            'implied_ml': implied stellar M/L at 3.6 um
            'extrapolation_warning': True if galaxy is outside training range
            'model': which model was used
    """
    logV = np.log10(max(vflat, 1e-3))
    logL = np.log10(max(luminosity, 1e-6))
    fg = np.clip(f_gas, 0.0, 1.0)

    # Check extrapolation
    z_logV = abs(logV - SPARC_STATS['logV_mean']) / SPARC_STATS['logV_std']
    z_logL = abs(logL - SPARC_STATS['logL_mean']) / SPARC_STATS['logL_std']
    z_fg = abs(fg - SPARC_STATS['f_gas_mean']) / SPARC_STATS['f_gas_std']
    extrapolation = max(z_logV, z_logL, z_fg) > 3.0

    if model == '6var':
        if c_V is None:
            raise ValueError("c_V is required for the 6-var model")
        c = COEFF_6VAR
        offset = (c['intercept'] + c['logV'] * logV + c['logL'] * logL
                  + c['c_V'] * c_V + c['f_gas'] * fg
                  + c['logV_x_cV'] * logV * c_V
                  + c['logL_x_fgas'] * logL * fg)
        unc = LOO_RMS_6VAR
    else:
        c = COEFF_3VAR
        offset = c['intercept'] + c['logV'] * logV + c['logL'] * logL + c['f_gas'] * fg
        unc = LOO_RMS_3VAR

    # Scale uncertainty for extrapolation
    if extrapolation:
        unc *= 1.5

    # Implied M/L: offset ≈ log(Upsilon_true / 0.5)
    # So Upsilon_true ≈ 0.5 * 10^offset
    implied_ml = 0.5 * 10**offset

    return {
        'offset': offset,
        'uncertainty': unc,
        'implied_ml': implied_ml,
        'extrapolation_warning': extrapolation,
        'model': model,
    }


def predict_corrected_rar(g_obs, g_bar, vflat, luminosity, f_gas,
                           c_V=None, model='3var'):
    """Apply the M/L correction to a galaxy's RAR data points.

    Removes the predicted M/L offset from the observed RAR, producing
    a tighter radial acceleration relation.

    Args:
        g_obs: observed centripetal acceleration [m/s^2] (array)
        g_bar: baryonic acceleration [m/s^2] (array)
        vflat: flat rotation velocity [km/s]
        luminosity: total luminosity at 3.6 um [SPARC units, 10^9 L_sun]
        f_gas: gas fraction
        c_V: RC concentration (for 6-var model)
        model: '3var' or '6var'

    Returns:
        dict with keys:
            'log_g_obs_corrected': corrected log10(g_obs) [m/s^2]
            'log_g_bar': log10(g_bar) [m/s^2]
            'offset_applied': the offset that was subtracted [dex]
            'residual_scatter': expected remaining scatter [dex]
    """
    result = predict_offset(vflat, luminosity, f_gas, c_V=c_V, model=model)
    offset = result['offset']

    log_gobs = np.log10(np.clip(np.asarray(g_obs, dtype=float), 1e-15, None))
    log_gbar = np.log10(np.clip(np.asarray(g_bar, dtype=float), 1e-15, None))

    # The offset is added to log(g_bar * nu), so subtract it from log(g_obs)
    log_gobs_corrected = log_gobs - offset

    return {
        'log_g_obs_corrected': log_gobs_corrected,
        'log_g_bar': log_gbar,
        'offset_applied': offset,
        'residual_scatter': result['uncertainty'],
    }


def predict_batch(vflat_arr, luminosity_arr, f_gas_arr, model='3var'):
    """Predict offsets for multiple galaxies at once.

    Args:
        vflat_arr: array of V_flat values [km/s]
        luminosity_arr: array of luminosities [L_sun]
        f_gas_arr: array of gas fractions
        model: '3var' or '6var' (6var not supported in batch without c_V)

    Returns:
        dict with arrays:
            'offsets': predicted offsets [dex]
            'uncertainties': per-galaxy uncertainties [dex]
            'implied_ml': implied M/L values
    """
    vflat_arr = np.asarray(vflat_arr, dtype=float)
    luminosity_arr = np.asarray(luminosity_arr, dtype=float)
    f_gas_arr = np.asarray(f_gas_arr, dtype=float)

    logV = np.log10(np.clip(vflat_arr, 1e-3, None))
    logL = np.log10(np.clip(luminosity_arr, 1e-6, None))
    fg = np.clip(f_gas_arr, 0.0, 1.0)

    c = COEFF_3VAR
    offsets = c['intercept'] + c['logV'] * logV + c['logL'] * logL + c['f_gas'] * fg

    # Check extrapolation per galaxy
    z_max = np.maximum(
        np.abs(logV - SPARC_STATS['logV_mean']) / SPARC_STATS['logV_std'],
        np.maximum(
            np.abs(logL - SPARC_STATS['logL_mean']) / SPARC_STATS['logL_std'],
            np.abs(fg - SPARC_STATS['f_gas_mean']) / SPARC_STATS['f_gas_std']
        )
    )
    uncertainties = np.where(z_max > 3.0, LOO_RMS_3VAR * 1.5, LOO_RMS_3VAR)

    implied_ml = 0.5 * 10**offsets

    return {
        'offsets': offsets,
        'uncertainties': uncertainties,
        'implied_ml': implied_ml,
    }


# ============================================================================
# CONVENIENCE: FROM ROTATION CURVE DATA
# ============================================================================

def compute_galaxy_features(v_obs, v_gas, v_disk, v_bul, radius, vflat, luminosity):
    """Compute all model features from rotation curve data.

    This reproduces the exact pipeline used in SPARC analysis (Sessions #372-587).

    Args:
        v_obs: observed rotation velocity [km/s] (array, per radius)
        v_gas: gas contribution velocity [km/s] (array)
        v_disk: disk contribution velocity [km/s] (array, at M/L=1)
        v_bul: bulge contribution velocity [km/s] (array, at M/L=1)
        radius: galactocentric radius [kpc] (array)
        vflat: flat rotation velocity [km/s] (scalar)
        luminosity: total luminosity at 3.6 um [SPARC units, 10^9 L_sun] (scalar)

    Returns:
        dict with:
            'logV': log10(V_flat)
            'logL': log10(luminosity)
            'f_gas': gas fraction (kinematic)
            'c_V': RC concentration
            'offset_outer': measured outer offset [dex]
            'n_points': number of valid data points
    """
    v_obs = np.asarray(v_obs, dtype=float)
    v_gas = np.asarray(v_gas, dtype=float)
    v_disk = np.asarray(v_disk, dtype=float)
    v_bul = np.asarray(v_bul, dtype=float)
    radius = np.asarray(radius, dtype=float)

    valid = (v_obs > 0) & (radius > 0)
    v_obs = v_obs[valid]
    v_gas = v_gas[valid]
    v_disk = v_disk[valid]
    v_bul = v_bul[valid]
    radius = radius[valid]

    # Accelerations
    g_obs = (v_obs * KMS_TO_MS)**2 / (radius * KPC_TO_M)
    g_bar = (np.abs(v_disk * KMS_TO_MS)**2 / (radius * KPC_TO_M)
             + np.abs(v_gas * KMS_TO_MS)**2 / (radius * KPC_TO_M))
    if np.any(v_bul != 0):
        g_bar += np.abs(v_bul * KMS_TO_MS)**2 / (radius * KPC_TO_M)
    g_bar = np.clip(g_bar, 1e-15, None)

    # MOND offset
    x = g_bar / A0_MOND
    nu_val = nu_mcgaugh(x)
    offset_pts = np.log10(g_obs) - np.log10(g_bar * nu_val)

    # Outer region average
    r_frac = radius / np.max(radius)
    outer = r_frac > 0.5
    if outer.sum() < 2:
        outer = r_frac > 0.3
    offset_outer = np.mean(offset_pts[outer]) if outer.sum() >= 2 else np.mean(offset_pts)

    # RC concentration
    mid = len(v_obs) // 2
    c_V = np.mean(v_obs[:mid]) / np.mean(v_obs[mid:]) if np.mean(v_obs[mid:]) > 0 else 1.0

    # Gas fraction (kinematic)
    gas_m = np.sum(np.abs(v_gas)**2)
    tot_m = gas_m + np.sum(np.abs(v_disk)**2)
    if np.any(v_bul != 0):
        tot_m += np.sum(np.abs(v_bul)**2)
    f_gas = gas_m / tot_m if tot_m > 0 else 0.0

    return {
        'logV': np.log10(max(vflat, 1e-3)),
        'logL': np.log10(max(luminosity, 1e-6)),
        'f_gas': f_gas,
        'c_V': c_V,
        'offset_outer': offset_outer,
        'n_points': len(v_obs),
    }


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Predict MOND RAR offset for a disk galaxy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --vflat 120 --luminosity 1e10 --f_gas 0.15
  %(prog)s --vflat 50 --luminosity 1e8 --f_gas 0.6 --model 6var --c_V 0.7
        """)
    parser.add_argument('--vflat', type=float, required=True,
                        help='Flat rotation velocity [km/s]')
    parser.add_argument('--luminosity', type=float, required=True,
                        help='Total luminosity at 3.6 um [L_sun]')
    parser.add_argument('--f_gas', type=float, required=True,
                        help='Gas fraction M_gas/M_baryonic (0-1)')
    parser.add_argument('--c_V', type=float, default=None,
                        help='RC concentration <V_inner>/<V_outer> (for 6var)')
    parser.add_argument('--model', choices=['3var', '6var'], default='3var',
                        help='Model to use (default: 3var)')

    args = parser.parse_args()

    result = predict_offset(args.vflat, args.luminosity, args.f_gas,
                             c_V=args.c_V, model=args.model)

    print(f"Galaxy: V_flat={args.vflat} km/s, L={args.luminosity:.2e} [10^9 L_sun], f_gas={args.f_gas:.3f}")
    print(f"Model: {result['model']}")
    print(f"Predicted offset: {result['offset']:+.4f} dex")
    print(f"Uncertainty (LOO): {result['uncertainty']:.4f} dex")
    print(f"Implied M/L (3.6um): {result['implied_ml']:.3f}")
    if result['extrapolation_warning']:
        print("WARNING: Galaxy is outside the SPARC training range")
