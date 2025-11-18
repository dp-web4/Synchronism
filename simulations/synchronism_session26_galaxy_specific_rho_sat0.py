#!/usr/bin/env python3
"""
Synchronism Research - Session #26
Galaxy-Specific ρ_sat,0 Model Test

Following Nova Session #25 recommendation (Priority 1): Test if ρ_sat,0 is
galaxy-specific rather than universal.

PARADOX (Session #25): NGC has stronger B-fields (12 µG) than F/DDO (4 µG),
yet lower ρ_sat. Magnetic screening model predicts INVERSE.

HYPOTHESIS: ρ_sat,0 varies with galaxy properties (morphology, SFR, dynamics)
- ρ_sat,0,NGC << ρ_sat,0,F explains NGC low ρ_sat despite strong B-field
- Physical mechanism: Synchronism coherence threshold depends on dynamics
  - Ordered rotation (NGC) → Low coherence threshold → Low ρ_sat,0
  - Chaotic motion (F) → High coherence threshold → High ρ_sat,0

TEST: Refit Session #22 magnetic screening model with galaxy-type-specific ρ_sat,0
"""

import numpy as np
import pickle
from pathlib import Path
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# ============================================================================
# PART 1: LOAD SESSION #20 DATA
# ============================================================================

def load_session20_data():
    """Load empirical ρ_sat data from Session #20."""

    data_path = Path(__file__).parent / "session20_rho_sat_data.pkl"

    if not data_path.exists():
        raise FileNotFoundError(f"Session #20 data not found: {data_path}")

    with open(data_path, 'rb') as f:
        data = pickle.load(f)

    print(f"✓ Loaded Session #20 data from: {data_path}")
    print(f"  {data['n_galaxies']} galaxies ({data['n_success']} successful fits)")
    print()

    # Rename for consistency
    data['rho_central'] = np.array(data['rho_centrals'])
    data['rho_sat'] = np.array(data['rho_sats_fitted'])

    return data


# ============================================================================
# PART 2: GALAXY TYPE CLASSIFICATION
# ============================================================================

def classify_galaxy_types(data):
    """
    Classify galaxies by type for galaxy-specific ρ_sat,0 fitting.

    SPARC naming conventions:
    - NGC: Spirals (typically)
    - UGC: Spirals
    - F: Irregulars
    - DDO: Dwarfs
    - IC: Mixed (irregulars/spirals)
    - M: Messier objects (often NGC)
    - Other: Various

    For this test, we simplify to two groups:
    - SPIRAL: NGC + UGC (ordered rotation, evolved systems)
    - IRREGULAR: F + DDO + other (chaotic motion, young systems)
    """

    print("CLASSIFYING GALAXIES BY TYPE")
    print("="*80)
    print()

    galaxy_names = data['galaxy_names']

    # Classification
    spiral_indices = []
    irregular_indices = []

    for i, name in enumerate(galaxy_names):
        if name.startswith('NGC') or name.startswith('UGC'):
            spiral_indices.append(i)
        else:  # F, DDO, IC, M, other
            irregular_indices.append(i)

    print(f"SPIRAL group (NGC + UGC): n = {len(spiral_indices)}")
    print(f"IRREGULAR group (F + DDO + other): n = {len(irregular_indices)}")
    print()

    # Create type arrays
    galaxy_types = np.array(['SPIRAL' if i in spiral_indices else 'IRREGULAR'
                             for i in range(len(galaxy_names))])

    return {
        'galaxy_types': galaxy_types,
        'spiral_indices': spiral_indices,
        'irregular_indices': irregular_indices
    }


# ============================================================================
# PART 3: UNIVERSAL MODEL (SESSION #22 BASELINE)
# ============================================================================

def fit_universal_model(data):
    """
    Refit Session #22 universal magnetic screening model as baseline.

    Model: ρ_sat = ρ_sat,0 / [1 + (ρ_central/ρ_mag)^δ]

    Universal parameters: ρ_sat,0, ρ_mag, δ (same for all galaxies)
    """

    print("="*80)
    print("UNIVERSAL MODEL (SESSION #22 BASELINE)")
    print("="*80)
    print()

    # Extract data
    rho_central = data['rho_central']
    rho_sat = data['rho_sat']

    # Filter valid data (positive values)
    valid = (rho_central > 0) & (rho_sat > 0) & np.isfinite(rho_central) & np.isfinite(rho_sat)
    rho_c = rho_central[valid]
    rho_s = rho_sat[valid]

    print(f"Valid galaxies: {np.sum(valid)}/{len(rho_central)}")
    print()

    # Define model
    def rho_sat_model(rho_c, rho_sat_0, rho_mag, delta):
        """Magnetic screening model."""
        return rho_sat_0 / (1.0 + (rho_c / rho_mag)**delta)

    # Initial guess
    p0 = [1e6, 300, 2.0]  # rho_sat_0, rho_mag, delta

    # Fit (log space for better convergence across orders of magnitude)
    def log_model(rho_c, log_rho_sat_0, log_rho_mag, delta):
        rho_sat_0 = 10**log_rho_sat_0
        rho_mag = 10**log_rho_mag
        return np.log10(rho_sat_model(rho_c, rho_sat_0, rho_mag, delta))

    # Convert to log space
    log_rho_s = np.log10(rho_s)
    p0_log = [np.log10(p0[0]), np.log10(p0[1]), p0[2]]

    # Fit
    try:
        popt, pcov = curve_fit(log_model, rho_c, log_rho_s, p0=p0_log,
                              bounds=([-np.inf, -np.inf, 0], [np.inf, np.inf, 10]),
                              maxfev=10000)

        # Convert back
        rho_sat_0 = 10**popt[0]
        rho_mag = 10**popt[1]
        delta = popt[2]

        # Uncertainties
        perr = np.sqrt(np.diag(pcov))
        rho_sat_0_err = rho_sat_0 * perr[0] * np.log(10)
        rho_mag_err = rho_mag * perr[1] * np.log(10)
        delta_err = perr[2]

        # Predictions
        rho_sat_pred = rho_sat_model(rho_c, rho_sat_0, rho_mag, delta)

        # R²
        ss_res = np.sum((np.log10(rho_s) - np.log10(rho_sat_pred))**2)
        ss_tot = np.sum((np.log10(rho_s) - np.mean(np.log10(rho_s)))**2)
        r_squared = 1 - (ss_res / ss_tot)

        # RMS in log space
        rms_log = np.sqrt(np.mean((np.log10(rho_s) - np.log10(rho_sat_pred))**2))

        # Correlation
        correlation = np.corrcoef(np.log10(rho_s), np.log10(rho_sat_pred))[0, 1]

        print("✅ FIT SUCCESSFUL")
        print()
        print("Best-fit parameters:")
        print(f"  ρ_sat,0 = {rho_sat_0:.2e} ± {rho_sat_0_err:.2e} M_☉/pc³")
        print(f"  ρ_mag   = {rho_mag:.2e} ± {rho_mag_err:.2e} M_☉/pc³")
        print(f"  δ       = {delta:.2f} ± {delta_err:.2f}")
        print()
        print("Fit quality:")
        print(f"  R² = {r_squared:.3f}")
        print(f"  RMS(log) = {rms_log:.3f}")
        print(f"  Correlation: r = {correlation:.3f}")
        print()

        return {
            'rho_sat_0': rho_sat_0,
            'rho_mag': rho_mag,
            'delta': delta,
            'rho_sat_0_err': rho_sat_0_err,
            'rho_mag_err': rho_mag_err,
            'delta_err': delta_err,
            'r_squared': r_squared,
            'rms_log': rms_log,
            'correlation': correlation,
            'rho_sat_pred': rho_sat_pred,
            'valid_mask': valid
        }

    except Exception as e:
        print(f"❌ FIT FAILED: {e}")
        return None


# ============================================================================
# PART 4: GALAXY-SPECIFIC ρ_sat,0 MODEL
# ============================================================================

def fit_galaxy_specific_rho_sat0(data, type_data):
    """
    Fit magnetic screening model with galaxy-type-specific ρ_sat,0.

    Model: ρ_sat = ρ_sat,0(type) / [1 + (ρ_central/ρ_mag)^δ]

    Parameters:
    - ρ_sat,0,SPIRAL: Coherence threshold for spiral galaxies
    - ρ_sat,0,IRREGULAR: Coherence threshold for irregulars
    - ρ_mag: Universal magnetic screening scale (same for all)
    - δ: Universal screening exponent (same for all)

    Hypothesis: ρ_sat,0,SPIRAL << ρ_sat,0,IRREGULAR
    """

    print("="*80)
    print("GALAXY-SPECIFIC ρ_sat,0 MODEL (SESSION #26)")
    print("="*80)
    print()

    # Extract data
    rho_central = data['rho_central']
    rho_sat = data['rho_sat']
    galaxy_types = type_data['galaxy_types']

    # Filter valid data
    valid = (rho_central > 0) & (rho_sat > 0) & np.isfinite(rho_central) & np.isfinite(rho_sat)
    rho_c = rho_central[valid]
    rho_s = rho_sat[valid]
    types = galaxy_types[valid]

    print(f"Valid galaxies: {np.sum(valid)}/{len(rho_central)}")
    print(f"  SPIRAL: {np.sum(types == 'SPIRAL')}")
    print(f"  IRREGULAR: {np.sum(types == 'IRREGULAR')}")
    print()

    # Define model with type-specific ρ_sat,0
    def rho_sat_model_typespecific(rho_c, types, rho_sat_0_spiral, rho_sat_0_irregular, rho_mag, delta):
        """Magnetic screening with type-specific ρ_sat,0."""
        rho_sat_0 = np.where(types == 'SPIRAL', rho_sat_0_spiral, rho_sat_0_irregular)
        return rho_sat_0 / (1.0 + (rho_c / rho_mag)**delta)

    # For curve_fit, need to flatten types into parameters
    def log_model_typespecific(rho_c, log_rho_sat_0_spiral, log_rho_sat_0_irregular, log_rho_mag, delta):
        rho_sat_0_spiral = 10**log_rho_sat_0_spiral
        rho_sat_0_irregular = 10**log_rho_sat_0_irregular
        rho_mag = 10**log_rho_mag

        result = rho_sat_model_typespecific(rho_c, types, rho_sat_0_spiral,
                                            rho_sat_0_irregular, rho_mag, delta)
        return np.log10(result)

    # Initial guess (from universal model, split ρ_sat,0)
    p0 = [5.0, 6.5, 2.5, 2.0]  # log(rho_sat_0_spiral), log(rho_sat_0_irregular), log(rho_mag), delta

    # Convert to log space
    log_rho_s = np.log10(rho_s)

    # Fit
    try:
        popt, pcov = curve_fit(log_model_typespecific, rho_c, log_rho_s, p0=p0,
                              bounds=([-np.inf, -np.inf, -np.inf, 0],
                                     [np.inf, np.inf, np.inf, 10]),
                              maxfev=10000)

        # Convert back
        rho_sat_0_spiral = 10**popt[0]
        rho_sat_0_irregular = 10**popt[1]
        rho_mag = 10**popt[2]
        delta = popt[3]

        # Uncertainties
        perr = np.sqrt(np.diag(pcov))
        rho_sat_0_spiral_err = rho_sat_0_spiral * perr[0] * np.log(10)
        rho_sat_0_irregular_err = rho_sat_0_irregular * perr[1] * np.log(10)
        rho_mag_err = rho_mag * perr[2] * np.log(10)
        delta_err = perr[3]

        # Predictions
        rho_sat_pred = rho_sat_model_typespecific(rho_c, types, rho_sat_0_spiral,
                                                   rho_sat_0_irregular, rho_mag, delta)

        # R²
        ss_res = np.sum((np.log10(rho_s) - np.log10(rho_sat_pred))**2)
        ss_tot = np.sum((np.log10(rho_s) - np.mean(np.log10(rho_s)))**2)
        r_squared = 1 - (ss_res / ss_tot)

        # RMS in log space
        rms_log = np.sqrt(np.mean((np.log10(rho_s) - np.log10(rho_sat_pred))**2))

        # Correlation
        correlation = np.corrcoef(np.log10(rho_s), np.log10(rho_sat_pred))[0, 1]

        print("✅ FIT SUCCESSFUL")
        print()
        print("Best-fit parameters:")
        print(f"  ρ_sat,0,SPIRAL     = {rho_sat_0_spiral:.2e} ± {rho_sat_0_spiral_err:.2e} M_☉/pc³")
        print(f"  ρ_sat,0,IRREGULAR  = {rho_sat_0_irregular:.2e} ± {rho_sat_0_irregular_err:.2e} M_☉/pc³")
        print(f"  ρ_mag              = {rho_mag:.2e} ± {rho_mag_err:.2e} M_☉/pc³")
        print(f"  δ                  = {delta:.2f} ± {delta_err:.2f}")
        print()
        print(f"  ρ_sat,0,IRREGULAR / ρ_sat,0,SPIRAL = {rho_sat_0_irregular / rho_sat_0_spiral:.1f}")
        print()
        print("Fit quality:")
        print(f"  R² = {r_squared:.3f}")
        print(f"  RMS(log) = {rms_log:.3f}")
        print(f"  Correlation: r = {correlation:.3f}")
        print()

        return {
            'rho_sat_0_spiral': rho_sat_0_spiral,
            'rho_sat_0_irregular': rho_sat_0_irregular,
            'rho_mag': rho_mag,
            'delta': delta,
            'rho_sat_0_spiral_err': rho_sat_0_spiral_err,
            'rho_sat_0_irregular_err': rho_sat_0_irregular_err,
            'rho_mag_err': rho_mag_err,
            'delta_err': delta_err,
            'r_squared': r_squared,
            'rms_log': rms_log,
            'correlation': correlation,
            'rho_sat_pred': rho_sat_pred,
            'valid_mask': valid,
            'types': types
        }

    except Exception as e:
        print(f"❌ FIT FAILED: {e}")
        return None


# ============================================================================
# PART 5: MODEL COMPARISON AND INTERPRETATION
# ============================================================================

def compare_models(universal_result, specific_result):
    """Compare universal vs galaxy-specific ρ_sat,0 models."""

    print("="*80)
    print("MODEL COMPARISON")
    print("="*80)
    print()

    if universal_result is None or specific_result is None:
        print("❌ Cannot compare - one or both fits failed")
        return

    # R² comparison
    r2_universal = universal_result['r_squared']
    r2_specific = specific_result['r_squared']
    delta_r2 = r2_specific - r2_universal

    print("Fit Quality Comparison:")
    print(f"  Universal model:        R² = {r2_universal:.3f}")
    print(f"  Galaxy-specific model:  R² = {r2_specific:.3f}")
    print(f"  Improvement:            ΔR² = {delta_r2:+.3f}")
    print()

    # RMS comparison
    rms_universal = universal_result['rms_log']
    rms_specific = specific_result['rms_log']
    delta_rms = rms_specific - rms_universal

    print(f"  Universal model:        RMS(log) = {rms_universal:.3f}")
    print(f"  Galaxy-specific model:  RMS(log) = {rms_specific:.3f}")
    print(f"  Improvement:            ΔRMS = {delta_rms:+.3f}")
    print()

    # Parameter comparison
    print("Parameter Comparison:")
    print(f"  Universal ρ_sat,0:          {universal_result['rho_sat_0']:.2e} M_☉/pc³")
    print(f"  Galaxy-specific ρ_sat,0,SPIRAL:    {specific_result['rho_sat_0_spiral']:.2e} M_☉/pc³")
    print(f"  Galaxy-specific ρ_sat,0,IRREGULAR: {specific_result['rho_sat_0_irregular']:.2e} M_☉/pc³")
    print()

    ratio = specific_result['rho_sat_0_irregular'] / specific_result['rho_sat_0_spiral']
    print(f"  ρ_sat,0,IRREGULAR / ρ_sat,0,SPIRAL = {ratio:.1f}×")
    print()

    # Interpretation
    print("="*80)
    print("INTERPRETATION")
    print("="*80)
    print()

    if delta_r2 > 0.05:
        print("✅ GALAXY-SPECIFIC MODEL SIGNIFICANTLY BETTER")
        print()
        print(f"  Improvement: ΔR² = {delta_r2:.3f} (> 0.05 threshold)")
        print()

        if ratio > 1:
            print("  ✅ HYPOTHESIS SUPPORTED")
            print()
            print("  Finding: ρ_sat,0,IRREGULAR > ρ_sat,0,SPIRAL")
            print(f"    Irregular galaxies have {ratio:.1f}× higher coherence threshold")
            print()
            print("  Physical Interpretation (Synchronism):")
            print("    - SPIRAL (NGC/UGC): Ordered rotation → Low turbulence")
            print("        → Low coherence threshold → Low ρ_sat,0")
            print("    - IRREGULAR (F/DDO): Chaotic motion → High turbulence")
            print("        → High coherence threshold → High ρ_sat,0")
            print()
            print("  Resolves Session #25 Paradox:")
            print("    - NGC has strong B-field (12 µG) BUT low ρ_sat,0")
            print("      → Overall ρ_sat = ρ_sat,0 / [screening] is LOW")
            print("    - F has weak B-field (4 µG) BUT high ρ_sat,0")
            print("      → Overall ρ_sat = ρ_sat,0 / [screening] is HIGH")
            print("    - ρ_sat,0 variation DOMINATES over B-field screening")
            print()
        else:
            print("  ⚠️ HYPOTHESIS REJECTED")
            print()
            print("  Finding: ρ_sat,0,SPIRAL > ρ_sat,0,IRREGULAR (opposite!)")
            print("    Does NOT explain paradox")
            print()

    elif delta_r2 > 0.01:
        print("⚠️ GALAXY-SPECIFIC MODEL SLIGHTLY BETTER")
        print()
        print(f"  Improvement: ΔR² = {delta_r2:.3f} (marginal, 0.01-0.05)")
        print("  May indicate weak galaxy-type dependence")
        print()

    else:
        print("❌ GALAXY-SPECIFIC MODEL NO BETTER")
        print()
        print(f"  Improvement: ΔR² = {delta_r2:.3f} (< 0.01, negligible)")
        print("  ρ_sat,0 appears UNIVERSAL (not galaxy-specific)")
        print("  Session #25 paradox remains unresolved")
        print()


# ============================================================================
# PART 6: BY-TYPE ANALYSIS
# ============================================================================

def analyze_by_type(data, specific_result):
    """Analyze predictions separately for SPIRAL vs IRREGULAR."""

    if specific_result is None:
        return

    print("="*80)
    print("BY-TYPE ANALYSIS")
    print("="*80)
    print()

    # Extract data
    rho_sat = data['rho_sat'][specific_result['valid_mask']]
    rho_sat_pred = specific_result['rho_sat_pred']
    types = specific_result['types']

    # Separate by type
    spiral_mask = types == 'SPIRAL'
    irregular_mask = types == 'IRREGULAR'

    # SPIRAL analysis
    print("SPIRAL GROUP (NGC + UGC):")
    print(f"  n = {np.sum(spiral_mask)}")

    rho_s_spiral = rho_sat[spiral_mask]
    rho_pred_spiral = rho_sat_pred[spiral_mask]

    # R² for spirals only
    ss_res_spiral = np.sum((np.log10(rho_s_spiral) - np.log10(rho_pred_spiral))**2)
    ss_tot_spiral = np.sum((np.log10(rho_s_spiral) - np.mean(np.log10(rho_s_spiral)))**2)
    r2_spiral = 1 - (ss_res_spiral / ss_tot_spiral)

    print(f"  R² = {r2_spiral:.3f}")
    print(f"  Median observed:  {np.median(rho_s_spiral):.2e} M_☉/pc³")
    print(f"  Median predicted: {np.median(rho_pred_spiral):.2e} M_☉/pc³")
    print()

    # IRREGULAR analysis
    print("IRREGULAR GROUP (F + DDO + other):")
    print(f"  n = {np.sum(irregular_mask)}")

    rho_s_irregular = rho_sat[irregular_mask]
    rho_pred_irregular = rho_sat_pred[irregular_mask]

    # R² for irregulars only
    ss_res_irregular = np.sum((np.log10(rho_s_irregular) - np.log10(rho_pred_irregular))**2)
    ss_tot_irregular = np.sum((np.log10(rho_s_irregular) - np.mean(np.log10(rho_s_irregular)))**2)
    r2_irregular = 1 - (ss_res_irregular / ss_tot_irregular)

    print(f"  R² = {r2_irregular:.3f}")
    print(f"  Median observed:  {np.median(rho_s_irregular):.2e} M_☉/pc³")
    print(f"  Median predicted: {np.median(rho_pred_irregular):.2e} M_☉/pc³")
    print()


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print()
    print("╔" + "="*78 + "╗")
    print("║" + " "*78 + "║")
    print("║" + "  Synchronism Research - Session #26".center(78) + "║")
    print("║" + "  Galaxy-Specific ρ_sat,0 Model Test".center(78) + "║")
    print("║" + " "*78 + "║")
    print("╚" + "="*78 + "╝")
    print()

    # Load data
    data = load_session20_data()

    # Classify galaxy types
    type_data = classify_galaxy_types(data)

    # Fit universal model (Session #22 baseline)
    universal_result = fit_universal_model(data)

    # Fit galaxy-specific ρ_sat,0 model (Session #26 test)
    specific_result = fit_galaxy_specific_rho_sat0(data, type_data)

    # Compare models
    compare_models(universal_result, specific_result)

    # By-type analysis
    analyze_by_type(data, specific_result)

    print()
    print("SESSION #26 COMPLETE")
    print("="*80)
    print()

    if specific_result is not None and universal_result is not None:
        delta_r2 = specific_result['r_squared'] - universal_result['r_squared']

        if delta_r2 > 0.05:
            ratio = specific_result['rho_sat_0_irregular'] / specific_result['rho_sat_0_spiral']
            if ratio > 1:
                print("✅ HYPOTHESIS SUPPORTED")
                print(f"  Galaxy-specific ρ_sat,0 improves fit: ΔR² = {delta_r2:+.3f}")
                print(f"  ρ_sat,0,IRREGULAR / ρ_sat,0,SPIRAL = {ratio:.1f}×")
                print("  → Explains Session #25 B-field paradox")
            else:
                print("⚠️ MODEL IMPROVED BUT PARADOX NOT RESOLVED")
                print(f"  ΔR² = {delta_r2:+.3f} but ρ_sat,0,SPIRAL > ρ_sat,0,IRREGULAR")
        else:
            print("❌ HYPOTHESIS NOT SUPPORTED")
            print(f"  Galaxy-specific ρ_sat,0 shows minimal improvement: ΔR² = {delta_r2:+.3f}")
            print("  → ρ_sat,0 appears universal")

    print()
    print("Ready for Session #26 documentation and commit.")
    print()
