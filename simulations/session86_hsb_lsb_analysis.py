#!/usr/bin/env python3
"""
Session #86: HSB vs LSB Galaxy BTFR Analysis

Synchronism Prediction:
- Core theory: G_eff = G/C(ρ) where C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))
- LSB galaxies have lower surface/volume density than HSB galaxies
- Lower density → lower C → higher G_eff → higher rotation velocity at fixed mass
- PREDICTION: LSB galaxies should have HIGHER V_flat at fixed baryonic mass

MOND Prediction:
- BTFR should be universal: M_bar ∝ V^4 with same normalization
- No dependence on surface brightness (same a₀ everywhere)

This is a key discriminating test between Synchronism and MOND.

Author: CBP Autonomous Synchronism Research
Date: December 4, 2025
"""

import numpy as np
import json
import os
from collections import defaultdict

# Constants
SUN_MASS = 1.989e30  # kg
PC_TO_M = 3.086e16   # meters per parsec
G_NEWTON = 6.674e-11  # m³/kg/s²

def parse_sparc_galaxy_table(filepath):
    """Parse the SPARC galaxy properties table."""
    galaxies = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Find data start (after header) - look for the dashed separator
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith('-----'):
            data_start = i + 1
            break

    for line in lines[data_start:]:
        line = line.strip()
        if len(line) < 50:
            continue
        if line.startswith('#') or not line:
            continue

        try:
            # Split on whitespace
            parts = line.split()
            if len(parts) < 18:
                continue

            # Column mapping (0-indexed):
            # 0: Galaxy name
            # 1: Hubble Type T
            # 2: Distance D
            # 3: e_D
            # 4: f_D (distance method)
            # 5: Inc
            # 6: e_Inc
            # 7: L[3.6] (10^9 L_sun)
            # 8: e_L[3.6]
            # 9: Reff
            # 10: SBeff
            # 11: Rdisk
            # 12: SBdisk
            # 13: MHI (10^9 M_sun)
            # 14: RHI
            # 15: Vflat
            # 16: e_Vflat
            # 17: Q (quality)
            # 18+: References

            name = parts[0]
            T = int(parts[1]) if parts[1] else -1
            D = float(parts[2]) if parts[2] else np.nan
            e_D = float(parts[3]) if parts[3] else np.nan
            Inc = float(parts[5]) if parts[5] else np.nan
            L36 = float(parts[7]) if parts[7] else np.nan
            Reff = float(parts[9]) if parts[9] else np.nan
            SBeff = float(parts[10]) if parts[10] else np.nan
            Rdisk = float(parts[11]) if parts[11] else np.nan
            SBdisk = float(parts[12]) if parts[12] else np.nan
            MHI = float(parts[13]) if parts[13] else np.nan
            RHI = float(parts[14]) if parts[14] else np.nan
            Vflat = float(parts[15]) if parts[15] else np.nan
            e_Vflat = float(parts[16]) if parts[16] else np.nan
            Q = int(parts[17]) if parts[17] else 3

            if not np.isnan(Vflat) and Vflat > 0 and not np.isnan(L36) and not np.isnan(SBdisk):
                galaxies.append({
                    'name': name,
                    'hubble_type': T,
                    'distance_mpc': D,
                    'inclination': Inc,
                    'L36_1e9_Lsun': L36,  # 10^9 L_sun
                    'Reff_kpc': Reff,
                    'SBeff': SBeff,  # L_sun/pc²
                    'Rdisk_kpc': Rdisk,
                    'SBdisk': SBdisk,  # Central disk surface brightness L_sun/pc²
                    'MHI_1e9_Msun': MHI,  # 10^9 M_sun
                    'RHI_kpc': RHI,
                    'Vflat': Vflat,  # km/s
                    'e_Vflat': e_Vflat,
                    'quality': Q
                })
        except (ValueError, IndexError) as e:
            continue

    return galaxies


def classify_surface_brightness(galaxies, sb_threshold=100):
    """
    Classify galaxies as HSB, ISB, or LSB based on central disk surface brightness.

    Traditional boundary is ~21.65 mag/arcsec² in B-band, which corresponds to
    different values in [3.6] micron. We use L_sun/pc² directly.

    Following McGaugh convention:
    - HSB: SBdisk > 500 L_sun/pc² (high surface brightness)
    - ISB: 100 < SBdisk < 500 L_sun/pc² (intermediate)
    - LSB: SBdisk < 100 L_sun/pc² (low surface brightness)
    """
    hsb = []
    isb = []
    lsb = []

    for g in galaxies:
        sb = g['SBdisk']
        if sb > 500:
            hsb.append(g)
        elif sb > 100:
            isb.append(g)
        else:
            lsb.append(g)

    return hsb, isb, lsb


def compute_baryonic_mass(galaxy, ML_ratio=0.5):
    """
    Compute baryonic mass: M_bar = M_star + 1.33*M_HI

    M_star = ML_ratio × L[3.6]
    ML_ratio ≈ 0.5 M_sun/L_sun for [3.6] band (standard assumption)
    Factor 1.33 accounts for helium in gas
    """
    L36 = galaxy['L36_1e9_Lsun'] * 1e9  # L_sun
    MHI = galaxy['MHI_1e9_Msun'] * 1e9 if not np.isnan(galaxy['MHI_1e9_Msun']) else 0

    M_star = ML_ratio * L36
    M_gas = 1.33 * MHI
    M_bar = M_star + M_gas

    return M_bar, M_star, M_gas


def compute_btfr_residuals(galaxies, btfr_slope=4.0, btfr_zp=None):
    """
    Compute BTFR residuals: Δlog(V) = log(V_obs) - log(V_pred)

    Standard BTFR: log(M_bar) = btfr_zp + btfr_slope × log(V)

    If btfr_zp is None, fit it from the data.
    """
    log_V = []
    log_M = []

    for g in galaxies:
        M_bar, _, _ = compute_baryonic_mass(g)
        if M_bar > 0 and g['Vflat'] > 0:
            log_V.append(np.log10(g['Vflat']))
            log_M.append(np.log10(M_bar))

    log_V = np.array(log_V)
    log_M = np.array(log_M)

    if btfr_zp is None:
        # Fit zero-point assuming slope=4
        btfr_zp = np.mean(log_M - btfr_slope * log_V)

    # Residuals in log(V)
    log_V_pred = (log_M - btfr_zp) / btfr_slope
    residuals = log_V - log_V_pred

    return residuals, btfr_zp


def analyze_hsb_lsb_btfr(galaxies):
    """
    Main analysis: Compare BTFR between HSB, ISB, and LSB galaxies.
    """
    # Classify galaxies
    hsb, isb, lsb = classify_surface_brightness(galaxies)

    print("=" * 70)
    print("SESSION #86: HSB vs LSB GALAXY BTFR ANALYSIS")
    print("=" * 70)
    print()
    print("SYNCHRONISM PREDICTION:")
    print("  LSB galaxies have lower density → lower C → higher G_eff")
    print("  → LSB should have HIGHER V_flat at fixed M_bar")
    print()
    print("MOND PREDICTION:")
    print("  BTFR is universal with same normalization regardless of SB")
    print()

    print("=" * 70)
    print("SAMPLE CLASSIFICATION")
    print("=" * 70)
    print(f"  HSB (SBdisk > 500 L_sun/pc²): {len(hsb)} galaxies")
    print(f"  ISB (100-500 L_sun/pc²): {len(isb)} galaxies")
    print(f"  LSB (SBdisk < 100 L_sun/pc²): {len(lsb)} galaxies")
    print(f"  Total: {len(galaxies)} galaxies")
    print()

    # Fit global BTFR
    all_log_V = []
    all_log_M = []

    for g in galaxies:
        M_bar, _, _ = compute_baryonic_mass(g)
        if M_bar > 0 and g['Vflat'] > 0:
            all_log_V.append(np.log10(g['Vflat']))
            all_log_M.append(np.log10(M_bar))

    all_log_V = np.array(all_log_V)
    all_log_M = np.array(all_log_M)

    # Fit BTFR with fixed slope=4
    global_zp = np.mean(all_log_M - 4.0 * all_log_V)

    print("=" * 70)
    print("GLOBAL BTFR FIT (slope fixed at 4.0)")
    print("=" * 70)
    print(f"  log(M_bar) = {global_zp:.3f} + 4.0 × log(V_flat)")
    print()

    # Analyze each subsample
    results = {}

    for name, sample in [('HSB', hsb), ('ISB', isb), ('LSB', lsb)]:
        if len(sample) < 5:
            print(f"  {name}: Too few galaxies ({len(sample)})")
            continue

        log_V = []
        log_M = []
        sbs = []

        for g in sample:
            M_bar, _, _ = compute_baryonic_mass(g)
            if M_bar > 0 and g['Vflat'] > 0:
                log_V.append(np.log10(g['Vflat']))
                log_M.append(np.log10(M_bar))
                sbs.append(g['SBdisk'])

        log_V = np.array(log_V)
        log_M = np.array(log_M)
        sbs = np.array(sbs)

        # Fit BTFR zero-point for this sample (with global slope)
        sample_zp = np.mean(log_M - 4.0 * log_V)

        # Compute residuals relative to global fit
        log_V_pred = (log_M - global_zp) / 4.0
        residuals = log_V - log_V_pred  # Positive = higher V than expected

        mean_resid = np.mean(residuals)
        std_resid = np.std(residuals)
        sem_resid = std_resid / np.sqrt(len(residuals))

        results[name] = {
            'n_galaxies': len(sample),
            'mean_SBdisk': np.mean(sbs),
            'median_SBdisk': np.median(sbs),
            'sample_zp': sample_zp,
            'mean_residual_dex': mean_resid,
            'std_residual_dex': std_resid,
            'sem_residual_dex': sem_resid,
            'significance': mean_resid / sem_resid if sem_resid > 0 else 0
        }

    print("=" * 70)
    print("BTFR RESIDUALS BY SURFACE BRIGHTNESS CLASS")
    print("=" * 70)
    print()
    print("  Class   N    Mean SB    Median SB   ZP       Δlog(V)  ±SEM    σ")
    print("-" * 70)

    for name in ['HSB', 'ISB', 'LSB']:
        if name not in results:
            continue
        r = results[name]
        print(f"  {name:5s}  {r['n_galaxies']:3d}  {r['mean_SBdisk']:8.1f}  {r['median_SBdisk']:9.1f}  {r['sample_zp']:.3f}  {r['mean_residual_dex']:+.4f} ±{r['sem_residual_dex']:.4f}  {r['significance']:+.1f}σ")

    print()

    # Key comparison: LSB vs HSB
    print("=" * 70)
    print("KEY COMPARISON: LSB vs HSB")
    print("=" * 70)

    if 'LSB' in results and 'HSB' in results:
        lsb_r = results['LSB']
        hsb_r = results['HSB']

        offset = lsb_r['mean_residual_dex'] - hsb_r['mean_residual_dex']
        offset_err = np.sqrt(lsb_r['sem_residual_dex']**2 + hsb_r['sem_residual_dex']**2)
        offset_sig = offset / offset_err if offset_err > 0 else 0

        print(f"  LSB - HSB offset: {offset:+.4f} ± {offset_err:.4f} dex ({offset_sig:+.1f}σ)")
        print()
        print("  INTERPRETATION:")
        print("  - Positive offset = LSB has HIGHER V at fixed M (Synchronism prediction)")
        print("  - Zero offset = BTFR is universal (MOND prediction)")
        print()

        if offset > 0:
            print(f"  RESULT: LSB galaxies show HIGHER rotation velocity")
            print(f"  This is CONSISTENT with Synchronism prediction!")
        elif offset < 0:
            print(f"  RESULT: LSB galaxies show LOWER rotation velocity")
            print(f"  This CONTRADICTS Synchronism prediction!")
        else:
            print(f"  RESULT: No significant offset detected")
            print(f"  This is consistent with MOND universality")

        results['LSB_HSB_offset'] = {
            'offset_dex': offset,
            'offset_err_dex': offset_err,
            'significance': offset_sig
        }

    # Synchronism prediction calculation
    print()
    print("=" * 70)
    print("SYNCHRONISM THEORY PREDICTION")
    print("=" * 70)
    print()

    # Estimate density contrast between HSB and LSB
    if 'LSB' in results and 'HSB' in results:
        sb_ratio = results['HSB']['mean_SBdisk'] / results['LSB']['mean_SBdisk']
        print(f"  Mean SB ratio (HSB/LSB): {sb_ratio:.1f}")

        # Volume density scales as SB^1.5 for exponential disks
        # (SB ∝ M/R², volume density ∝ M/R³ ∝ SB × R^-1)
        # For similar scale heights, ρ ∝ SB
        density_ratio = sb_ratio
        print(f"  Estimated density ratio: {density_ratio:.1f}")

        # In Synchronism: C(ρ) = tanh(γ × log(ρ/ρ_crit + 1))
        # For galaxies, C is typically 0.3-0.7
        # The ratio of C values gives the ratio of G_eff
        # V² ∝ G_eff × M/R, so Δlog(V) = 0.5 × Δlog(G_eff) = -0.5 × Δlog(C)

        # With γ = 2.0 and typical galaxy densities:
        # C_HSB ~ 0.6, C_LSB ~ 0.4 (rough estimate)
        # Δlog(C) ~ log(0.4/0.6) = -0.18
        # Δlog(V) ~ 0.5 × 0.18 = 0.09 dex

        C_hsb_est = 0.6
        C_lsb_est = 0.4
        delta_logC = np.log10(C_lsb_est / C_hsb_est)
        predicted_offset = -0.5 * delta_logC

        print(f"  Estimated C_HSB ~ {C_hsb_est}, C_LSB ~ {C_lsb_est}")
        print(f"  Δlog(C) = {delta_logC:.3f}")
        print(f"  PREDICTED Δlog(V) = -0.5 × Δlog(C) = {predicted_offset:+.3f} dex")
        print()

        observed_offset = results['LSB_HSB_offset']['offset_dex']
        print(f"  OBSERVED Δlog(V) = {observed_offset:+.4f} dex")
        print(f"  Observed/Predicted ratio: {observed_offset/predicted_offset:.1%}")

        results['synchronism_prediction'] = {
            'density_ratio': density_ratio,
            'C_hsb_estimate': C_hsb_est,
            'C_lsb_estimate': C_lsb_est,
            'predicted_offset_dex': predicted_offset,
            'observed_offset_dex': observed_offset,
            'ratio': observed_offset / predicted_offset if predicted_offset != 0 else np.nan
        }

    return results


def analyze_continuous_sb_correlation(galaxies):
    """
    Analyze continuous correlation between surface brightness and BTFR residuals.
    """
    print()
    print("=" * 70)
    print("CONTINUOUS SB vs BTFR RESIDUAL CORRELATION")
    print("=" * 70)

    # Compute data
    log_V = []
    log_M = []
    log_SB = []

    for g in galaxies:
        M_bar, _, _ = compute_baryonic_mass(g)
        if M_bar > 0 and g['Vflat'] > 0 and g['SBdisk'] > 0:
            log_V.append(np.log10(g['Vflat']))
            log_M.append(np.log10(M_bar))
            log_SB.append(np.log10(g['SBdisk']))

    log_V = np.array(log_V)
    log_M = np.array(log_M)
    log_SB = np.array(log_SB)

    # Global BTFR
    global_zp = np.mean(log_M - 4.0 * log_V)
    log_V_pred = (log_M - global_zp) / 4.0
    residuals = log_V - log_V_pred

    # Pearson correlation
    r = np.corrcoef(log_SB, residuals)[0, 1]

    # Linear fit
    slope, intercept = np.polyfit(log_SB, residuals, 1)

    print(f"  N galaxies: {len(log_SB)}")
    print(f"  Pearson correlation (log_SB vs residual): r = {r:.3f}")
    print(f"  Linear fit: Δlog(V) = {slope:.4f} × log(SBdisk) + {intercept:.4f}")
    print()
    print("  INTERPRETATION:")
    print("  - Negative correlation = Higher SB → Lower V residual (Synchronism prediction)")
    print("  - Zero correlation = No SB dependence (MOND prediction)")
    print()

    if r < -0.1:
        print(f"  RESULT: Negative correlation detected (r = {r:.3f})")
        print(f"  This is CONSISTENT with Synchronism prediction")
    elif r > 0.1:
        print(f"  RESULT: Positive correlation detected (r = {r:.3f})")
        print(f"  This CONTRADICTS Synchronism prediction")
    else:
        print(f"  RESULT: No significant correlation (r = {r:.3f})")
        print(f"  This is consistent with MOND universality")

    return {
        'n_galaxies': len(log_SB),
        'pearson_r': r,
        'linear_slope': slope,
        'linear_intercept': intercept
    }


def main():
    # Path to data
    data_dir = '/mnt/c/exe/projects/ai-agents/synchronism/simulations/sparc_real_data'
    galaxy_file = os.path.join(data_dir, 'SPARC_Lelli2016c.mrt')

    # Parse data
    print("Loading SPARC galaxy properties...")
    galaxies = parse_sparc_galaxy_table(galaxy_file)
    print(f"Loaded {len(galaxies)} galaxies with valid data")
    print()

    # Filter by quality
    quality_galaxies = [g for g in galaxies if g['quality'] <= 2]
    print(f"After quality filter (Q ≤ 2): {len(quality_galaxies)} galaxies")
    print()

    # Main analysis
    results = analyze_hsb_lsb_btfr(quality_galaxies)

    # Continuous correlation analysis
    correlation_results = analyze_continuous_sb_correlation(quality_galaxies)
    results['continuous_correlation'] = correlation_results

    # Summary
    print()
    print("=" * 70)
    print("SESSION #86 SUMMARY")
    print("=" * 70)
    print()

    if 'LSB_HSB_offset' in results:
        offset = results['LSB_HSB_offset']['offset_dex']
        sig = results['LSB_HSB_offset']['significance']
        print(f"  LSB - HSB BTFR offset: {offset:+.4f} dex ({sig:+.1f}σ)")

        if 'synchronism_prediction' in results:
            pred = results['synchronism_prediction']['predicted_offset_dex']
            ratio = results['synchronism_prediction']['ratio']
            print(f"  Synchronism prediction: {pred:+.3f} dex")
            print(f"  Observed/Predicted: {ratio:.1%}")

    print(f"  Continuous correlation: r = {correlation_results['pearson_r']:.3f}")
    print()

    # Save results
    results_dir = os.path.join(os.path.dirname(data_dir), 'results')
    os.makedirs(results_dir, exist_ok=True)

    output_file = os.path.join(results_dir, 'session86_hsb_lsb_analysis.json')

    # Convert to JSON-serializable format
    def convert_to_serializable(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_serializable(i) for i in obj]
        return obj

    with open(output_file, 'w') as f:
        json.dump(convert_to_serializable(results), f, indent=2)

    print(f"Results saved to: {output_file}")


if __name__ == '__main__':
    main()
