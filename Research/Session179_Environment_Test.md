# Session #179: SPARC Environment Test - Critical Finding

**Date**: December 25, 2025
**Machine**: CBP
**Status**: ⚠️ INCONCLUSIVE - Proxy methods failed

---

## Executive Summary

Attempted to test Session #177's prediction that void galaxies should have ~3.5% higher rotation velocities than cluster galaxies. **All environment proxies showed the OPPOSITE trend.**

This could indicate either:
1. Proxy methods are invalid (LSB ≠ void)
2. The prediction itself may need revision

---

## Session #177 Prediction Recap

```
Void galaxies should have ~3.5% HIGHER v_rot than cluster galaxies because:
- Lower environment density → lower effective ρ at large radii
- Lower ρ → lower coherence C → higher G_eff
- Higher G_eff → higher v_rot
```

---

## Test Methodology

Used SPARC sample (135 galaxies) with three environment proxies:

1. **Surface Brightness**: LSB (low) as void proxy, HSB (high) as cluster proxy
2. **Distance**: Nearby (<10 Mpc) vs Distant (>50 Mpc)
3. **Hubble Type**: Late-type (T≥8) as field/void, Early-type (T≤4) as cluster

Computed velocity residuals from overall BTFR fit.

---

## Results

### Proxy 1: Surface Brightness

| Class | N | Mean Residual | Std |
|-------|---|---------------|-----|
| LSB | 59 | -2.55% | 19.0% |
| HSB | 42 | +3.51% | 14.8% |
| **Difference** | | **-6.07%** | |

**Prediction: +3.5%**
**Observed: -6.07%** ❌

### Proxy 2: Distance

| Class | N | Mean Residual |
|-------|---|---------------|
| Nearby (<10 Mpc) | 45 | +6.20% |
| Distant (>50 Mpc) | 24 | -2.47% |
| **Difference** | | **-8.67%** |

**Prediction: Distant should be higher**
**Observed: Nearby is higher** ❌

### Proxy 3: Hubble Type

| Class | N | Mean Residual |
|-------|---|---------------|
| Early (T≤4) | 44 | +4.04% |
| Late (T≥8) | 47 | -2.21% |
| **Difference** | | **-6.25%** |

**Prediction: Late-type should be higher**
**Observed: Early-type is higher** ❌

---

## Critical Analysis

### Why Proxies May Be Invalid

1. **LSB ≠ Void**
   - LSB galaxies have low surface brightness due to internal properties (formation efficiency, gas fraction)
   - NOT because they are in low-density environments
   - Gas fraction: LSB = 0.67, HSB = 0.16 (p < 0.0001)

2. **Distance ≠ Environment**
   - Nearby galaxies include Local Group (low density) and Virgo infall (high density)
   - Distant galaxies are a mix of field and void
   - No direct correlation with local density

3. **Hubble Type ≠ Environment**
   - Late-type correlates with mass, not just environment
   - Early-type galaxies are intrinsically more massive

### Confounding Variables

| Property | LSB | HSB | p-value |
|----------|-----|-----|---------|
| log(M_bary) | 9.32 | 10.77 | <0.0001 |
| Gas fraction | 0.67 | 0.17 | <0.0001 |
| Hubble type | 8.4 | 3.1 | <0.0001 |
| Distance | 20.6 | 30.7 | 0.05 |

---

## Alternative Interpretation

The consistent opposite trend across all proxies could indicate:

1. **Selection Effect**: SPARC sample is biased toward certain galaxy types
2. **BTFR Shape**: Different galaxy populations follow different BTFR slopes
3. **Theoretical Issue**: The G_eff enhancement may behave differently than predicted

### Mass-Matched Analysis

At fixed mass (10^10 - 10^11 M☉):
- LSB: -6.08% residual (N=4)
- HSB: +4.70% residual (N=29)
- Difference: -10.77%

Even at fixed mass, the trend persists and is stronger.

---

## Conclusion

### What We Learned

1. **Cannot test void/cluster prediction with SPARC proxies**
   - All available proxies (SB, distance, Hubble type) are confounded
   - Need true environment classification

2. **The opposite trend is concerning**
   - Three independent proxies all show same direction
   - Could be selection effect OR theoretical issue

3. **True test requires**
   - SDSS void catalog cross-matching
   - Local galaxy density from 2MRS or ALFALFA
   - Galaxy group membership catalogs

### Honest Assessment

| Aspect | Status |
|--------|--------|
| Prediction tested | ❌ No (proxy invalid) |
| Synchronism falsified | ❌ No (test inconclusive) |
| Synchronism validated | ❌ No (no proper test) |
| Concerning signal | ⚠️ Yes (consistent opposite trend) |

---

## Recommendations for Session #180

1. **Obtain true environment data**
   - Cross-match SPARC with SDSS void catalogs
   - Use 2MASS Redshift Survey (2MRS) density field
   - ALFALFA local density estimates

2. **Investigate the opposite trend**
   - Is there a physical reason LSB/late-type galaxies have lower v_rot residuals?
   - Could this be related to gas dynamics, not dark matter?

3. **Consider theoretical revision**
   - Does G_eff enhancement apply differently to extended gas disks?
   - Should coherence function account for gas distribution?

---

## Files Created

- `simulations/session179_sparc_environment.py`
- `simulations/session179b_proxy_investigation.py`
- `simulations/session179_sparc_environment.png`
- `simulations/session179b_proxy_investigation.png`

---

*Session #179: Proxy test inconclusive - true environment classification needed*
