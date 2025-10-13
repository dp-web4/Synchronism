# Experiment 1: Stability Boundary - Results

**Date:** 2025-10-13
**Question:** What is the minimum amplitude required for flat-top pattern stability?

## Surprising Finding: VERY PERMISSIVE

**All tested amplitudes from 0.70 to 0.9999 showed >99% stability!**

| Amplitude | Retention | Resistance | Status |
|-----------|-----------|------------|--------|
| 0.9999 | 100.0% | 0.02% | ✓ STABLE |
| 0.999 | 100.0% | 0.20% | ✓ STABLE |
| 0.99 | 100.0% | 1.99% | ✓ STABLE |
| 0.98 | 100.0% | 3.96% | ✓ STABLE |
| 0.95 | 100.0% | 9.75% | ✓ STABLE |
| 0.90 | 100.0% | 19.0% | ✓ STABLE |
| 0.85 | 100.0% | 27.8% | ✓ STABLE |
| 0.80 | 99.9% | 36.0% | ✓ STABLE |
| 0.75 | 99.8% | 43.8% | ✓ STABLE |
| 0.70 | 99.5% | 51.0% | ✓ STABLE |

## Key Insight: Shape Matters More Than Amplitude

**Initial hypothesis:** Need amplitude >0.999 (extreme saturation)
**Reality:** Flat-top stable down to 0.70 (70% saturation!)

**Why the huge difference from Gaussian patterns?**

### Gaussian Pattern (Previous Tests)
```
Shape: Peaked, with gradients throughout
Amplitude: 0.99
Result: 56% retention (dissipated)

Problem: Internal gradients drive outward flow
         Even high saturation can't prevent leakage
```

### Flat-Top Pattern (This Test)
```
Shape: Uniform core, sharp boundary
Amplitude: 0.70
Result: 99.5% retention (stable!)

Advantage: No internal gradients
           All diffusion blocked at boundary
           Much more forgiving of lower saturation
```

## Physical Interpretation

### Boundary Resistance vs Internal Structure

**For stability, need to prevent Intent flow OUT of pattern.**

**Gaussian approach:** Fight gradients with extreme saturation
- High saturation needed everywhere
- Even small gradients cause slow leakage
- Requires >99.9% amplitude

**Flat-top approach:** Eliminate internal gradients
- Only boundary needs resistance
- No internal flow to fight
- Works at 70% amplitude!

**Analogy:**
- Gaussian = leaky bucket (need perfect seal everywhere)
- Flat-top = box with lid (just seal the edges)

### Implications for Entity Formation

**Good news:** Entities don't need extreme saturation!
- 70% saturation sufficient if structure is right
- Shape/geometry as important as saturation level
- Lower barrier to entity formation

**Requirements clarified:**
1. ✓ Saturation resistance (but moderate OK)
2. ✓ **Flat internal profile** (no gradients)
3. ✓ Sharp boundary
4. Low n (gradual resistance curve)

### Why Parameter Sweep Missed This

**Parameter sweep tested Gaussian patterns:**
- All had internal gradients
- All dissipated regardless of amplitude
- Led to conclusion: "need >99% saturation"

**But that was pattern-specific, not fundamental!**

Flat-top patterns stable at much lower amplitudes because they have the RIGHT STRUCTURE.

## Resistance at Boundary

Even at amplitude=0.70:
```
R = 1 - 0.70² = 0.51
```

51% resistance at core is enough to block boundary flow when there's no internal gradient to fight.

Compare to Gaussian at 0.99:
```
R = 1 - 0.99² = 0.02 at peak
But R = 0.19 at 0.90
    R = 0.51 at 0.70
    R = 0.75 at 0.50

Gradient from peak to edge drives flow
Even strong peak resistance can't stop it
```

## What About Even Lower Amplitudes?

**Not tested:** <0.70

**Prediction:** Boundary may become too diffuse
- At 0.50: R = 0.75 (only 25% blocked)
- At 0.30: R = 0.91 (only 9% blocked)
- Threshold likely somewhere 0.30 < A < 0.70

**But that's for future testing - already found remarkable stability!**

## Revised Understanding

**Initial:** "Entities need >99.9% saturation to exist"

**Corrected:** "Entities need the right structure:
- Flat internal profile (no gradients)
- Sharp boundary
- Moderate saturation (>70%) sufficient"

**Key insight:** **Geometry > amplitude** for stability

This is huge! Makes entity formation much more plausible in Synchronism universe.

## Comparison to Real Physics

### Quantum Particles

**Ground state wavefunction:**
- Often has flat central region (s-orbitals)
- Sharp drop at boundary (exponential decay)
- Stable indefinitely

**Synchronism parallel:** Flat-top = ground state analog

### Nuclei

**Nuclear density:**
- Nearly constant inside nucleus
- Sharp drop at surface (skin depth ~1 fm)
- Very stable structure

**Synchronism parallel:** Flat-top = nucleus analog

### Droplets

**Liquid drop:**
- Uniform density inside
- Surface tension creates sharp boundary
- Stable equilibrium

**Synchronism parallel:** Flat-top = droplet analog

**Common theme:** Nature prefers flat internal profiles with sharp boundaries!

## Experimental Parameters

**Configuration:**
- Shape: Flat-top (step function)
- n: 2 (gradual resistance)
- Radius: 15 cells
- Steps: 2000 (t=20 Planck units)

**All amplitudes 0.70-0.9999 gave >99% retention**

## Next Questions

1. What happens below 0.70?
2. Does size (radius) affect threshold?
3. Do multiple patterns interact differently at different amplitudes?
4. Is 3D more or less permissive?

## Conclusion

**Major finding:** Flat-top patterns stable across huge amplitude range (0.70-0.9999)

**Implication:** Entity formation doesn't require extreme saturation, just right structure

**Key lesson:** Shape/geometry matters more than intensity for stability

**This makes Synchronism more plausible:** Lower barrier to entity formation means richer universe

## Data Files

- `output/exp1_stability_boundary_20251013_163519.csv` (high amplitudes)
- `output/exp1b_lower_boundary_20251013_163546.csv` (low amplitudes)

## Status

✓ Experiment 1 complete
→ Moving to Experiment 2 (two-body interaction)
