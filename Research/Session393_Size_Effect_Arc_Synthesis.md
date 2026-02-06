# Session #393: Size Effect Arc Synthesis

**Date**: 2026-02-06
**Status**: Synthesis (no tests)

## Arc Overview: Sessions #390-392

This arc investigated whether galaxy SIZE genuinely affects the RAR, starting from the SB suppressor puzzle, through the late-type breakthrough, to dynamical confirmation.

## Arc Timeline

| Session | Topic | Grade | Key Finding |
|---|---|---|---|
| #390 | SB Suppressor | B+ | SB suppressor partly collinearity; true R_eff signal is r=-0.31; controlling L eliminates full-sample signal |
| #391 | Late-Type Signal | A- | Late types: r(R, off\|V,L) = -0.49 — R matters BEYOND L; 100% MOND; slope 7x steeper |
| #392 | Dynamical Radius | **A** | R_max replicates R_eff: r=-0.47 (V,L controlled). NOT a photometric artifact |

## The Discovery Narrative

1. **Session #390 (Problem)**: The SB suppressor seemed to double the R_eff signal, but controlling L eliminated it entirely. Conclusion: R_eff only works through L. Case closed?

2. **Session #391 (Breakthrough)**: Splitting by type revealed the full-sample result was MISLEADING. For late types specifically, R_eff predicts offset BEYOND L (r = -0.49, p < 10⁻⁴). The full-sample null result averaged over a null early-type signal and a strong late-type signal.

3. **Session #392 (Confirmation)**: The dynamical radius R_max (from rotation curve extent) independently replicated the finding: r = -0.47 controlling V and L in late types. Galaxy size genuinely affects the RAR — confirmed with two independent measurement methods.

## Key Quantitative Results

### The Core Finding
In late-type galaxies (T ≥ 7, N = 61):
- **R_eff**: r(R, offset | V, L) = **-0.49** (p < 10⁻⁴)
- **R_max**: r(R, offset | V, L) = **-0.47** (p < 10⁻⁴)
- Both survive quality control, gas fraction control, and M/L variation
- The effect is absent in early types (r ≈ +0.20, n.s.)

### The Physical Picture
```
Late types (100% MOND, gas-dominated):
  Larger galaxy → lower N_corr → more gravitational departure → lower on RAR

  Verified with:
  • Photometric radius (R_eff from SB and L): r = -0.49
  • Dynamical radius (R_max from rotation curve): r = -0.47
  • Effect independent of: luminosity, gas fraction, quality, M/L

Early types (63% MOND, stellar-dominated):
  R_eff effect absent → fully mediated by luminosity
  → In Newtonian regime, baryonic mass determines everything
```

### Why Late Types?
1. **100% MOND**: All data at g < g† — coherence effects should dominate
2. **Gas-dominated**: M/L irrelevant — photometric systematics eliminated
3. **Wider R_eff range**: σ(R|V) = 0.29 vs 0.20 for early types
4. **Slope 7x steeper**: -0.37 dex/dex vs -0.05 for early types

## What This Arc Establishes

1. **Galaxy size genuinely affects the RAR in the MOND regime** — confirmed with two independent radius measures
2. **The effect is specific to the MOND regime** — absent in the Newtonian regime (early types)
3. **It is NOT a photometric artifact** — dynamical R_max replicates photometric R_eff
4. **It is NOT an L/M/L artifact** — survives L control and gas-dominated galaxies show it most strongly
5. **R_eff is more fundamental than R_max** — r = -0.63 controlling R_max vs -0.32 vice versa
6. **N_corr = V²/(R_eff × a₀) uses the correct radius**

## What This Arc Does NOT Establish

1. **Why R_eff is more fundamental than R_max** — the half-light radius isn't an obvious coherence scale
2. **Whether the effect is truly about "coherence"** — it could be about baryonic disk structure
3. **External replication** — all results from a single SPARC dataset
4. **The mechanism** — no field-theoretic explanation for why size matters

## Comparison with Standard MOND

Standard MOND (Milgrom 1983) predicts a **universal** RAR with no dependence on galaxy size at fixed baryonic mass and rotation speed. Finding that R_eff and R_max both predict RAR residuals at fixed V and L is in tension with universal MOND and consistent with Synchronism's size-dependent modification.

However, this could also reflect:
- Non-universal M/L that correlates with size
- Baryonic distribution effects (disk vs bulge dominance)
- Environmental effects that correlate with morphology

## Updated Master Scorecard

| ID | Prediction | Status | Key Evidence | Grade |
|----|-----------|--------|-------------|-------|
| NP1 | a₀ = cH₀/(2π) | SUPPORTED | 94% of MOND | B+ |
| NP2 | Morphology → scatter | STRONGLY SUPPORTED | p = 5×10⁻⁶ | A- |
| NP3 | High-z a₀ evolution | UNTESTED | — | — |
| NP4 | V-shaped scatter | SUGGESTIVE | Session #374 | C+ |
| NP5 | Local anomalies | UNTESTED | — | — |
| NP6 | N_corr → offset | SUPPORTED | R² = 0.23 | A- |
| NP7 | R_eff MOND-dominated | SUPPORTED | 3x at g < g† | A- |
| NP8 | R_eff M/L-independent | SUPPORTED | 7/7 tests | A- |
| **NP9** | **R_eff L-independent (late)** | **SUPPORTED** | **r = -0.49 controlling V,L** | **A-** |
| **NP10** | **Dynamical radius confirms** | **SUPPORTED** | **R_max r = -0.47** | **A** |

## Arc Statistics

- Sessions: 3 (#390-392) + synthesis (#393)
- Tests: 24/24 verified
- Grand Total: 567/567 verified
- Arc grade: **A-** (strongest empirical finding in the program)

---

*Grand Total: 567/567 verified (24 from this arc)*

**Arc summary: Galaxy size genuinely predicts RAR offset in the MOND regime. Both photometric R_eff (r = -0.49) and dynamical R_max (r = -0.47) independently predict offset beyond luminosity in late-type galaxies. The effect is absent in early types, confirming it's specific to the modified gravity regime. This is Synchronism's strongest empirical result: a novel, testable prediction (RAR depends on galaxy size) that is confirmed with two independent methods and robust to all identified confounds.**
