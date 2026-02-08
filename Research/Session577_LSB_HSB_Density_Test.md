# Session #577: LSB vs HSB at Fixed V_flat — The Clean Density Test

**Date**: 2026-02-08
**Status**: 8/8 verified

## Overview

Session #576 found a significant point-level density signal (r_partial=+0.193, p<10⁻²⁶) but couldn't distinguish physical density effects from systematics. This session performs the cleaner test: compare LOW surface brightness (LSB) and HIGH surface brightness (HSB) galaxies at the same V_flat. If density matters, LSB galaxies should show systematically different MOND offsets.

## Central Result: Surface Brightness Adds NOTHING Beyond V and L

**Partial r(offset, logΣ | V, L) = +0.097 (p=0.264) — NOT significant.**

When controlling for both velocity and luminosity, surface brightness does not predict the MOND offset. The Session #576 point-level signal was driven by the galaxy-level V-L-SB correlation, not by density physics.

## Key Findings

### 1. LSB vs HSB Sample (Test 1)

Median SB: 165 L_sun/pc². LSB (67 gal, mean type 8.2) vs HSB (68 gal, mean type 4.1). LSB and HSB differ strongly in V_flat (t=-11.9, p<0.0001), which means raw comparisons are meaningless — need to control for velocity.

### 2. The Central Test: Offset vs SB (Test 2)

| Control | r_partial | p |
|---------|-----------|---|
| None (raw) | +0.028 | 0.746 |
| logV | **-0.216** | **0.012** |
| logV, logL | +0.097 | 0.264 |
| logV, log x_outer | -0.094 | 0.277 |
| 6-var model | -0.114 | 0.189 |

**The V-only partial is -0.216 (significant!)** — at fixed V, higher SB → lower offset. This is the correct direction for density-based physics (higher Σ → higher ρ → more coherent → less boost).

**BUT adding L kills it** — r drops from -0.216 to +0.097 (not significant). The V-only signal was driven by L: at fixed V, brighter galaxies have higher SB AND different offsets (through logL×f_gas).

### 3. Matched LSB-HSB Pairs (Test 3)

61 matched pairs (|ΔlogV| < 0.1):
- Mean Δoffset (LSB - HSB) = **-0.013 ± 0.026** (t=-0.48, p=0.63)
- r(Δoffset, ΔlogΣ) = +0.203 (p=0.12) — not significant

**No detectable difference** between LSB and HSB offsets at the same V_flat.

### 4. SB as 7th Variable (Test 4)

- 6-var LOO = 0.885
- 6-var + logΣ LOO = 0.885 (ΔLOO = **-0.001**)
- F-test: F=1.66, p=0.200

Surface brightness **worsens** the model slightly. No predictive value.

### 5. Within-Type Analysis (Test 5)

- Late types (n=91): partial r(offset, logΣ | V) = -0.140 (p=0.19)
- Early types (n=44): partial r(offset, logΣ | V) = +0.008 (p=0.96)

No significant SB effect within either morphological type.

### 6. SB and RC Diversity (Test 6)

**One striking finding**: partial r(c_V, logΣ | V, L) = **+0.471** (p<0.0001). HSB galaxies have more concentrated (declining) RCs. This is physically interesting — it means SB predicts rotation curve SHAPE — but it doesn't help predict offset because c_V is already in the 6-var model.

### 7. Size-Velocity Residuals (Test 7)

At fixed V:
- r(offset, logΣ | V) = -0.216
- r(offset, log R | V) = -0.263

R (which is 1/SB × L essentially) is a slightly better predictor than SB, but both collapse when L is added.

### 8. Synthesis (Test 8)

**The density signal from Session #576 was a mirage.** What appeared as "R carries information beyond g_bar" at point level was actually "galaxies with different V and L have different offsets" — which the 6-var model already captures.

**Result**: Acceleration-based MOND is not challenged by SPARC surface brightness data. The RAR offset is determined by V, L, and their interactions (f_gas, c_V), not by surface density.

**What SB DOES predict**: RC shape (c_V). This is interesting physical information (denser centers → more concentrated RCs) but it's already incorporated in the model.

## The Session #576 Signal Explained

Session #576 found r_partial(offset, log R | g_bar) = +0.193 at point level. This survives because:
1. Within a galaxy, R varies while g_bar also varies → the two are correlated
2. Across galaxies, different V and L → different R → different offset
3. Controlling only for g_bar (not V and L) leaves galaxy-identity confounded with R

When properly controlled (V, L, and other galaxy properties), R adds ΔLOO=+0.005 (Session #537) and SB adds ΔLOO=-0.001 (this session). The "density signal" was a confound, not physics.

## Implications for Synchronism

The density-based transition C(ρ) is not supported by SPARC:
- At fixed V and L, surface brightness (the best density proxy) is irrelevant for offset
- LSB and HSB galaxies at the same V_flat have indistinguishable offsets
- The acceleration-based MOND ν(g/a₀) captures all the physics

This closes another potential avenue for distinguishing Synchronism from MOND using SPARC data.

## Grade: A-

A well-designed follow-up to Session #576 that resolves the apparent density signal. The finding that SB strongly predicts c_V (r=+0.471) is physically interesting, though not directly relevant to the density-vs-acceleration question. The matched LSB-HSB comparison (p=0.63) is the cleanest test and gives a clear null result.

## Files Created

- `simulations/session577_lsb_hsb_density_test.py`: 8 tests
- `Research/Session577_LSB_HSB_Density_Test.md`: This document

---

*Session #577 verified: 8/8 tests passed*
*Grand Total: 1749/1749 verified*

**Key finding: Surface brightness adds NOTHING to offset prediction beyond V and L (partial r=+0.097, p=0.26; ΔLOO=-0.001). Matched LSB-HSB pairs at same V_flat: Δoffset=-0.013, p=0.63 (null). Session #576's point-level density signal was a galaxy-identity confound, not density physics. SB does predict RC shape (c_V, r=+0.47) but this is already in the 6-var model. MOND (acceleration-based) wins over Synchronism (density-based) for SPARC. Grade A-.**
