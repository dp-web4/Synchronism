# Session 633: C(ρ) One-Decade Saturation — Fourth Site Claim with Same Failure Mode

**Date**: 2026-04-25
**Type**: Targeted archive + analytical audit responding to 2026-04-24 back-annotation
**Grade**: A- (clean math, clean disposition)

---

## Trigger

`Research/proposals/coherence_function_saturation_one_decade.md` (filed 2026-04-24, surfaced via S632 commit). Visitor researcher used the public Coherence Explorer at γ=2 and observed:

```
C(ρ_crit)     = 0.8824
C(10·ρ_crit)  = 0.9999
C(100·ρ_crit) = 1.0000
```

The function saturates within ~1 decade of ρ. The public site states: *"Synchronism maps density to coherence across 80 orders of magnitude, from quantum to cosmic."* The proposal flags these as structurally incompatible.

## Numerical Verification

For C(ρ) = tanh(γ·log(ρ/ρ_crit + 1)), the transition window (C: 0.05 → 0.95) versus γ:

| γ | Decade-width of transition |
|---|----------------------------|
| 0.5 | 2.56 |
| 1.0 | 2.01 |
| 2.0 | 1.77 |
| 5.0 | 1.64 |
| 10.0 | 1.60 |

The transition asymptotically bounds at ~1.6 decades for any sharp γ. **Even an arbitrarily sharp γ cannot make the function vary across more than ~2 decades of ρ around its specific ρ_crit.** This is a structural property of tanh-of-log(1+x).

The proposal's claim is verified: meaningful variation in C spans ~2 decades, not 80.

## Critical Exponent Test (Interpretation 1)

The proposal's Interpretation 1 reframes C(ρ) as a phase-transition order parameter. Phase transitions have critical exponents (β, ν, η, …) characterizing power-law singularities at T_c. Mean-field gives β=1/2; 2D Ising β=1/8; 3D Heisenberg β≈0.367.

For C(ρ) = tanh(γ·log(ρ/ρ_crit + 1)), expand near ρ_crit (ε ≡ (ρ-ρ_crit)/ρ_crit):
- log(ρ/ρ_crit + 1) = log(2 + ε) = log(2) + ε/2 − ε²/8 + …
- C = tanh(γ·log(2) + γ·ε/2 + …) = C(ρ_crit) + (γ/2)·sech²(γ·log(2))·ε + O(ε²)

C is **analytic** in ε with regular Taylor expansion. There is no power-law singularity. The function has **no critical exponent** β in the phase-transition sense.

**Interpretation 1 fails analytically.** C(ρ) is a smooth crossover, not a phase transition. The MIPT analogy in the proposal cannot be sustained — MIPTs have actual singularities, while tanh-of-log is everywhere smooth. The proposal's connection to "BKT scaling on tree graphs" (cited 2026-04-12 explorer note) is closer to the truth: C(ρ) is too smooth to be either a thermodynamic phase transition or a MIPT.

## Archive Search

The "80 orders of magnitude" phrasing is **not in the archive**. It originates on the public site's homepage metadata (`<title>Synchronism | One Equation, Every Scale</title>` and `<meta name="description" content="...density to coherence across 80 orders of magnitude...">`).

The archive talks about "across all scales" in Hermetic-principle prose but does not make the specific quantitative interpolation claim. The site's framing introduces a precision the archive doesn't support.

## What the Function Actually Does (Interpretation 3, Made Honest)

ρ_crit = A·V_flat² is per-galaxy. γ = 2/√N_corr is per-system. So each system gets its own ρ_crit and its own γ. The function C(ρ) is **one functional form applied per-system with system-specific parameters**, not one curve spanning all densities.

The 80 orders of magnitude refers to the **range of ρ_crit values across systems** (quantum nuclear at ~10²⁸ kg/m³ down to galactic at ~10⁻²² kg/m³). It does NOT refer to the smoothness window of any single C(ρ) curve.

A genuinely scale-invariant C function would have a power-law argument like (ρ/ρ_crit)^α with α << 1, which would actually broaden the transition. tanh-of-log doesn't have that property; the +1 inside the log forces local behavior near ρ_crit.

## Verdict

| Claim | Status |
|-------|--------|
| C(ρ) saturates within 1-2 decades regardless of γ | ✅ Verified |
| Site's "smooth interpolation across 80 orders" framing | ❌ Misleading — conflates ρ_crit range across systems with smoothness window of one curve |
| Interpretation 1: C is a phase-transition order parameter | ❌ Refuted — function is analytic, no critical exponents |
| Interpretation 3: per-system ρ_crit | ✅ Operationally honest, but site doesn't communicate it |

## What This Adds

**Four site-claim audits, all same failure mode:**

| Site claim | Archive source | Failure |
|------------|---------------|---------|
| BTFR n ≈ 2.2 (S631) | Session #48 Track B | Self-labeled "not rigorous"; observed n = 3.85 ± 0.09 |
| A = 4π/(α²GR₀²) (S631) | Session #66 | α = 1.0 (fiducial), not fine-structure constant |
| TEST-07 λ ~ 500 Mpc (S632) | Session #4 Track C | Dimensionally inconsistent; cosmetic factor adjustments |
| C(ρ) "80 orders of magnitude" (S633) | Site metadata, not archive | Conflates ρ_crit-range with smoothness window |

The pattern is robust. Three external-feedback proposals in one week have all surfaced site-archive disconnects. Internal review never did.

## Recommended Site-Side Actions (Operator: Site Is Reference-Only)

1. **Replace homepage tagline** "across 80 orders of magnitude" with something like: *"Each system has its own coherence transition (ρ_crit, γ); the framework provides one functional form, applied per-system. The 80-orders-of-magnitude span refers to the range of ρ_crit across physical systems, not the smoothness window of any single curve."*
2. **Coherence Explorer page**: add a note that the visible saturation within 1 decade is a structural property of tanh, not a numerical artifact. State explicitly that C(ρ) is a smooth crossover, not a phase transition.
3. **Drop the MIPT analogy** in any framing that suggests C(ρ) has critical exponents — it doesn't, by direct Taylor expansion.

## Connection to Prior Findings

- 2026-04-11 explorer note: "C(ρ) is a mean-field caricature of MIPTs." This session strengthens that — even mean-field has β=1/2; tanh-of-log has β=∞ (no critical behavior). C(ρ) is below mean-field.
- 2026-04-12 explorer note: "fails on tree graphs (BKT scaling)." Consistent with this session's finding.
- 2026-03-30 explorer note: "wrong deep-MOND asymptotics, saturates too quickly for galaxies." This is the same saturation problem, observed at galactic scale.

The C(ρ) story is now closed: it's a per-system local crossover, mathematically smooth, with no scale-bridging interpolation property despite the site's framing.

## So What?

S631–633 establish that the archive contains evidence against four specific public-site claims. The disconnect is structural — internal sessions never cross-audited. External feedback closed it three times in three days.

For the operator: the audit channel is reproducible. If desired, the same shape applies to TEST-02, TEST-04, and any other Tier-1 prediction. Each takes ~30 minutes of focused work plus an archive search.

For my role: protocol is calibrated correctly. Silence on stale CBP firings; respond to genuine new content. Three productive audits in three days, with silence in between.

## Files

- `Research/Session633_Coherence_Saturation_Audit.md` (this document)
- Updates to SESSION_FOCUS.md
