# Proposal: retract the DATED SCOPE NEGATIVE — the framework has a dark-energy sector, and it forbids the DESI quadrant

**Origin**: synchronism-site explorer track, 2026-08-10.
**Full finding**: `synchronism-site/explorer/findings/w-eff-the-archive-has-a-dark-energy-sector-and-it-forbids-the-desi-quadrant.md`
**Executable**: `synchronism-site/explorer/findings/scripts/w_eff_from_C_rho_cosmic.py`
**Refutation count: UNCHANGED.** This retracts a scope negative and adds a prospective registration.

---

## 1. What is being retracted

`PREDICTIONS.md:152` — **📌 DATED SCOPE NEGATIVE (registered 2026-07-22)**:

> "nothing in the framework sources, modifies, or couples to the expansion history, so there is no
> mechanism producing w ≠ −1 or w(z) evolution. **Verified before registering** (grep of
> SPINE/FUNDAMENTALS/PREDICTIONS/STATUS): no DE machinery exists."

This is false, and was false when registered. `Research/Session100_Modified_Friedmann.md` (2025-12-08)
derives `H² = 8πGρ_m/(3C)`, identifies `ρ_DE = ρ_m(1−C)/C`, and tabulates a w(z).
`Research/Session101_Cosmic_Coherence.md` follows it. `Research/Session107_DESI_Forecasts.md`
(2025-12-10) issues bin-by-bin DESI predictions.

**The verification grepped the four compilation documents** — SPINE, FUNDAMENTALS, PREDICTIONS, STATUS
— **and never touched `Research/`.**

### Rule this yields (recommend adding to registration discipline)

> **A negative existence claim must be verified against the primary derivation layer
> (`Research/`, `manuscripts/`, `explorations/`). Compilation documents can establish that something
> IS present; they can never establish that something is ABSENT.**

This is the mirror of the 2026-08-08 Appendix-D finding ("grep `manuscripts/`, not just `Research/`").
Two instances in three days; the failure is the *layer*, not the directory.

**Second-order harm.** The registration states that "any later claim that coherence 'explains' w₀wₐ is
to be judged against this statement that no such coupling exists." It was written as insurance against
retro-fitting. As it stands it would cause a **correct** rediscovery of the archive's own December 2025
work to be scored as a retro-fit. The anti-overclaiming machinery produced an under-claim and then
protected it.

---

## 2. Erratum against Sessions #100 and #101

Both state `w_eff = −1 + (1/3) d(ln ρ_DE)/d(ln a)`. The continuity equation gives
`w = −1 − (1/3) d(ln ρ_DE)/d(ln a)`. Unit test: the published form returns **−2 for matter** and
**−7/3 for radiation** (correct: 0 and 1/3). It is right for Λ only because the derivative vanishes.

Worse, the tabulated `w_galactic` column follows *neither* formula — it equals
`(1/3) d ln ρ_DE/d ln a`, i.e. the stated expression **with the leading −1 dropped**. Two independent
errors.

Inputs reproduce exactly (γ = 2, x₀ = ρ₀/ρ_crit = 0.16738 fixed by C(0) = Ω_m), so the diagnosis is
unambiguous:

| z | published | correct |
|---|---|---|
| 0.1 | +0.32 | **−1.3186** |
| 0.5 | +0.73 | **−1.7329** |
| 1.0 | +1.37 | **−2.3690** |
| 2.0 | +2.28 | **−3.2788** |

Closed form: `w(z) = −γ(1+C)x/[C(1+x)]`, `x = x₀(1+z)³` (analytic ↔ numeric to 7×10⁻⁹).

**Consequence for Session #101.** Its entire premise — "w_eff > 0 contradicts w ≈ −1, therefore the
galactic C form is a category error at cosmic scales" — is an artifact. The true value is w(0) = −1.24.
And its replacement is not a replacement: SymPy confirms, exactly,

```
C(γ=½) = tanh(½ ln(1+x)) = x/(x+2)
C_galactic(γ=½) with C₀=Ω_m  ≡  Ω_m(1+z)³/(Ω_m(1+z)³+Ω_Λ)  ≡  Ω_m(z)  ≡  C_cosmic   [difference = 0]
```

**`C_cosmic` IS the galactic function at γ = 1/2.** Session #101's "derivation" back-solved
`d ln ρ_DE/d ln a = 0` for C and rediscovered a member of the family it was rejecting; its
"verification" table (w = −1.00 at all z) is the assumption read back. Recommend a dated erratum on
both sessions and withdrawal of the "category error" framing.

---

## 3. What replaces the scope negative — a live, disfavoured prediction

**No free calibration exists.** Ω_m ≡ 8πGρ_m0/(3H₀²) by definition and the modified Friedmann gives
H₀² = 8πGρ_m0/(3C₀), so **C₀ = Ω_m identically — forced**. γ is the only parameter. Two exact limits
(w monotone throughout):

```
z → ∞ :  w → −2γ           a → ∞ :  w → −1   (every γ, independent of x₀)
```

⇒ **sign(w₀ + 1) = sign(wₐ) for all γ.** `w = −1` is a one-sided attractor; **the model cannot produce
a phantom crossing.**

DESI DR2 (arXiv:2503.14738) prefers a crossing — w₀ > −1 *with* wₐ < 0 — in all four combinations
(BAO+CMB: −0.42 ± 0.21, −1.75 ± 0.58; +Pantheon+: −0.838 ± 0.055, −0.62; +Union3: −0.667 ± 0.088,
−1.09; +DESY5: −0.752 ± 0.057, −0.86). Scanning γ ∈ [0.05, 20]: **0 of 16 reach that quadrant.** Forcing
w₀ to match forces wₐ to **+0.21 … +0.32** — wrong sign in every case.

So the framework is **not silent** on the dominant live anomaly in cosmology. It **forbids** it.

---

## 4. Honest bounding (please preserve these caveats verbatim if this is propagated)

- **The naive H(z) comparison over-refutes by ~5.6×.** Fixed Ω_m and scale: −9.8% at z ≈ 0.72 vs ~1%
  BAO ("10σ"). **Marginalising Ω_m and r_d H₀** — the two nuisances BAO cannot separate from shape —
  gives rms shape residual **1.76% at γ = 2**, 0.97% at γ = 1, **0.031% at γ = 0.489**, **0.000% at
  γ = 1/2**. Real tension at γ = 2; *not* a kill. Reporting the fixed-parameter number would have been
  this program's fifth self-inflicted over-refutation.
- **The γ ≈ 1/2 branch is not refuted at all.** It *is* ΛCDM in the background and inherits ΛCDM's
  2.8–4.2σ DESI tension exactly — no more.
- **No quadrature "σ" is claimed.** w₀ and wₐ are anti-correlated at ρ ≈ −0.9; only the
  covariance-free *quadrant* statement is load-bearing.
- **The no-go is conditional on the substitution.** Session #100 substitutes G_eff = G/C into Friedmann;
  it does not derive the 00-component covariantly. A covariant Appendix-D completion generates Ċ terms
  absent here. **The theorem holds for the model as specified.**
- **In the framework's favour:** the source branch (L1, ∇²Φ = 4πGρ/C) was excluded a priori in the
  galaxy sector by a vacuum source floor (C → 0 as ρ → 0). That objection **does not apply to FRW**,
  where ρ̄ > 0 everywhere. Cosmology is the one arena where L1 is well-defined.

---

## 5. Two structural results worth carrying into the core

**(a) γ = 1/2 is a double-degeneracy point, and it is ONE algebraic fact, not two coincidences.**
γ = 1/2 is the unique member of the tanh-log family that is Möbius in x. Möbius C is exactly simple-μ
MOND in the galaxy sector (`C = x/(x+2) = μ_simple(x/2)`, on record) and exactly Λ in the cosmic sector
(proved here). The free-γ SPARC fit returns **0.489** — 2.2% away.
⚠ **But** the SPARC γ is fitted in *acceleration* space and this one lives in *density* space; they are
the same parameter only if the g_bar→ρ substitution is exact, which is open. Suggestive, not licensed.

**(b) The cosmic sector has the same architecture as the galaxy sector.** The galaxy sector is on record
as "MOND ∩ {B ≤ 3.17} — a strict submodel that could only tie or lose." The cosmic sector is a
one-parameter deformation of ΛCDM whose optimum *is* ΛCDM and whose every other member is worse. Two
sectors, one architecture. That is a statement about how the framework is built, not about any one test.

**(c) ρ_crit is unanchored across sectors by 1.5 × 10¹⁰.** Since C₀ = Ω_m is forced, ρ_crit,cosmic is
determined: 4.41 × 10⁻⁸ M⊙/pc³ at γ = 1/2, against NGC 3198's 652.3. The `ρ_crit = A V_flat²` law would
need V_flat ≈ **1.2 m/s** for the universe. ρ_crit is not a constant of the framework; it is a
per-system fitted quantity with no cross-sector relation. (Exact aside: ρ_crit,cosmic = ρ_Λ/2 at γ = ½.)

---

## 6. Requested decisions (gate on dp)

1. **Retract** the DATED SCOPE NEGATIVE at `PREDICTIONS.md:152`; keep date and reasoning as a record of
   the failure mode, replace the conclusion, and correct the "Verified before registering" line to name
   what was actually grepped.
2. **Adopt the verification rule** in §1 into registration discipline.
3. **Erratum** on Sessions #100 and #101 per §2; withdraw the "category error" framing.
4. **Adopt or decline** the prospective registration below.

> **PROSPECTIVE REGISTRATION (proposed, 2026-08-10, before DR3).** The source-branch cosmic sector
> requires **sign(w₀ + 1) = sign(wₐ)** and cannot produce a w = −1 crossing, for any γ or calibration.
> **Refuted if** a DESI DR3 / Stage-IV (w₀, wₐ) contour excludes the sign-locked locus at > 3σ using the
> **full published covariance**, with a single SNe compilation fixed in advance, and Ω_m and r_d H₀
> marginalised. **Nuisances marginalised**: Ω_m, r_d H₀. **Fixed**: flatness, no radiation, no massive-ν
> freedom. **Explicitly recorded**: this test can only kill or tie — γ = 1/2 *is* ΛCDM, so a confirming
> outcome is not discrimination and must not later be read as one. **Timeline**: ~2027–2028.

---

## 7. Highest-value follow-on

**Derive the 00-component covariantly from Appendix D and re-run.** It is the only identified route that
could move the locus off the sign lock. If the Ċ terms do not break it, the no-go becomes unconditional
and is a small, well-posed, potentially publishable result. Secondary: **Session #107's DESI forecasts
(3.1–3.2σ fσ₈) have never been audited** and predate the TEST-04a correction that found the registered
fσ₈ criterion met at only ~1.5σ — same over-forecast shape, unopened file.
