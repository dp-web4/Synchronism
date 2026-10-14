# Proposal — the Jeans-knee "10σ" is not a MOND test; RG Eq. 4.1 has two in-repo readings (2026-10-14)

*Site maintainer, from simulated-visitor graduate and researcher personas (same model family; internal review, not
outside review). No bucket moves. Count 6. Bucket 0 = 0.*

## 1. The ρ_crit(V) row's "MOND's V⁻² is excluded at ~10σ too" names the wrong object

The Bucket-2 row "C(ρ) galaxy threshold ρ_crit(V) — wrong exponent" (restated 2026-08-28) says SPARC gives
ρ_crit ∝ V^(−0.16 ± 0.19), so "the framework's V⁺² is excluded at ~11σ — but MOND's V⁻² is excluded at ~10σ too".

The measured quantity is the **Session 53 Jeans-construction knee**, ρ_crit ∝ V^(s−2) with s = dlog Σ_c / dlog V from
the forward regression of disc central surface density on V. MOND's transition density is a different object: the
mean enclosed density where g_bar = a₀, 3a₀/(4πG r_t) with r_t ∝ V², which scales as V⁻² **by construction** on any
sample that obeys the BTFR. The site's own pages say "mean-enclosed and isothermal estimates both give exactly −2".

So the 10σ excludes *a V⁻² law for the Jeans-construction knee*, i.e. it shows that the Jeans construction cannot be
the a₀ knee. It does not test MOND. Read literally, the current sentence claims SPARC excludes MOND's transition
scaling at 10σ on data MOND fits to a few percent. A referee would read it that way.

**Recommended ledger wording (row, bracketed 2026-08-28 text):** "…The framework's V⁺² is excluded at ~11σ. A V⁻²
law *for this knee estimator* is excluded at ~10σ too, so the Jeans construction cannot track a₀; MOND's own
transition density (mean enclosed, at g_bar = a₀) is a different quantity and scales as V⁻² by construction."
Verdict unchanged. Applied on the site today (/parameter-derivations, /honest-assessment, /for-researchers) and
added to the site lint.

Also on the site: the /for-researchers "one-line citable form" still read "inverted sign, 240×–300,000× magnitude
error", which the same 08-28 restatement withdrew. A withdrawn figure cannot be the citable form; the site now cites
the measured exponent.

## 2. "Identical to RG's permittivity" is true for C_Ω, not for the headline form

`simulations/publisher_20260827_rg_floor_is_not_universal.py` verified C_Ω ≡ RG's ε. That C_Ω is the **floored Hill
form** Ω + (1−Ω)·u^p/(1+u^p). The site's /coherence-function said "a permittivity identical to C_Ω in closed form"
directly under the headline tanh(γ ln(1+x)), and a persona correctly checked it against the headline. The headline
equals the Hill form only at γ = ½ (p = 1, knee doubled). At γ = 2, x = 1: 0.882 vs 0.059 (p = 4, knee doubled) or 0.5.
Small-x: linear vs x^p. Site corrected.

Consequence for the GC window (P611.2 execution, 09-07/09-08 blocks): every "density-keyed" row used the floored
**headline** form, `OMEGA_M + (1−OMEGA_M)·tanh(γ·log1p(ρ/ρ_c))`. The row "γ = 2, knee at RG's 0.0083 → passes" is
therefore *this framework's function at RG's knee*, not RG. The site row is relabelled.

## 3. RG's own form: two readings in-repo, and the GC verdict depends on which

The explorer's `gc_knee_bound.py::fine()` evaluated RG "with its own published permittivity" as
tanh((log₁₀ρ − log₁₀ρ_c)/Q). The archive's verbatim read of M&D Eq. (4.1) (publisher, 08-27) is
tanh[log (ρ/ρ_c)^q] = tanh(q·log(ρ/ρ_c)), base unstated. For Q = 0.47 the steepness per decade is 2.13 (/Q), 1.08
(Q·ln), 0.47 (Q·log₁₀).

Executed today (`synchronism-site/maintainer/scripts/gc_rg_permittivity_form_check.py` + `_output.txt`; predictions
in the header before running, all three held; identity control: the /Q row reproduces the published −0.070):

| RG reading (ε₀ = 0.089, Q = 0.47, ρ_c = 0.0083) | obs − pred | verdict |
|---|---|---|
| /Q (explorer, 09-07) | −0.070 | ok |
| Q·ln (verbatim, ln) | −0.083 (+1σ ρ_c: −0.101 marginal) | ok |
| Q·log₁₀ (verbatim, log₁₀) | −0.114 (−0.097 … −0.127 across rows) | marginal |

MOND+EFE = −0.093, Newton −0.057, stat ±0.027 per row. **RG is never excluded under any reading**, so the 09-08
"Refracted Gravity is not excluded" stands. "Passes" depends on the log base. The /Q form is not the published form
under either base. The differences are inside about 1σ stat, so nothing here is a verdict on RG.

Asks:
- (a) Fix or annotate `gc_knee_bound.py::fine()` (explorer lane): the /Q form is not Eq. 4.1.
- (b) Fix the log base by reading Cesare+2020/2022's equation verbatim. The 08-27 note already withholds the exponent
  comparison for this reason. That base decides whether this row reads ok or marginal.
- (c) **Before any of this is stated as a claim about RG**, the RG authors are the right check (persona's
  suggestion). Contacting them is outward-facing and **gates on the project lead**.

## 4. Two research directions the personas converged on (seeded as explorer topics, not proposals)

- **A general LLR bound on density-keyed G_eff.** The 09-24 (i) Ġ/G argument was presented as "68/74 grid laws". The
  researcher persona's general form: the Sun's vertical oscillation through the disc gives |dlnρ/dt| ~ 10⁻⁸ yr⁻¹ even
  at kpc smoothing, so LLR bounds |dlnG_eff/dlnρ| ≲ 10⁻⁵ at the solar position for *any* local-density-keyed G_eff.
  If that is not prior art, it is the most portable result in the record. Topic:
  `synchronism-site/explorer/topics/general-llr-bound-on-density-keyed-g.md`.
- **A shared-scale bet.** Both personas independently: every sector sits at γ = ½ with a separately fitted scale (Λ:
  ρ_DE = 2ρ_crit; galaxies: a₀ → 2a₀′). The only novelty left in the functional form would be a relation between
  ρ_Λ/2 and a₀′ that the framework fixes rather than borrows from Milgrom's cH₀/2π. Topic:
  `synchronism-site/explorer/topics/shared-scale-between-rho-lambda-and-a0.md`.

## 5. Not new, recorded

The post-cutoff A2ACW arm was asked for a sixth time, this time by a researcher persona as "the cheapest
experiment on the site". It already gates on the project lead (10-12 block); CBP-Claude's 10-12 confound (add a "what does this
reduce to?" arm) still applies.
