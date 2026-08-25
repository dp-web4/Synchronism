# The argument of C: three coherence functions, three claims of novelty, one incommensurable ledger

**Date**: 2026-08-24 · **Origin**: synchronism-site explorer track · **Status**: EXECUTED on real SPARC
**Site finding**: `explorer/findings/the-argument-of-C-three-functions-each-killed-by-its-own-distinguishing-feature.md`
**Scripts**: `explorer/findings/scripts/two_pillars_{argument_of_C,head_to_head_fit,neff_and_epsilon}.py`
**Supersedes the enumeration in**: `three-C-problem-existential-ambiguity.md` (2026-03-20)

---

## Claim

The galaxy sector runs **three different coherence functions**, each carrying a *different* claim to
being the framework's structural difference from MOND. **No function has more than one of those
features**, and each is refuted by exactly the feature that distinguishes it.

| | C_ρ (headline) | C_g (the numbers) | C_Ω (the ceiling) |
|---|---|---|---|
| form | `tanh(γ ln(1+ρ/ρ_crit))` | `tanh(γ ln(1+g_obs/a₀))`, implicit | `Ω_m+(1−Ω_m)x/(1+x)`, `x=(g_bar/a₀)^(1/φ)` |
| novelty feature | **EFE = 0** | **return exponent q = 2γ** | **bounded boost B ≤ 1/Ω_m** |
| refuted by | that feature (null lever; ΔBIC +2843) | that feature (γ=2 ⇒ RAR; γ≈½ ⇒ Cassini) | that feature (BTFR n→2; f_DM ≤ 0.685) |

Provenance of the keying, verbatim from this repository:

- `simulations/sparc_tanhlog_profile.py:84` — inverts `g_bar = g_obs·tanh(γ·ln(1+g_obs/a₀))`
- `simulations/sparc_cassini_q2.py:43` — `mu_tanh_log`, *"registered tanh-log family in the mu convention"*
- `simulations/session131..152` — the bounded form, run in **ρ**: `x = (rho/rho_t)^(1/φ)`
- the site's TEST-09/10 scripts run the bounded form in **g_bar** — a fourth cell

The 2026-08-04 exploration recorded the acceleration keying and drew one conclusion (γ ≈ 0.489 is a
property of MOND's μ). This is the second conclusion, which it did not draw.

## What is new

1. **The implicitness that produces flat rotation curves is the field-dependence that destroys
   EFE = 0.** Measured deep-limit slopes `dlog g_obs/dlog g_bar`: explicit **+0.0014** (⇒ `v ∝ √r`,
   never flat), implicit **+0.5091** (⇒ MOND's √ law, flat). Appendix D's completion
   `∇·[C(ρ)∇Φ] = 4πGρ` is linear in Φ *because* C is keyed to ρ; keying it to `|∇Φ|` makes it AQUAL.
   There is no implicit version of C_ρ — ρ does not depend on `g_obs`.

2. **At γ = ½ the EFE is MOND's identically.** `tanh(½ln(1+x)) = x/(x+2) = μ_simple(x/2)`, so the
   algebraic-μ EFE factor `μ(x)[1+L(x)]` matches MOND-simple at `x/2` to **2.2×10⁻¹⁶** over six
   decades. The γ=½≡MOND identity was known for the *fit* since 2026-08-02; extending it to the EFE
   removes the last structural difference.

3. **The a₀ factor-2 is that identity.** The identity forces `a₀ = g†/2 = 6.0×10⁻¹¹`; an independent
   head-to-head fit here returns **6.06×10⁻¹¹** (1.0%). `2×a₀_frozen` vs `g† = 1.20±0.24 (sys)` is
   **0.56σ**. Retire it as an anomaly.

4. **C_ρ has now been fitted head-to-head** — same points, same likelihood, 3 parameters each, so
   ΔBIC = Δ(−2lnL). N = 2438, 122 galaxies.

   | model | γ | σ_int | ΔBIC vs C_g |
   |---|---|---|---|
   | C_g | 0.5093 | 0.1215 dex | — |
   | C_ρ, global ρ_crit | 0.0462 | 0.2261 dex | **+2843** |
   | C_ρ, `ρ_crit = A·V_flat²` | 0.0389 | 0.2496 dex | **+3309** |

   Sign holds across 6 (h, ϒ) conventions (+2843 … +3582). N_eff-deflated: **+142 / +166**, against
   the site's headline +184 → +11.5. **γ_ρ → 0.04 means the likelihood is switching the density
   dependence off** — an independent route to the 2026-08-19 "local density carries ≤0.7% of RAR
   variance" result. σ_int of 0.23–0.25 dex is ~2× the RAR's *total* observed scatter (0.13 dex).

5. **ε = 2γ − 1 is the galaxy sector's only deformation from MOND, and SPARC does not constrain it.**
   Galaxy-level bootstrap (120×): **γ = 0.545 ± 0.129 ⇒ ε = +0.090 ± 0.258 (0.35σ)**. The quoted
   γ = 0.489 ± 0.02 is the same interval inflated ~4.5× by treating 20 points/galaxy as independent.
   This is the same ε the 2026-08-11 dark-energy work found controls the ΛCDM deviation.

6. **The ledger of six is not commensurable.** #1,#2 → C_Ω; #3,#5 → C_g; #4 → C_ρ; #6 → QM sector.
   **The largest coherent sub-ledger is 2.**

## Consequences for the archive

- **Appendix D §D.2's `∇²Φ = 4πGρ/C`** is a C_ρ statement. If the galaxy sector's numbers are C_g's,
  Appendix D does not describe the model that was fitted. The 2026-08-08 L1/L2/L3 fork analysis
  enumerated only division-family laws in ρ; the implicit-in-g branch was not in it.
- **The DF2 repairs** (`manuscripts/arXiv_preprint_draft_v1.md` §5.2 formation coherence; Session #97
  tidal stripping) both modify C away from `C(ρ_local)`. ~~which the 2026-08-23 note observed
  invalidates the EFE = 0 premise.~~ **[CORRECTED 2026-08-25, Publisher lane]** The cited warrant is
  refuted and the conclusion survives on this finding's own reason. The 2026-08-23 claim — that any
  non-local C-contribution invalidates `C local ⇒ EFE = 0` — was refuted by computation on
  2026-08-24 (commit `7f0d208a`, `simulations/efe_locality_vs_phi_dependence.py`), 4h49m before this
  proposal was filed: the derivation turns on **Φ-independence**, not locality, and a fully non-local
  but Φ-independent C leaves the EFE at 5.6×10⁻¹³ while a ∇Φ-keyed C gives 4.6×10⁻². Moving C off
  `ρ_local` therefore does *not*, by itself, cost EFE = 0. **But this document's own Result 1 supplies
  the reason that does**: the fits key C to `g_obs`, and keying C to `g_obs` *is* Φ-dependence. So the
  sentence that follows stands, and stands on a stronger footing than the one it cited —
  the *fits* invalidated EFE = 0 five months earlier, by field-dependence rather than by non-locality.
  *(Defect shape: verified conclusion, refuted warrant. A grep cannot find it, because the conclusion
  is correct; only reading the citation catches it. Flagged independently by the Archivist, 2026-08-25.)*
- **Sessions 131–152's bounded-ρ form** is a fourth cell with no site presence and no test attached.
  Either retire it or state what it is for.

## Not claimed

- C_g being MOND is not a success: MOND-simple is excluded by Cassini at +17.95σ, and that constraint
  is Desmond, Hees & Famaey 2024's.
- The EFE used is the algebraic-μ linearisation `B_∥ = 1/(μ[1+L])`. The γ=½ identity is algebra on C
  and is unaffected; the *size* of the O(ε) residual would change under a full AQUAL solve.
- C_ρ could be defended by stipulating that `g_obs` in the fit means the internal field only. That is
  a new, unstated, non-covariant premise and would mean the fitted model is not the stated one.

## Recommended archive action

Register the three functions as distinct models with distinct IDs, and attach every existing galaxy
result to one of them. Until that is done, "what does the framework predict for X" has no answer,
which is what `three-C-problem-existential-ambiguity.md` said on 2026-03-20 and what five months of
non-propagation has left unchanged.
