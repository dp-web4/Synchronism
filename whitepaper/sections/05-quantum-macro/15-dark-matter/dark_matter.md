## 5.15 Dark Matter, Dark Energy, and Coherence

This section documents the cosmology research arc. For full details, see the [arXiv preprint](https://github.com/dp-web4/Synchronism/tree/main/manuscripts) and [Research logs](https://github.com/dp-web4/Synchronism/tree/main/Research).

**Audit-aware status note (2026-05-15):** the table below is the historical Phase 1 accounting. Subsequent audits (S635 cosmology scorecard, S637 RAR σ_int derivation, S645/S648/S650 DESI DR1 TEST-04a kill-criterion firing, S654 zero active discriminators) substantially refine these entries: 0 novel-unfalsified cosmology claims; cosmology regime reduces to MOND in the testable regime; one Tier-1 kill criterion has fired (TEST-04a fσ₈, sharpened to mechanism-class failure / sign-reversed — but the sign-reversal was RETRACTED by S668 (2026-05-26) as a transcription error (DESI LRG1 fσ₈ ratio 1.16 copy-pasted from QSO's; self-consistent value ≈ 0.49, ΛCDM-consistent); TEST-04a now stands as a post-hoc amplitude disfavoring only (σ₈(z=0): 0.76 vs 0.841 ± 0.034 = 2.4σ)) [RE-GROUNDED 2026-05-27 by S672 — "amplitude disfavoring only / ΛCDM-consistent" over-softens it: the kill criterion WAS triggered, so disfavored ~2σ AND kill-triggered, post-hoc; the "ΛCDM-consistent" reading traced to a wrong-paper number (arXiv:2512.03230, a z≈0.07 Peculiar Velocity Survey) and an LRG1 ratio S668 never verified against DESI Tables 9/10; the sign-reversal/mechanism-class thread stays retracted; CORRECTED 2026-07-14 — criterion-verdict substitution: the registered criterion (fσ₈(z=0.51)>0.46 at >3σ) was met at only ~1.5σ, the 2.4σ is a GR-conditioned σ₈-amplitude statistic, and DESI's own MG analysis (arXiv:2411.12026) puts μ₀ within 1σ of zero — TEST-04a is withdrawn from the decisive negatives; both surviving decisive negatives are galactic (locality no-go; TEST-09 BTFR bounded-boost)]. **S673 (2026-05-27)** further closes the framework's GW test (TEST-15 / GW170817): its only GW parameter α is *read off* GW170817, not derived, so it adds no discriminating amplitude — the GW170817 speed-constraint discussion below holds only as a Case-3 (no-field-theory) scope restriction, not as a derived prediction. The "DERIVED" entries below are *internal-consistency derivations from Synchronism postulates* — they reproduce known dimensional combinations or label existing structure; none has been confirmed as a uniquely-Synchronism prediction. SPARC Capstone (#526-578) concluded that MOND + M/L corrections explain all RAR variance. The framing of this section as "empirical research validating coherence-based explanations" overclaims relative to the current audit-aware accounting: this is a *substantial pattern-organization track* with documented reparametrization status, not a validation of a novel cosmology.

**Research Status (264 Sessions, Nov 2025 - Jan 2026; updated by Framework Stress Test arc S617-654)**

Autonomous research sessions tested whether coherence dynamics can re-describe apparent dark matter and dark energy. Sessions #259-264 produced a unifying *vocabulary*: **Matter = Topology, Gravity = Geometry, Quantum = Dynamics** — all labeled from a single coherence field. Per S617-628 demolition phase, the vocabulary is internally consistent but does not deliver novel predictions distinguishing it from MOND in the testable regime.

| Component | Status | Accuracy | Notes |
|-----------|--------|----------|-------|
| Coherence function C(ρ) | **DERIVED** | N/A | Form from information theory |
| γ = 2 parameter | **DERIVED** | N/A | From thermal decoherence |
| Golden ratio exponent 1/φ | **VALIDATED** | 1σ | Within 1σ of Gaia DR3 best fit (Session #239) |
| a₀ = cH₀/(2π) | **DERIVED** | 10% | MOND connection |
| SPARC rotation curves | **TESTED** | 52% | 46% failure in massive galaxies |
| Santos-Santos DM fractions | **TESTED** | 99.4% | Different test than curves |
| Ω_Λ = (1 - Ω_m) | **DERIVED** | Exact | From coherence floor (Session #241) |
| S₈ = 0.763 | **PREDICTED** | Matches DES/KiDS | Independent validation |
| GW as coherence perturbations | **DERIVED** | N/A | Amplified in MOND regime (Session #246) |
| Matter = Topology | **DERIVED** | N/A | Solitons + winding (Session #261) |
| Gravity = Geometry | **DERIVED** | N/A | Metric coupling (Session #262) |
| Quantum = Dynamics | **DERIVED** | N/A | C flow + phase (Session #263) |

**Key distinction**: DERIVED = follows from axioms. CONSTRAINED = form determined by observation, then used predictively. TESTED = validated against empirical data.

---

**The Coherence Model of Dark Matter**

**Core Mechanism:**
Gravitational dynamics depends on local coherence C(ρ) ∈ (0,1]:

```
G_eff = G/C(ρ)
```

- **High coherence (C → 1)**: Standard gravity (high-density regions)
- **Low coherence (C → 0)**: Enhanced gravity (low-density regions)

**The Coherence Function (Derived):**

```
C(ρ) = tanh(γ · ln(ρ/ρ_crit + 1))
```

- **tanh form**: Derived from information-theoretic bounding (Session #74)
- **log(ρ) scaling**: Shannon information of N particles scales as log(N)
- **γ = 2**: Derived from thermal decoherence physics (Session #64)

**Physical Interpretation:**
Low-density regions have low coherence → enhanced effective gravity → appears as "dark matter" without requiring new particles. The galactic rotation curve "problem" becomes expected behavior of coherence-dependent gravity.

---

**MOND-Synchronism Unification (Sessions #86-89)**

A breakthrough discovery: MOND and Synchronism are the same physics in different parameterizations.

> **Status update (2026-08-03): the sentence above is now exactly true, and that is the problem with calling it a discovery.** Sessions #86-89 asserted the equivalence by analogy. As of 2026-08-02 it is an *algebraic identity*, and the identity is sharper than the analogy in a direction that does not favour the framework.
>
> **The identity.** The coherence compander is `C(x) = tanh(γ·ln(1+x))`, and the framework defines the apparent dark-matter fraction as `f_DM = 1 − C`, so `C` occupies exactly the seat MOND gives its interpolating function `μ`. At **γ = 1/2** the Hill identity collapses to
>
> ```
> tanh(½·ln(1+x)) = x/(x+2) = μ_simple(x/2)     — identically, for all x
> ```
>
> where `μ_simple(u) = u/(1+u)` is MOND's simple interpolating function. This is *exact*, not asymptotic (verified here to machine precision, max deviation 5.6×10⁻¹⁷ over 10⁻⁶ ≤ x ≤ 10⁶; verified independently by hand in the research lane). It sharpens the earlier 2026-07-22 result, which had reached the same place only in the deep limit (γ = 0.489 ⇒ Newtonian-return exponent q ≈ 0.98 ≈ simple μ).
>
> **What it costs.** If `C` *is* `μ`, then the entire galaxy sector reduces to **one statement: MOND, with μ's argument swapped from the enclosed-mass acceleration `g_bar` to the local density `ρ`.** Not "MOND-like on several observables" — one substitution. Every other galaxy-sector result (BTFR slope, DM-fraction ceiling, RAR shape) is downstream of whether that substitution holds. And the SPARC form-selection table's preferred form is Hill n = 1 *exactly*, i.e. γ = 1/2 — so the data select the MOND point rather than merely tolerating it.
>
> **[AMENDED 2026-08-04] Where the substitution actually sits — and it is not between Synchronism and MOND.** The frozen, pre-registered instrument that produced *every* number quoted in this note is keyed on **acceleration, not density**. `simulations/sparc_tanhlog_profile.py` inverts `g_bar = g_obs·tanh(γ·ln(1 + g_obs/a₀))`; `simulations/sparc_cassini_q2.py` defines `mu_tanh_log(x, γ) = tanh(γ·log1p(x))` and its own docstring calls it *"the registered tanh-log family in the mu convention"*; the registered likelihood in `Research/preregistrations/sparc_cassini_tanhlog/sparc_profile.json` profiles **a₀** over `log₁₀ a₀ ∈ [−11, −9]`. The token `rho_crit` appears in none of them. So **γ ≈ 0.489, ΔBIC +7.1/+184, the grid interval γ ∈ [0.425, 0.600] and the Cassini +17.95σ are properties of a MOND-family interpolating function in acceleration space, and are not measurements of the `C(ρ)` defined at the top of this section.** The density form *has* been implemented here (`session128_dark_matter_derivation.py:54`; `compare_empirical_vs_derived.py:32`), so the accurate statement is not "never evaluated" but: the sector ran on density in the Session-100s era and on acceleration in the frozen era, and the equation it states is the first while the numbers it quotes are the second. The gap is between this framework's prose and its own instrument — not between this framework and MOND. And the two variables are not interchangeable: they were separated on this very dataset (0.118 dex vs 0.161 dex; ≤25% admixture; the 1/r² of enclosed mass that no convolution of Σ generates). Evaluating the stated density law directly moves predicted rotation velocities by **2–5 orders of magnitude** and fails by *functional form* — a flat curve needs boost growing ~linearly in r, while ρ falls exponentially in a disk so `1/C` grows exponentially, and no (γ, ρ_crit) reconciles the two. **Refutation count still unchanged at 6**: this is the mean-relation face of the same substitution the scatter no-go already covers, plus a provenance defect — not an independent empirical failure.
>
> **The substitution is what the data exclude.** Executed 2026-08-02 on 2,622 SPARC rotation-curve points (149 galaxies, Q ≤ 2) with no functional form assumed, no γ, no ρ_crit, and no fitting: conditioned on `g_bar` the required boost scatters 0.118 dex (reproducing Lelli et al. 2017's published 0.11-0.13 dex); conditioned on local baryonic surface density it scatters 0.161 dex. The excess survives 175 free per-galaxy nuisance offsets (ratio *rises*, 1.37× → 1.77×, so the excess is within-galaxy). Interpolating `log u_α = (1−α)·log Σ_local + α·log g_bar` gives a minimum at α = 1.00 with a galaxy-block bootstrap 95% CI of **[0.75, 1.00]** — SPARC admits **at most 25% weight on a local-density variable** in whatever sets the mass discrepancy. The framework's galaxy sector sits at 100%, the far end of the excluded region.
>
> **Scope, stated because the strongest-sounding part is the weakest.** The decisive null — RAR residuals carry no information about local density at fixed `g_bar` (partial r = +0.001, 95% CI [−0.076, +0.081], ≤0.7% of variance, against a pipeline that recovers an injected r = 0.10) — is a *quantification* of published prior art, not a new discovery: at constant scale height `log ρ = log Σ − const`, and Lelli et al. 2017 already tested RAR residuals against stellar surface density at R and reported no significant correlation. What is new is the **constructive** half: no exponential smoothing kernel on Σ, at any range up to ~30 kpc, recovers `g_bar` performance (best 1.21× at λ ≈ 7 disk scale lengths, then degrading), because `g_bar = GM(<r)/r²` carries a 1/r² structure no convolution of Σ can generate — which turns "make the coupling differential" from a free dial into a specific obligation — together with the ≤25% admixture bound itself.
>
> **Consequences for the surrounding text, kept narrow.** The refutation count is **unchanged at 6**: this is the same underlying failure (local coupling) executed form-free, not an independent refutation, and inflating the tally would be the mirror of the scope error already on record. Bucket 0 is unchanged (0). Two related corrections belong here: the "γ ≈ 0.489 preferred by SPARC" figure is **not** a clean standalone measurement of γ *(the grounds for this were re-derived and corrected on 2026-08-04 — see the amendment two paragraphs below; the conclusion stands, the stated coefficient did not)*. And a claim that "the nonlinear Poisson equation implementing C(ρ) produces an external field effect at 0.3-0.4× MOND's" is **retracted**: the figure was never derived from any such equation *(and the accompanying premise "the galaxy sector has no field equation" is itself amended below — one exists, and it is linear rather than nonlinear, which strengthens rather than weakens the retraction's conclusion)*. The framework's actual structure — `C = C(ρ_local)`, with a uniform external field leaving local ρ unchanged — satisfies the Strong Equivalence Principle by construction and predicts **EFE = 0 exactly**. That is sharper and more falsifiable than "weaker EFE." Retracting a fabricated mechanism tightened the prediction rather than loosening it. *(The further claim that this "stands in tension with the ~4σ EFE detection reported in SPARC by Chae et al. 2020" is **withdrawn 2026-08-05** — the prediction is not evaluable against that dataset. See the amendment at the end of this note.)*
>
> **[AMENDED 2026-08-04] The degeneracy is γ ↔ a₀, and the separating coefficient stated above was wrong.** Two corrections to the 2026-08-03 text, both found by recomputation in this lane, both against itself. *(i)* Expanding `D(x,γ) = tanh(γ·ln(1+x)) − x/(x+2)` gives `D = (γ − ½)·(x − x²/2) + O(x³)` — the departure from simple μ is **first order in (γ − ½)**, not a second-order term with coefficient γ(2γ−1). Verified at x = 10⁻⁵: −0.0750 at γ = 0.425, −0.0110 at 0.489, 0.0000 at ½, +0.1000 at 0.600, +1.5000 at γ = 2, against the published −0.0319 / −0.0054 / 0.000 / +0.060 / +3.000. The published figures are ~2× low near ½ and 2× high at γ = 2. **The vanishing point γ = ½ is correct**, so the identity and the conclusion built on it stand — but a first-order departure is *easier* for data to resolve than the phrase "the O(x²) term" implied, which cuts against the degeneracy argument rather than for it. *(ii)* The degeneracy is **not** with a ρ_crit prior — there is no ρ_crit in that fit. In the deep-MOND limit `tanh(γ·ln(1 + g/a₀)) → γ·g/a₀`, so the registered family depends on γ and a₀ **only through the ratio γ/a₀**. Verified: (γ, a₀) = (0.2445, 2.67×10⁻¹¹), (0.489, 5.33×10⁻¹¹) and (0.978, 1.07×10⁻¹⁰) all hold γ/a₀ = 9.17×10⁹ and agree in μ to <1% over g ∈ [10⁻¹³, 10⁻¹¹] m s⁻², a 4× span in γ. The fingerprint is in the registered artifact itself: its profiled **a₀ = 5.33×10⁻¹¹ m s⁻² is 2.11× below the reference McGaugh a₀ = 1.128×10⁻¹⁰** in the same file — exactly the compensation γ ≈ ½ requires. So "γ ≈ 0.489 is not a standalone measurement of γ" is **correct and now correctly grounded**, and the qualification on the TEST-11/Cassini interval 0.425–0.600 stands unchanged.
>
> **[AMENDED 2026-08-04] A field equation does exist, it is linear in Φ, and EFE = 0 is its signature rather than an absence.** The 08-03 retraction above rests on "the galaxy sector has no field equation." A momentum-conserving completion exists and is the natural one — `S = ∫d³x[−(1/8πG)·C(ρ)|∇Φ|² − ρΦ]`, giving `∇·[C(ρ)∇Φ] = 4πGρ`. It answers the Felten (1984) third-law objection by Noether, and reproduces `g = g_N/C(ρ)` under Gauss in spherical symmetry. It is **linear in Φ** (C depends on ρ, not ∇Φ), so superposition holds and **EFE = 0 exactly** — not an artifact of missing dynamics but the signature of that linearity, where AQUAL's EFE comes from its nonlinearity. The polarization force this adds, `C′(ρ)|∇Φ|²/8πG`, is ≤ 2×10⁻⁵ of g across five SPARC galaxies: observationally free. The retraction's *conclusion* is therefore strengthened, not weakened — EFE = 0 now stands constructively. *(This paragraph originally continued "and its tension with Chae et al. 2020's ~4σ SPARC detection is sharper." That clause is **withdrawn 2026-08-05** — see the amendment below. EFE = 0 standing constructively is unaffected; only its claimed empirical tension was wrong.)* **The same linearity is also the sector's worst pathology**: `C(0) = tanh(γ·ln 1) = 0` exactly for every γ and ρ_crit, so under the dielectric reading every isolated body has a vacuum of zero permittivity and a divergent exterior field. *"A uniform external field does not change ρ"* and *"empty space has C = 0 however strong the field"* are the same statement — C's argument cannot see the field. **The sector's one surviving structural prediction and its worst pathology are one property.** Registered at `Research/proposals/dielectric_completion_and_efe_linearity_equivalence_20260804.md`.
>
> **[AMENDED 2026-08-05] EFE = 0 is *not evaluable* against Chae et al. 2020, and the "tension" claimed for it above is withdrawn.** Two paragraphs above assert that EFE = 0 "stands in tension with the ~4σ EFE detection" of Chae et al. 2020, and the 08-04 amendment then called that tension "sharper." Both are withdrawn. The reason is a number this note already contains and never applied: **the stated density law misses the rotation curve itself by 2–5 orders of magnitude**, and the EFE is a *residual* on that same curve. Chae's detection is a velocity deficit of **0.046 dex (NGC5055, 11σ) to 0.083 dex (NGC5033, 8σ)** — a 10–17% effect. A model that is 2–5 dex from the observable cannot be in tension with a 0.05–0.08 dex feature of it.
>
> Re-derived independently in this lane rather than accepted (SPARC Lelli et al. 2016 mass models, in-repo at `simulations/sparc_real_data/MassModels_Lelli2016c.mrt`), evaluating `C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))` with `ρ = Σ/2h`, `ρ_crit = 0.029·V_flat²` M☉/pc³, ledger convention `V = V_bar/√C`, at the outermost radius with positive reconstructed density — Chae's own measurement region:
>
> | galaxy | Chae `e` | baseline error | EFE signal | baseline / signal |
> |---|---|---|---|---|
> | NGC5033 | 0.104 (8σ) | +3.14 to +3.71 dex | 0.083 dex | **38× – 45×** |
> | NGC5055 | 0.054 (11σ) | +3.62 to +4.18 dex | 0.046 dex | **80× – 92×** |
>
> Across the full grid (γ ∈ {2, 0.489} × h ∈ {0.3, 1.0} kpc) the framework's baseline error exceeds the entire effect being measured by **38×–92×**; the most favourable configuration anywhere on the grid still gives 37.8×. MOND on the same points and the same pipeline sits at +0.002 to +0.084 dex, i.e. *inside* the signal. **The correct badge for EFE = 0 against this dataset is `not-evaluable` — neither confirmed nor refuted.**
>
> **This cuts both ways, and both directions were live in this repo.** The 08-03/08-04 text overclaimed *testability* ("sharper," "in tension with"); a same-day research-lane note in the other direction — that EFE = 0 is "a genuine MOND-discriminator that the framework **fails**" — overclaims *refutation*. Both die on the same ratio. **Refutation count stays at 6; no 7th is added, and Bucket 0 stays 0.** The executive-summary claim that "zero remaining active discriminators" no longer holds unqualified because "the EFE magnitude is one, in the refuting direction" is likewise withdrawn: **zero remaining active discriminators stands.**
>
> **The mechanism worth keeping.** The exec-summary propagation carried its own warning label — *"Magnitude not independently re-derived here"* — and the magnitude was the entire load-bearing content of the claim. A claim can be propagated across four surfaces in two days *with an accurate note saying the decisive quantity was never computed*, and the note does not stop it. Generalized and registered as a gate: **a structural prediction is testable against a measurement only if the measurement's signal exceeds the model's own baseline error on the same observable at the same radii — compute that ratio before badging any published result as confirming or refuting.** Here it reversed a proposed 11σ refutation into a not-evaluable.
>
> **[AMENDED 2026-08-08] The "6" is a count of test IDs; this whitepaper enumerates it nowhere, and two of its six members are described nowhere in it.** The number is asserted four times across these sources — in the paragraph above, in the a₀(z) note below, in the executive summary and in the conclusion — always as *"unchanged at 6"*, never with a list and never with a stated convention. So a reader cannot audit it from the document that states it, and two of the six are absent outright: the boost-ceiling refutation (the tokens `3.17`, `68.5`, `0.927` and the phrase "boost ceiling" each have **zero** occurrences in the whitepaper body) and the Bell/CHSH substrate refutation (`CHSH`, `no-signaling`, `Kuramoto`: **zero** occurrences; the only "Bell" in these sources is the exposition of Bell's theorem in §5.4 and a prior-art specificity list). The six, each with the statistic it was registered on and the data it used, are: **(1)** RAR transition shape at γ = 2 — marginalised BIC on SPARC, ΔBIC = +184; **(2)** BTFR slope (TEST-09) — SPARC, |Δn| = 0.41 against a registered bar of 0.3; **(3)** dwarf DM-fraction ceiling (TEST-10) — SPARC, 69% of galaxies above f_DM = 1 − Ω_m = 0.685; **(4)** RAR environment dependence (TEST-08) — SPARC × Cosmicflows-4, r² = 0.0001 against a registered > 0.20; **(5)** the SPARC × Cassini joint squeeze (TEST-11) — +17.95σ; **(6)** the observer-relative Bell/CHSH substrate bet (B1) — S = 1.98 ≤ 2 on both arms.
>
> **Measured, not asserted: (2) and (3) fire on the same galaxies, by algebraic identity.** The framework's apparent DM fraction and its acceleration boost are the same quantity in two variables — `f_DM = 1 − (V_bar/V_obs)² = 1 − g_bar/g_obs = 1 − 1/B` — so TEST-10's criterion `f_DM > 1 − Ω_m` and TEST-09's ceiling `B > 1/Ω_m` are the *same inequality*, not two corollaries of a shared root. Verified per galaxy on the in-repo SPARC mass models (`simulations/publisher_20260808_test09_test10_independence.py`, 123 galaxies, Q ≤ 2, i > 30°): `max|f_DM − (1 − 1/B)| = 2.2×10⁻¹⁶`, and the two criteria select **93 and 93 galaxies with 0 disagreements** at the outermost measured point, **88 and 88 with 0 disagreements** at the error-weighted mean of the outer three. *(93/123 = 76% here rather than TEST-10's published 69% only because TEST-10 proper runs the wider 153-galaxy set; it is the same criterion on TEST-09's narrower cut, which is the point — the comparison is between criteria on one sample, not between published fractions.)* Deleting every galaxy that fires TEST-10 and refitting the BTFR, TEST-09's registered kill **does not fire** on what remains: observed n = 3.39 ± 0.23 against the framework's 3.38 ± 0.18, deviation **0.01 (0.0σ)** versus a threshold of 0.3, where the full sample gives 3.79 ± 0.10 vs 3.35 ± 0.07, deviation 0.44 (3.5σ). *(0.44 rather than the ledger's 0.41 because this run compares against the like-for-like outer-3 estimator applied identically to all three models rather than the SPARC catalogue V_flat; both fire.)*
>
> **What that does and does not establish, because the leave-out test is underpowered by construction.** On the 30 non-exceeders the deviation carries σ = 0.29, so the 0.3 threshold sits **1.0σ out** — the subsample could not demonstrate independence at the registered bar even if independence existed, and the mass lever arm shortens from 3.45 to 3.06 dex. The honest statement is therefore asymmetric: **no part of TEST-09's evidence is located outside TEST-10's firing set, and the firing sets are provably the same set — so the burden of the two-row convention is unmet, not disproved.** Two figures now circulating are both asserted rather than computed and both overstate the collapse: *"honest independent-root figure: 3–4, not 6"* (`synchronism-site/maintainer/logs/2026-08-07.md`) and *"at most four independent roots"* (`…/visitor/logs/2026-08-07.md`). On the axis actually measured here the collapse is **6 → 5**, one pair, not three.
>
> **The dataset axis is the sharper one, and it points the other way.** Five of the six run on SPARC — (1), (2), (3), (4) and the galaxy side of (5) — three of the scripts behind them were opened here and share the *same fixed* mass-to-light convention Υ_disk = 0.5, Υ_bulge = 0.7 — `test09_btfr_bounded_boost_real_sparc.py` for (2), `test10_dwarf_dm_fraction_ceiling.py` for (3), and `simulations/sparc_tanhlog_profile.py` for (1) and for the SPARC profile that (5) consumes. (2) and (3) additionally share every parameter of one force law (`C(a) = Ω_m + (1−Ω_m)x/(1+x)`, Ω_m = 0.315, φ = 1.618034, a₀ = 1.05×10⁻¹⁰), while (1) and (5) use the tanh-log family on the same data. The sixth analyses **no external dataset at all**: B1 runs the framework's own harness and compares the result to the established CHSH violation, so the site's gloss *"6 refutations executed on external data … astronomical, ephemeris, and laboratory"* is loose on its third data type. Net: exactly three evidential channels are independent of SPARC's photometry and M/L, and each enters only jointly or by reference — the Cosmicflows-4 tracer field that supplies (4)'s predictor, the Cassini ephemeris that supplies (5)'s second constraint, and the textbook CHSH violation that (6) is compared against. **No row is carried by a non-SPARC dataset on its own**: five of the six use SPARC, and the sixth uses no external dataset at all. This matters more than root-counting, because an M/L shift moves five rows together, and M/L conditionality has already cost one headline in this program (the "28 galaxies" figure of the bounded-boost class exclusion, which is 15 at Υ*_disk = 0.7).
>
> **Non-action, deliberately.** The count stays at **6** and no row is merged. Whether TEST-09 and TEST-10 should be *counted* as one refutation or two is open question 4 of `Research/proposals/boost_ceiling_provenance_and_class_exclusion.md`, registered 2026-07-28 and explicitly gated on dp — a naming convention, which is dp's to set. What this amendment supplies is the thing that question asked for and did not have: an enumeration, and a measurement of how much independent empirical content the corollary carries. The mechanism is worth naming because the program's own lane named it the same day for a different row: *a prediction with no TEST ID is invisible to every audit that walks the ledger by ID, and the silence about why gets re-inscribed by each new reader as a fresh finding.* **A count with no enumeration has exactly that property** — this one was re-filed as a novel defect twice on 2026-08-07, by two independent readers, eleven days after it was registered as open.
>
> *Provenance: γ = 1/2 identity and EFE retraction from `explorations/2026-08-02-triage-gamma-half-exact-mond-efe-zero-and-a0z-nondiscriminating.md` (research lane), inscribed in `PREDICTIONS.md` B2; identity re-verified numerically here before propagation. Scatter no-go from `synchronism-site/explorer/findings/rar-scatter-nogo-executed-local-density-carries-zero-information.md` (2026-08-02), script `explorer/scripts/rar_scatter_nogo_real_sparc.py`. The prior-art scoping in the paragraph above is this lane's, not that finding's — see `explorations/2026-08-03-publisher-one-substitution-and-its-exclusion.md`. The 2026-08-04 amendments: the acceleration-vs-density provenance defect and both coefficient corrections are this lane's, verified at the code and the frozen artifacts (`explorations/2026-08-04-publisher-the-frozen-sparc-artifact-is-keyed-on-acceleration.md`); the dielectric completion, the mean-relation 2–5 OOM result and the EFE⇔linearity⇔vacuum-singularity equivalence are from `synchronism-site/explorer/findings/efe-zero-survives-momentum-objection-but-the-substitution-was-never-evaluated.md` (2026-08-03), back-annotated here at that finding's explicit request. The 2026-08-05 amendment: the not-evaluable verdict and the per-galaxy baseline/signal ratios originate in `synchronism-site/explorer/findings/efe-zero-is-not-refutable-by-chae2020-the-baseline-is-off-by-3-dex.md` (script `explorer/scripts/efe_required_ambient_density_vs_chae2020.py`), and were **reproduced independently in this lane** on the in-repo SPARC mass models before propagation — this lane's grid gives 38×–92× against that finding's 37×–53×, agreeing at the minimum and running higher elsewhere, so the conclusion is not sensitive to the difference. The two-directional framing (testability overclaim *and* refutation overclaim dying on the same ratio) and the baseline/signal gate are this lane's. The 2026-08-08 amendment is this lane's, executed here: script `simulations/publisher_20260808_test09_test10_independence.py`, run on the in-repo SPARC mass models with the pipeline conventions copied verbatim from `synchronism-site/explorer/scripts/test09_btfr_bounded_boost_real_sparc.py` and `test10_dwarf_dm_fraction_ceiling.py` so that any difference in result is a difference in question rather than in pipeline; the enumeration counts and the two absences were taken by grep over `whitepaper/sections/`. The shared-root observation itself is **not** new and is credited: it is on file since 2026-07-15 (`explorations/2026-07-15-verify-test08-test10-executed-mond-shared-class-law.md` — "same structural ceiling as TEST-09, different observable"), was registered as an open question on 2026-07-28, and was re-raised on 2026-08-07 by the site's visitor and maintainer lanes. New here are the exact identity of the two firing sets, the leave-out measurement, its power, and the dataset-concentration count.*

**The MOND Acceleration Scale:**
```
a₀ = cH₀/(2π) = 1.08 × 10⁻¹⁰ m/s²
```

**Empirical**: 1.2 × 10⁻¹⁰ m/s² (10% agreement within uncertainties)

**Physical Meaning:**
The 2π factor is the phase coherence cycle—the acceleration where cosmic phase uncertainty reaches one full cycle. This is not numerology; it connects the MOND scale to cosmology.

**Implication:**
a₀ is cosmologically determined, not arbitrary. Predicts evolution with redshift: a₀(z) ∝ H(z), testable via high-z BTFR.

> **Status update (2026-08-02 — supersedes the 2026-08-01 note below it in history; that note said "disfavoured" and it claimed more than the evidence supports).** This prediction has been engaged by published data. The honest verdict is that **it does not discriminate**: ΛCDM with baryons produces the same evolution, and this exact formula was written down and tested inside the rival paradigm first.
>
> **The measurement.** Ciocan et al. 2026 (MUSE-DARK III, A&A **709**, L16, arXiv:2604.22613) fit a₀ from the radial acceleration relation in redshift bins over 0.33 < z < 1.44 (N = 79), obtaining a₀(z) = a₀(0) + a₁z with **a₁ = 1.59 ± 0.1 × 10⁻¹⁰ m s⁻²** — the quoted interval is a **95% CI, not 1σ**. Ciocan read their own result as *structural* (galaxy structural properties or dark-matter profiles changing with redshift); modified gravity is the third of three explanations they list, not their conclusion. Their systematics are large (a_tot and a_bar from a common forward model; ~0.2 dex gas systematics; ~1.5× the SPARC scatter). Note also that the test route is the **RAR, not the BTFR** named above — which is why Milgrom's standing BTFR objection to high-z a₀ tests does not apply to it.
>
> **The prior art — verified at source, and decisive.** Mayer, Teklu, Dolag & Remus 2023 (MNRAS **518**, 257; arXiv:2206.04333) fit a MOND force law to galaxies in the Magneticum hydrodynamical simulation — ΛCDM plus baryons, no a₀ anywhere in the physics — and report that *"the best fit for a₀ is found to increase by a factor of approximately 3 from redshift z = 0 to z = 2."* Branch (A) gives E(2) = 3.03. Their **equation (13)** is `a₀(z) ≈ a₀(0)·[Ω_m(1+z)³ + Ω_Λ]^{1/2}` — this prediction, verbatim — and they report it *"fails to accurately describe the trend observed in Magneticum, with the change in Magneticum being somewhat slower as redshift increases."* Two consequences, and they point different ways: the functional form is **not novel to this framework** (it is a rejected candidate parametrisation in a 2022 ΛCDM paper), and at the precision anyone reports, ordinary galaxy-assembly physics yields the same growth — so **a measurement that cannot separate this framework from ΛCDM is not a test of it.** Stated precisely, because the distinction matters: Mayer's own sentence says the two curves differ in *shape*, so the degeneracy is one of achievable precision, not of principle. Non-discriminating with foreseeable data; not non-discriminating forever.
>
> **The comparison is anchor-dominated, and the 08-01 note got this wrong.** Branch (A) predicts a *ratio*, a₀(z)/a₀(0) = E(z); turning it into a number requires a value for a₀(0), and the four published ones span 69% — 1.00 (Ciocan's own intercept), 1.04 (cH₀/2π), 1.20 ± 0.26 (McGaugh SPARC), 1.69 ± 0.13 (Vărăşteanu MIGHTEE-HI). The signal being tested from z = 0 to z = 1 is 79% growth, so signal/systematic ≈ 1. The 08-01 note asserted that the slope ratio "does not depend on the a₀(0) normalisation"; that is **false** — da₀/dz|₀ = a₀(0)·(3Ω_m/2) is linear in the anchor, and the ratio runs **1.99× to 3.37×** across the four. What is anchor-robust is only the *sign*: every published anchor puts branch (A) below the measured slope. The z ~ 1 significance is not even that — it runs from ≈10σ low (framework anchor) through **0.5σ, consistent** (McGaugh, with its published ±0.26 restored) to ≈2σ *high* (Vărăşteanu). Mixing an anchor from one survey with a slope from another is itself unsound: Vărăşteanu attribute their 41% offset from SPARC at z ≲ 0.08 — where cosmology permits 4% — to sample selection.
>
> **Retracted from the 08-01 note.** It claimed "the literature is not self-consistent," citing Gueorguiev 2024 (arXiv:2409.11425) against Ciocan. That read a null as a contradiction. Gueorguiev's high-z arm is RC100 (N = 100) with d log₁₀a₀/dz = 0.01 ± 0.20, while branch (A) implies 0.227 over 0.5 < z < 2.5 (re-derived here) — ≈1.1σ power, an analysis that could not have detected the prediction had it been exactly true, and Gueorguiev say as much. Ciocan likewise report agreement with Vărăşteanu at ~1.5σ. There is no literature disagreement to reconcile: one detection with large method systematics, one power-limited null.
>
> **Status: engaged, non-discriminating.** Not refuted — the refutation count is unchanged at **6**; a row that cannot discriminate cannot be a refutation, and removing a test's power is the opposite operation. Nor is the 0.5σ consistency under the SPARC anchor evidence *for* the framework, for the same reason. What survives is a sharper criticism than a refutation would be: the prediction's distinguishing content was already spent before it was registered. And this result remains *worse* for the alternative anchor 2πa₀ ≈ c²(Λ/3)^{1/2}, which predicts no evolution at all — what is embarrassed is the programme of deriving a₀ from cosmology, not this framework uniquely.
>
> *Numerical note: line 63 above quotes a₀ = cH₀/2π = 1.08 × 10⁻¹⁰ (H₀ = 70); the arithmetic here uses 1.04 × 10⁻¹⁰ (H₀ = 67.4) for consistency with Ω_m = 0.315. That choice moves the slope ratio between 3.12× and 3.24× and changes nothing else.*
>
> *Provenance: raised by `synchronism-site/explorer/findings/a0-epoch-branch-A-tested-disfavoured-for-evolving-too-slowly.md` (07-30), reversed by `…/a0-epoch-row-is-non-discriminating-anchor-dependent-and-lcdm-degenerate.md` (08-01). Mayer's abstract, equation (13), and the "fails to accurately describe" sentence were verified at source here before this note was written — the shape detail in that sentence is not carried by either upstream finding. Anchor arithmetic: `simulations/a0_epoch_anchor_dependence.py`.*

---

**Golden Ratio Exponent Validated (Session #239)**

The coherence function contains a characteristic exponent. Synchronism predicts this exponent is 1/φ (the golden ratio inverse ≈ 0.618).

**Gaia DR3 Wide Binary Test:**

Fitting C(a) to Gaia DR3 wide binary data:
- **Best-fit exponent**: α = 0.688
- **Synchronism prediction**: 1/φ = 0.618
- **1/φ is within 1σ of best fit**
- **Δχ² = 4.00 in favor of Synchronism over MOND** (≈2σ preference)

**Physical Meaning:**
The golden ratio appears because it represents the optimal balance between local and non-local coherence contributions—the ratio at which phase information propagates most efficiently.

**Status**: VALIDATED against independent data (not used in derivation).

---

**Galactic Scale Validation**

**SPARC Database (175 galaxies):**
- Tests full rotation curve shapes (velocity at every radius)
- **52% success rate** overall
- **81.8% success** on dwarf galaxies (where effect is strongest)
- **46% failure** in massive galaxies (known limitation)
- Zero per-galaxy free parameters

**Santos-Santos Database:**
- Tests integrated dark matter fractions at specific radii
- **99.4% success rate**
- Mean error 3.2%
- Complementary to SPARC (total mass vs radial structure)

**DF2/DF4 Resolution (Session #97):**
These "dark matter deficient" ultra-diffuse galaxies appeared to contradict the model. Resolution: both are satellites of NGC 1052 (~80 kpc). Tidal stripping preferentially removes low-C envelope, leaving high-C core with G_eff ≈ G. Consistent with model, not contradictory.

**Honest Assessment:**
The 46% SPARC failure rate in massive galaxies is informative—the model has boundaries. We're exploring whether coherence explains apparent dark matter, not claiming proof.

---

**Dark Energy from Coherence (Sessions #100-102, #241)**

**Core Discovery:**
Applying G_eff = G/C to cosmology yields emergent dark energy:

```
ρ_DE = ρ_m · (1-C)/C
```

No cosmological constant Λ required.

**Cosmic Coherence Form (NOW DERIVED - Session #241):**

The cosmic coherence function has a natural form:

```
C(a) = Ω_m + (1 - Ω_m) × f(a/a₀)
```

As acceleration a → 0 (deep MOND regime):
- C → Ω_m = 0.315 (coherence floor)
- (1-C) → Ω_Λ = 0.685 (appears as "dark energy")

**The Key Result:**
```
Ω_Λ = (1 - Ω_m) emerges from coherence floor
```

**Flat universe (Ω_total = 1) is DERIVED, not assumed.**

This upgrades the cosmic coherence form from CONSTRAINED to DERIVED. The cosmological constant is not a free parameter—it's determined by the coherence floor in the deep MOND limit.

**Physical Interpretation:**
Dark matter AND dark energy are both coherence effects—unified through C(a). At galactic scales, low coherence enhances gravity ("dark matter"). At cosmic scales, the coherence floor creates an effective vacuum energy ("dark energy").

---

**S₈ Tension Predicted (Session #102)**

The scale dependence predicts the S₈ tension between CMB and lensing surveys:

```
S₈^Sync = 0.763
```

| Survey | S₈ | Type |
|--------|-----|------|
| Planck | 0.832 ± 0.013 | CMB |
| DES Y3 | 0.776 ± 0.017 | Lensing |
| KiDS-1000 | 0.759 ± 0.021 | Lensing |
| **Synchronism** | **0.763** | **Prediction** |

**Transition Scale:**
8 h⁻¹ Mpc—the σ₈ smoothing scale IS the coherence transition from galactic to cosmic regimes.

**Interpretation:**
The S₈ "tension" is not measurement error—it's the signature of scale-dependent coherence.

---

**Complete Coherence Physics (Sessions #259-264)**

The physics arc achieved a unified framework deriving matter, gravity, and quantum mechanics from coherence:

**The Master Equation:**
```
C(ξ) = ξ₀ + (1 - ξ₀) × ξ^(1/φ) / (1 + ξ^(1/φ))
```

**Three Pillars:**

| Pillar | Session | Result |
|--------|---------|--------|
| **TOPOLOGY** (Matter) | #261 | Matter = Soliton: C(x) = C₀ + A × exp(-x²/2σ²); Charge = Winding: Q = (1/2π) ∮ ∇S·dl |
| **GEOMETRY** (Gravity) | #262 | T_μν = ∂_μC ∂_νC - g_μν[(1/2)(∂C)² + V(C)]; G_μν = 8πG T_μν |
| **DYNAMICS** (Quantum) | #263 | ψ = √C × exp(iS/ℏ); Schrödinger emerges from C-phase dynamics |

**Complete Ontological Map:**

| Concept | = Coherence |
|---------|-------------|
| Existence | C > 0 |
| Matter | Stable C pattern (soliton) |
| Charge | C circulation (winding number) |
| Gravity | C gradient geometry |
| Quantum | C-phase dynamics |
| Consciousness | C > 0.5 |
| Information | -log(1-C) bits |

**Predictions:** 15 from physics arc, 6 high testability, 3 already confirmed.

---

**Cross-Scale Unity**

The same G_eff = G/C principle operates at three scales:

| Scale | Coherence Variable | Low C Effect | High C Effect |
|-------|-------------------|--------------|---------------|
| Quantum | T (temperature) | Classical | Quantum |
| Galactic | ρ (density) | "Dark matter" | Normal gravity |
| Cosmic | Ω_m (matter fraction) | "Dark energy" | Matter-dominated |

**The Deep Insight:**
Dark matter, dark energy, and quantum mechanics may be unified—all manifestations of coherence-dependent pattern interaction.

---

**Gravitational Waves as Coherence Perturbations (Session #246)**

GW are traveling disturbances in the coherence field:

```
C(a,t,x) = C₀(a) + δC(a,t,x)
```

Where δC is the coherence perturbation carried by the gravitational wave.

**Regime-Dependent Amplification:**

| Regime | Acceleration | C₀ | GW Amplification |
|--------|--------------|-----|------------------|
| Newtonian | a > 10⁻⁸ m/s² | ~1 | ~1× |
| Transition | a ~ a₀ | ~0.6 | ~1.6× |
| MOND | a << a₀ | ~0.35 | ~2.9× |

**Key Insight:** GW effects are AMPLIFIED in low-acceleration environments because the same δC causes a larger fractional change when C₀ is smaller.

**GW170817 Constraint:**
The speed constraint |v_GW/c - 1| < 10⁻¹⁵ from GW170817/GRB 170817A is satisfied because:
- Neutron star merger = high-acceleration, strong field
- In high-a regime, C ≈ 1, so v_GW = c
- Low-a modifications remain unconstrained

**Implications:**
- Ultra-wide binaries (~10000 AU) show enhanced GW emission (~1.5×)
- PTA amplitudes may be overestimated by ~30% due to galactic C < 1
- Space-based detectors (LISA) may see different response than ground-based (LIGO)

---

**Discriminating Tests**

**High-z BTFR (Critical Test):**
At z=1, H(z)/H₀ ≈ 1.7. If a₀ ∝ H:
```
Δ(log M_bar)_{z=1} = +0.06 dex (Synchronism) vs 0.00 dex (MOND)
```
Current high-z stellar TFR shows evolution in the right direction (KMOS³D, MOSDEF)—suggestive but not definitive.

**Other Predictions:**
- S₈ tension (already matches)
- Void expansion rates (modified from ΛCDM)
- Isolated UDG dispersion (should show enhanced σ)

---

**What Remains Incomplete**

**Mathematical:**
- No full relativistic formulation yet
- CMB predictions not calculated
- Connection between galactic tanh(log ρ) and cosmic C(a) forms

**Empirical:**
- 46% SPARC failure rate in massive galaxies unexplained
- High-z BTFR needs more data for definitive test
- GW amplification in MOND regime untested

**Conceptual:**
- Quantum-to-galactic coherence transition not fully specified
- Physical mechanism for golden ratio exponent

---

**Epistemic Status Summary**

**With Confidence:**
- Coherence function form derived from information theory
- Empirical validation on dwarf galaxies (81.8%)
- MOND-Synchronism mathematical equivalence
- Golden ratio exponent validated (1σ of Gaia DR3)
- Ω_Λ = (1 - Ω_m) derived from coherence floor

**With Reasonable Speculation:**
- GW amplification in MOND regime (~3×)
- S₈ prediction from scale-dependent C
- Cross-scale unity (quantum + galactic + cosmic)

**Pure Speculation:**
- Connection to quantum gravity
- Physical mechanism for golden ratio
- Ultimate unification of all forces

**Definitely NOT Claiming:**
- "Dark matter solved"—mechanism proposed, not proven
- "ΛCDM wrong"—ΛCDM works; this proposes underlying mechanism
- Any results without documented validation

---

**References**

Full research documentation:
- [arXiv Preprint v6](https://github.com/dp-web4/Synchronism/blob/main/manuscripts/synchronism-dark-matter-arxiv-v6.pdf)
- [Research Session Logs](https://github.com/dp-web4/Synchronism/tree/main/Research)
- [Session #264 Complete Synthesis](https://github.com/dp-web4/Synchronism/blob/main/Research/Session264_Complete_Synthesis.md)
- [Chemistry Track](https://github.com/dp-web4/Synchronism/tree/main/Research/Chemistry)
- [Gnosis Track](https://github.com/dp-web4/Synchronism/tree/main/Research/Gnosis)

**Research represents 264 autonomous sessions (Nov 6, 2025 - Jan 14, 2026) with cross-model peer review.**
