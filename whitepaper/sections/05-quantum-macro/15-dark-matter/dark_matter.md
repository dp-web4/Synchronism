## 5.15 Dark Matter, Dark Energy, and Coherence

This section documents the cosmology research arc. For full details, see the [arXiv preprint](https://github.com/dp-web4/Synchronism/tree/main/manuscripts) and [Research logs](https://github.com/dp-web4/Synchronism/tree/main/Research).

**Audit-aware status note (2026-05-15):** the table below is the historical Phase 1 accounting. Subsequent audits (S635 cosmology scorecard, S637 RAR σ_int derivation, S645/S648/S650 DESI DR1 TEST-04a kill-criterion firing, S654 zero active discriminators) substantially refine these entries: 0 novel-unfalsified cosmology claims; cosmology regime reduces to MOND in the testable regime; one Tier-1 kill criterion has fired (TEST-04a fσ₈, sharpened to mechanism-class failure / sign-reversed — but the sign-reversal was RETRACTED by S668 (2026-05-26) as a transcription error (DESI LRG1 fσ₈ ratio 1.16 copy-pasted from QSO's; self-consistent value ≈ 0.49, ΛCDM-consistent); TEST-04a now stands as a post-hoc amplitude disfavoring only (σ₈(z=0): 0.76 vs 0.841 ± 0.034 = 2.4σ)) [RE-GROUNDED 2026-05-27 by S672 — "amplitude disfavoring only / ΛCDM-consistent" over-softens it: the kill criterion WAS triggered, so disfavored ~2σ AND kill-triggered, post-hoc; the "ΛCDM-consistent" reading traced to a wrong-paper number (arXiv:2512.03230, a z≈0.07 Peculiar Velocity Survey) and an LRG1 ratio S668 never verified against DESI Tables 9/10; the sign-reversal/mechanism-class thread stays retracted; CORRECTED 2026-07-14 — criterion-verdict substitution: the registered criterion (fσ₈(z=0.51)>0.46 at >3σ) was met at only ~1.5σ, the 2.4σ is a GR-conditioned σ₈-amplitude statistic, and DESI's own MG analysis (arXiv:2411.12026) puts μ₀ within 1σ of zero — TEST-04a is withdrawn from the decisive negatives; both surviving decisive negatives are galactic (locality no-go; TEST-09 BTFR bounded-boost)]. **S673 (2026-05-27)** further closes the framework's GW test (TEST-15 / GW170817): its only GW parameter α is *read off* GW170817, not derived, so it adds no discriminating amplitude — the GW170817 speed-constraint discussion below holds only as a Case-3 (no-field-theory) scope restriction, not as a derived prediction. The "DERIVED" entries below are *internal-consistency derivations from Synchronism postulates* — they reproduce known dimensional combinations or label existing structure; none has been confirmed as a uniquely-Synchronism prediction. SPARC Capstone (#526-578) concluded that MOND + M/L corrections explain all RAR variance. The framing of this section as "empirical research validating coherence-based explanations" overclaims relative to the current audit-aware accounting: this is a *substantial pattern-organization track* with documented reparametrization status, not a validation of a novel cosmology.

**Research Status (264 Sessions, Nov 2025 - Jan 2026; updated by Framework Stress Test arc S617-654)**

Autonomous research sessions tested whether coherence dynamics can re-describe apparent dark matter and dark energy. Sessions #259-264 produced a unifying *vocabulary*: **Matter = Topology, Gravity = Geometry, Quantum = Dynamics** — all labeled from a single coherence field. Per S617-628 demolition phase, the vocabulary is internally consistent but does not deliver novel predictions distinguishing it from MOND in the testable regime.

| Component | Status | Accuracy | Notes |
|-----------|--------|----------|-------|
| Coherence function C(ρ) | **DERIVED** | N/A | Form from information theory |
| γ = 2 parameter | **DERIVED** | N/A | From thermal decoherence — **but not the value the galaxy sector runs at**; SPARC selects γ ≈ 0.489 ≈ ½ and fitting γ is out-of-sample *harmful* (2026-08-21, below) |
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
- **the `+ 1`**: a regulator, adopted so the log stays finite at ρ = 0 — and, as of 2026-08-21, the term known to be carrying the most physics. It does not merely regularize: it *creates* the deep-limit power law and pins its index to 1, which is the whole of why this equation reproduces MOND's deep limit. The equation names no principle that fixes that index.
- **γ = 2**: derived from thermal decoherence physics (Session #64). **The galaxy sector does not run at this value.** SPARC selects γ ≈ 0.489 ≈ ½ — the exact MOND point (identity below) — and out of sample, *freezing* γ at ½ predicts held-out galaxies **better** than fitting it.

**Physical Interpretation:**
Low-density regions have low coherence → enhanced effective gravity → appears as "dark matter" without requiring new particles. The galactic rotation curve "problem" becomes expected behavior of coherence-dependent gravity.

> **[ADDED 2026-08-21] The regulator is the load-bearing term, and the framework's one fitted parameter is out-of-sample harmful to fit.** Two results, one executed in the research lane (`synchronism-site/explorer`, 2026-08-20, re-run to completion in-repo the same day) and one verified symbolically here.
>
> **(i) What the `+ 1` does.** Give the density ratio a general index, `C_p(x) = tanh(γ·ln(1 + x^p))` — the equation above is `p ≡ 1`, and nothing in this whitepaper or the archive names that choice. Then `C_p → γ·x^p` as `x → 0`, for **every** `p` (verified here by symbolic limit: `C_p/x^p → γ` at p = 0.7, 1, 1.5). Delete the regulator and `tanh(γ·ln(x^p)) → −1`: a saturated constant, *no power law at all*. So the `+ 1` does not exclude power-law behaviour as ρ → 0 — it is what **creates** it and what pins the index. And the index is the whole physics: `p = 1` ⇔ asymptotically flat rotation curves ⇔ BTFR slope 4 ⇔ MOND's deep limit, which Milgrom (2009, ApJ **698**, 1630) *derives* from spacetime scale invariance. **This framework inherits the right exponent from a term chosen for finiteness.** Two results elsewhere in this whitepaper are `p ≡ 1` statements rather than structural ones: §6.4's strict concavity of `C` (an inflection reappears at p = 1.5, x = 0.54, and p = 2, x = 0.82 — reproduced here) and TEST-09's "deep limit is BTFR slope → 2".
>
> **(ii) The index was freed, fitted, and the result is null.** Real SPARC (Lelli+2016; 149 galaxies, 2,700 points), power-gated first by injection–recovery (`p̂` unbiased, worst bias 0.038, scatter ≈ 0.05). In-sample `p̂ = 0.762`, `Δ(2 lnL) = 15.0` — but the two best fits differ by **≤ 0.038 dex** across the SPARC range, under the 0.106 dex per-galaxy nuisance floor for one-parameter RAR corrections (arXiv:2608.08945), and 10-fold galaxy-level cross-validation returns **+0.33σ**. `p ≡ 1` stands *on execution, not on assumption* — the honest reason the galaxy sector does not beat MOND. Note the sign: a free `p` is an **extension** containing MOND at p = 1, the opposite of the bounded-boost ceiling, so "strict submodel" was never an a-priori argument.
>
> **(iii) Fitting γ is worse than not fitting it.** Nested cross-validation on the same data: freeing γ from MOND's ½ costs **−0.00159 ± 0.00119 = −1.34σ** in held-out per-point lnL (held-out lnL/pt 0.44040 frozen vs 0.43881 fitted). The mechanism, computed here in this section's own units: across the *entire* registered interval γ ∈ [0.425, 0.600], a single rescaling of `a₀` (×0.83 to ×1.22) reproduces the whole γ-family to **≤ 0.012 dex RMS and ≤ 0.022 dex maximum** in log₁₀ of the boost — an order of magnitude below SPARC's σ_int = 0.122 dex and below the 0.106 dex nuisance floor. So γ is not merely unidentified at a factor of two: over its own confidence interval it is a direction `a₀` already spans to a precision ten times finer than the data can resolve. Fitting it can only buy noise, and out of sample it does. **A parameter that hurts you when you fit it is a reparametrization with a knob** — this is the sharpest empirical form yet of this section's γ ↔ a₀ degeneracy result, and it is independent support for γ = ½ being the exact MOND point rather than a coincidence.
>
> **(iv) Why γ is unfittable, derived rather than fitted — and what deleting the `+ 1` would cost.** The regulator is what keeps γ *out* of the deep-limit exponent. With it, `C → γ·x` and the deep-limit solution of `g_bar = g·C(g/a₀)` is `g = √(a₀·g_bar/γ)`: **γ and a₀ enter only through the ratio a₀/γ, an exact degeneracy**, and the deep index is 1 for *every* γ — the framework gets BTFR slope 4 for free, from a term adopted for finiteness. That degeneracy predicts, with no fitting, that rescaling `a₀` by exactly `2γ` should absorb any γ; imposing that rescale blind reproduces the γ-family to **≤ 0.013 dex RMS** over the SPARC range, matching the fitted optimum above. γ is therefore not a shape parameter of this equation at all — it is a relabelling of `a₀`, which is the whole content of the −1.34σ. **Now delete the regulator**, as this archive's own standing proposal recommends (`Research/proposals/rho_crit_asymmetry_saturation_knee.md`, Option B, 2026-05-08: `C = (1 + tanh(γ·ln(ρ/ρ_crit)))/2`, adopted to make ρ_crit a true half-maximum). γ moves **into** the exponent: `C → x^{2γ}`, deep index `2γ`, BTFR slope `2(2γ+1)`. At γ = ½ that is 4 and nothing changes; at **this whitepaper's own derived γ = 2 it is BTFR slope 10**, against an observed 3.75 ± 0.10. That proposal's open Research Question 4 — *"If we drop the '+1' … is there any physical problem?"* — is hereby answered: yes, unless γ is pinned to exactly ½, in which case the framework's derived value is the one that has to go. Naively deleting the `+ 1` without recentering is worse still: `tanh(γ·ln x) → −1`, a negative saturated constant, outside `C ∈ (0,1]` entirely.
>
> **Scope, and it must travel with this.** The `p` direction is alive only in the **acceleration** reading `x = g/a₀` (the reading under which the `C ≡ μ` identity below actually holds). In the literal `x = ρ/ρ_crit` reading the 2026-08-02 form-free bound — σ(log B | ρ) = 0.161 dex vs σ(log B | g_bar) = 0.118 dex, an infimum over **all** functions of ρ — excludes every `p`. This is a one-parameter generalization *of* MOND's deep limit, not an escape from MOND's variable. **Refutation count unchanged at 6; Bucket 0 unchanged at 0.** Registering the executed `p`-null with a TEST-ID is a catalog action and is not done here.

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
> **[AMENDED 2026-08-16] The obligation stated in the paragraph above is DISCHARGED, and it closes the local branch entirely — the classification that had scoped the no-go away was wrong, and the sector's one distinctive prediction turns out to be the signature of the property that kills the fit.** Three legs, all re-verified in this lane before propagation. **(1) The differential class is finite and it is empty.** Buckingham-π on `F(ρ, |∇ρ|, ∇²ρ; G, a₀)` gives a rank-3 dimensional matrix on 5 quantities ⇒ a null space of exactly **2**, spanned by `x_diff = Gρ²/(a₀|∇ρ|)` and `q = ρ∇²ρ/|∇ρ|²`; `∇³ρ` adds exactly one further group (rank 3 on 6). Both re-derived symbolically here. So "which differential `F`?" is a **closed** question answerable by 2-D conditioning, not an open search — which is what makes this a closure rather than a spot-check. Enumerated on real SPARC (2,614 points, 145 galaxies) via the symmetry-integrated requirement `F_req = g_bar/g_obs` — no assumed form, no γ, no ρ_crit, no fitting; a theory is viable iff `F_req` is single-valued in its π groups — the conditional scatter is σ = **0.117 dex on `g_bar`**, 0.161 on local ρ, 0.195 on `x_diff`, 0.299 on `q`, 0.180 on the full 2-D class jointly, against a **0.309** no-information ceiling. Controlling for `g_bar`, **every group in the complete class explains ≤ 0.16% of the RAR residual** (|r| < 0.04; local ρ itself managed 0.7%), and `q` adds nothing to MOND — (`g_bar`, `q`) gives 0.1171 vs 0.1174 for `g_bar` alone and 0.1210 for the 2-D binning-cost control. **Every group in the complete local differential class carries less information about the missing physics than the algebraic variable it was proposed to replace.** The degeneracy escape that the 2026-07-28 kill criterion demanded genuinely exists here (`q` is 95.3% within-galaxy) — but it is *empty*, which is a stronger closure than one leaning on the exponential-disk idealisation. Note the contrast with the ϒ result recorded further below: where ϒ_disk sweeps γ̂ across [0.27, 0.96] at flat rms, this closure is ϒ-flat — the best differential group never comes within **1.34×** of `g_bar` anywhere in ϒ ∈ [0.30, 0.80]. A conditional comparison does not inherit the degeneracy that an estimated shape parameter does. **(2) The 2026-07-27 counterexample was misclassified, and the retraction it licensed is withdrawn.** Burrage, Copeland & Millington (PRD 95, 064050, 2017) is a symmetron: `∇²φ = (ρ/M² − μ²)φ + λφ³` is a nonlinear elliptic PDE, so φ(x) solves a global boundary-value problem and is a **non-local** functional of ρ — the screening that lets it track galaxy structure *is* that environment-dependence — and its closed form `g_sym = g_bar/(exp√(g_bar/g†) − 1)` is written in `g_bar = GM(<r)/r²`, an integral of ρ, which no local function of ρ(r) can equal. BCM never populated the local class. "Local" is hereby operationalized as **algebraic-pointwise**, resolving the fork this whitepaper opened on 2026-07-28 and referred out; the 2026-07-10 fabricated-consensus retraction still stands, and only the misclassification and the "the axis is algebraic-vs-differential, not locality" reframe are withdrawn. **So no local density coupling — algebraic or differential at any derivative order — reproduces the RAR; only non-local survives** (Milgrom's non-locality theorem), now backed by execution rather than by the consensus that was correctly retracted. **[SHARPENED 2026-08-19/20, THEN THE SHARPENING REFUTED 2026-08-22/23 — the base verdict (“only non-local survives”) is unchanged; the *stated mechanism* is now wrong for the second time in four days and is replaced by a weaker, range-conditional statement.]** The 08-19/20 reading recorded here was: the surviving branch is not non-locality in general but specifically **directional / inward-cumulative** non-locality, because on real SPARC a *symmetric* kernel `f(|r−r′|)` failed at any range (σ = 0.193 dex at λ = ∞, worse than pointwise ρ's 0.161) while the *causal* kernel `W(r,r′)·Θ(r−r′)` at the same range reached 1.02× `g_bar` with radial weight at Newton's `p = 1` — hence “the discriminating axis is **symmetry, not range**,” hence every natural repair of a pointwise multiplier lies in the dead branch. **That rule is refuted.** This paragraph named its own cheapest falsifier and marked the rule not-yet-citable; the falsifier was executed on 2026-08-22 (site explorer, `yukawa_symmetric_kernel_self_check.py`) and it returned the branch the hedge said would break the rule. **The refutation is analytic before it is empirical, and was re-derived in closed form in this lane:** the screened linear scalar `(∇²−m²)h = 4πGρ` has Green's function `e^(−mr)/r`, which is symmetric — two-sided and isotropic — at *every* screening length, and recovers Newton as `m → 0` (`F_Yukawa/F_Newton` = 0.736 / 0.995 / 1.000 at `λ = 1 / 10 / 10³` in units of `r`). The symmetric branch therefore *contains* `g_bar`, and no SPARC number is needed to see it. On SPARC it also measures out: the unscreened member passes the validation gate (σ = 0.1080 dex vs `g_bar`'s 0.1107, reconstruction cost −0.0026 dex), the symmetric field `g_Y` overlaps `g_bar` under the galaxy-block bootstrap at **every** `λ_s ≥ 1 R_d`, and at matched range it beats the causal kernel on every row of the head-to-head. Re-verified in this lane: the addendum bootstrap script re-run at **exit 0** with its tables reproduced value-for-value, and the Yukawa→Newton limit re-derived symbolically. **What this costs the section.** (a) “Every natural repair of a pointwise multiplier is in the dead branch” and “the live branch reconstructs precisely the variable `C(ρ)` was introduced to replace” are **withdrawn**; the **escape taxonomy reopens**, which is verbatim what the hedge pre-registered. (b) The 08-19 claim that the bootstrap “separates every finite λ ≤ 4 R_d from `g_bar`” **does not reproduce** — its CI lower edge was +0.0003 dex, and independent re-implementation gives +0.0075 [−0.0040, +0.0223] on the common-validity mask and +0.0096 [−0.0046, +0.0215] on 08-19's own unmasked point set, i.e. **OVERLAPS** both ways. (c) In exchange there is a genuine constructive result, and it is the one durable output of the execution: **RAR scatter bounds any screened-scalar escape to `λ_s ≳ 1 R_d ≈ 2.4 kpc`** (separated for `λ_s ≤ 0.75 R_d`, overlapping from `1 R_d` up). Finite-range non-locality is a *measured, physically-motivated* escape with a lower bound on its range — not a closed branch. **The replacement axis, and why this whitepaper carries it weaker than its source states it.** The source replaces “symmetry” with *which derivative of Φ keys the modification* — prior art (Joyce, Jain, Khoury & Trodden, Phys. Rep. **568**, 2015, already cited on the project site), with `∇Φ` (k-mouflage = MOND's acceleration) viable, `Φ` (chameleon/symmetron) failing, and `∇²Φ ∝ ρ` (Vainshtein/Galileon — nominally `C(ρ)`'s own rung) failing. The **measurement** is the archive's contribution and it is real; the **ladder** is not yet a ladder, and this was checked here against the source's own printed bootstrap rather than taken on report. Two of the three rungs *swap verdict* between the only two ranges bootstrapped: `|Φ_Y|` **overlaps** `g_bar` at `λ_s = 4 R_d` (+0.0171 [−0.0023, +0.0376]) and separates only at λ = ∞ (+0.0327 [+0.0118, +0.0578]), while the normalised `⟨Σ⟩_Y` **separates** at `4 R_d` (+0.0190 [+0.0046, +0.0330]) and **overlaps** at ∞ (+0.0123 [−0.0094, +0.0415]). Further, the `∇²Φ` rung's 1.40× is scored at `λ_s = 0` (pointwise Σ) while the other two rungs are scored at `λ_s = ∞`, so the three rows are not a matched comparison. **Only `∇Φ` is range-robust.** The accompanying explanation — that the 08-19 symmetric family failed because it was a *normalised* (intensive) smoothing of Σ — is likewise **not established**, and the routed proposal says so in its own words (“my going-in hypothesis is NOT established”) because the normalised member at λ = ∞ is unresolved from `g_bar`; that caveat is present in the proposal and absent from the triage note this repository committed. So the honest statement of the mechanism, as of 2026-08-23, is narrower than any of its three predecessors: **pointwise local density is separated from `g_bar` (1.40×, +0.0390 [+0.0179, +0.0592]) and every extended symmetric functional tested overlaps `g_bar` somewhere in range; the only range-robust discriminator measured is `∇Φ` versus the rest.** **What is untouched.** `C(ρ)` as Synchronism proposes it is *pointwise* — the separated case — so the section's verdict does not move; the 2026-08-15/16 Buckingham-π closure of the local branch (leg 1 above) rests on its own enumeration and is unaffected; the mutual-exclusivity tension in leg (3) below is unaffected. **And an over-refutation guard in the other direction, adopted from the 08-19 execution and still standing:** the “≤ 0.16%” and “≤ 0.7%” figures above are shares of the **residual after conditioning on `g_bar`**, which is the correct comparison and is how they are stated — but read *raw*, local Σ explains **73.0%** of the variance of log B against `g_bar`'s 85.9%. Local density is a poor organising variable, not noise. **Refutation count HELD at 6, Bucket 0 = 0** — a stated mechanism was corrected, no closure was added or removed. *(Sources: `synchronism-site/explorer/findings/causal-kernel-scan-the-required-nonlocality-is-a-direction-not-a-length.md`, 2026-08-19, routed via `Research/proposals/causal_kernel_required_nonlocality_is_directional_20260819.md` — **superseded**; and `synchronism-site/explorer/findings/yukawa-self-check-the-sorting-rule-is-refuted-the-axis-is-grad-phi.md`, 2026-08-22, routed via `Research/proposals/yukawa_selfcheck_locality_axis_is_grad_phi_20260822.md`, with triage at `explorations/2026-08-22-yukawa-self-check-refutes-my-08-19-symmetry-rule-the-axis-is-grad-phi.md`. The Yukawa scan is source-executed; `yukawa_addendum_bootstraps.py` was re-run in this lane at exit 0 and the Yukawa→Newton limit re-derived symbolically here.)* **(3) ⚠ The real output is a named tension: `{EFE = 0} ⟺ locality` and `{RAR-viable} ⟹ non-locality` are MUTUALLY EXCLUSIVE.** Everything above about EFE = 0 — that it holds exactly, that it is the signature of the field equation's linearity, that it is *not evaluable* against Chae et al. 2020 — is unchanged and is presupposed here. What is new is what it *costs*: EFE = 0 follows from `C = C(ρ_local)`, i.e. from locality, and locality is precisely what leg (1) closes. The framework can keep its one distinctive discriminator (local, refuted on the RAR) or buy the fit (non-local, and the known non-local reproducer is environment-dependent ⇒ EFE ≠ 0, so no distinctive prediction) — **not both**. **EFE = 0 is the observable signature of the very locality the RAR refutes.** This is the sector's sharpest statement to date and it is strictly stronger than "the galaxy sector reduces to MOND," because it identifies *which* property must be surrendered and *what* is lost with it. **Un-sealed edge, stated because it is the one well-posed survival question and is NOT a theorem:** no result forbids a non-local-but-SEP-preserving coupling — the known reproducers forfeit EFE = 0 only *generically*. That is the sector's open question, and it is now a single sharp one rather than a diffuse "reduces to MOND." **Refutation count HELD at 6** — a constructive lead closed, not an independent refutation, and inflating the tally would be the mirror of the scope error already on record — and **Bucket 0 = 0**: EFE = 0 remains a real prediction, welded to a form that cannot fit galaxies. The preprint statement and the formal scope-restoration **gate on dp**. *Provenance:* source `Research/proposals/differential_coupling_pi_enumeration_local_branch_closed_20260815.md` (site explorer track) and triage `explorations/2026-08-15-differential-local-branch-closed-bcm-is-nonlocal-07-27-misclassification-corrected.md`; script `synchronism-site/explorer/findings/scripts/differential_coupling_pi_enumeration_real_sparc.py`, **re-run in this lane at exit 0** with the quoted table reproduced, and both π-groups plus both null-space dimensions re-derived symbolically here. The 07-27→08-16 propagation gap is recorded in the executive summary and conclusion, where the withdrawn scoping had been live for 18 days after this lane posed — and deferred — the exact question that settles it.
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
> **[AMENDED 2026-08-09] The 08-04 paragraph above says "a field equation does exist" and "one exists." Two exist, they had been on file in this repository for seven months, and they disagree by a median 0.821 dex.** `manuscripts/Appendix_D_Synchronism_in_General_Relativistic_Form.md` §D.2 and §D.6.1, committed **2025-12-01** (`4400d54f`), state `∇²Φ = 4πGρ/C(ρ)` — call it **L1** — as the framework's Newtonian limit. That is nonlinear in the *source*, not linear in Φ, and it is not the equation this note introduced on 08-04. So the clause "one exists, and it is **linear** rather than nonlinear" is corrected to: *two readings exist; the archive's is nonlinear in the source and the one written here is linear in Φ; they are different laws and only one of them carries any evidence.*
>
> **The 08-04 completion is not a rival to the archive's — it is the archive's repair.** §D.2 states L1 and then reports, "in this limit", `g_obs = g_bar/C` (**L3**). L3 does not follow from L1. L3 is the spherical Gauss solution of exactly the equation this note wrote on 08-04, `∇·[C(ρ)∇Φ] = 4πGρ` (**L2**): integrate over a ball and `C·g·r² = G·M_bar(<r)`. So the archive already contained L2's phenomenology, attached to the wrong field equation, and the 08-04 amendment reconstructed the right equation without knowing the wrong one was on file. **L1 and L3 coincide iff C is spatially constant** — and C is the sector's entire content.
>
> **Measured here, not accepted** (`simulations/publisher_20260809_appendixD_coupling_fork.py`; in-repo SPARC mass models, 113 galaxies at Q ≤ 2, i > 30°; self-consistent spherical construction `ρ_sph = (4πr²)⁻¹ d[V_bar²r/G]/dr`, so no scale-height or disk-geometry convention enters and both laws see the same ρ): at the outermost measured point `log₁₀(g_L3/g_L1)` has **median +0.821 dex**, IQR [+0.515, +1.416], range [+0.023, +3.203]; **105/113 galaxies exceed 0.3 dex, 47/113 exceed 1 dex**. The constant-C control returns 2.6×10⁻¹⁵ dex. The gap is **γ-invariant to 0.0000 dex** between γ = 2 and γ = 0.489 — *and the reason bounds the fact*: every SPARC point sits at ρ/ρ_crit ≤ 0.045 (median 6.1×10⁻⁶), deep in C's linear regime where C → γρ/ρ_crit to within 2.4%, so 1/C ∝ 1/γ in **both** laws and γ cancels identically. The invariance holds *at this ρ_crit* and is a consequence of the calibration being far enough off that no galaxy reaches C's knee — the same fact as the 2–5 dex miss recorded above, not an independent robustness result. *(The sibling explorer lane reached median 0.81 dex on five galaxies the same day; this is an independent reproduction on 113 with a different density construction, and it adds the reason for the γ-invariance, which that lane reported as a bare fact.)*
>
> **Provenance, which is the publication-relevant part.** L1 has exactly one implementation in this repository — `simulations/session72_spherical_toy_model.py`, dated 2025-12-01 and written to complete §D.6 — and **zero citations anywhere in this whitepaper**. L3 appears in 26 simulation scripts and carries every number quoted in this section. L2 has **zero implementations**: it exists only as the prose two paragraphs above. So **the framework's only stated field equation has never been evaluated against data, and the force law that carries all of the evidence never had a field equation implemented for it.**
>
> **The eliminated reading is the better-fitting one, and saying so is the point.** Against the observed accelerations at the same radii, **L1 is closer than L2 = L3 in 113/113 galaxies** at both γ = 2 and γ = 0.489 (median offsets +4.64 vs +5.65 dex, and +5.25 vs +6.26 dex). Both are catastrophically wrong and neither is viable; the elimination below is **structural, not a goodness-of-fit verdict**, and a reader should not take the survivor to be the better-supported law. *(These absolute offsets are computed on the spherical `ρ_sph` construction and are **not** comparable to the "2–5 orders of magnitude" recorded above, which uses `ρ = Σ/2h`: spreading a disk's mass over spheres lowers ρ and raises the boost. The fork amplitude itself is a ratio of two laws on the same ρ and is robust to the choice — 0.821 dex here against the sibling lane's 0.81 dex on `Σ/2h`.)*
>
> **L1 is eliminated a priori, and the reason does not stop at L1.** `C → γρ/ρ_crit` as ρ → 0, so `ρ_eff = ρ/C → ρ_crit/γ`, a constant (verified to <10⁻⁴ relative for γ ∈ [0.25, 5] down to ρ/ρ_crit = 10⁻¹⁰): L1's source never vanishes in vacuum. The sharper objection, which the divergence obscures: **ρ_crit = A·V_flat^B is a per-galaxy constant**, so L1 assigns the same point of empty space a different source density for every galaxy one asks about — no single-valued source, i.e. ill-posed rather than merely divergent. **That objection survives into L2 = L3**, which the divergence argument does not, and it is the sector's one structural problem that no reading escapes. Contracting the same limit with `u_μu_ν` for dust carries it into §D.3: `𝒯_μν → (ρ_crit/γ)u_μu_ν` in vacuum — a nonzero stress–energy filling all space whose direction requires a preferred vacuum 4-velocity.
>
> **Scope, and what does not change. Refutation count stays at 6; Bucket 0 stays 0.** Nothing was refuted here: a *reading* was eliminated, which is bookkeeping, not evidence. What does change is a claim elsewhere in this whitepaper. The executive summary and conclusion both carry, from S642, *"there is no Lagrangian, no action principle, no equation of motion for C(ρ)."* Read against Appendix D that is three conjuncts of different truth value, and the sibling lane's flat verdict that "Appendix D contains all three" over-corrects in the other direction. **(a)** "No field equation" is **false** — §D.2/§D.6.1 state one. **(b)** "No action principle" is **defensible but must be scoped**: §D.5 writes an effective *worldline* action for a test particle's centre of mass, which is not a field action, and **no action generates L1**. *(§5.15's L2 is the one object here that does come from a field action — the one written in this note on 08-04.)* **(c)** "No equation of motion for C(ρ)" is **true and is the load-bearing conjunct**: C is algebraic in local ρ in every reading, with no kinetic term anywhere in Appendix D. So the GW170817 conclusion **survives on (c) alone** — no propagating scalar mode exists to modify tensor dispersion — and the correct statement is *"C has no kinetic term,"* not *"the framework has no action principle."* Back-annotated to the executive summary and conclusion rather than left at the source, because this lane has now twice measured that corrections do not propagate on their own.
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

**Cosmic Coherence Form (Session #241 — "DERIVED" corrected below, 2026-08-11):**

The cosmic coherence function has a natural form:

```
C(a) = Ω_m + (1 - Ω_m) × f(a/a₀)
```

As acceleration a → 0 (deep MOND regime):
- C → Ω_m = 0.315 (coherence floor)
- (1-C) → Ω_Λ = 0.685 (appears as "dark energy")

**The Key Result [lead corrected 2026-08-11]:**
```
Ω_Λ = (1 - Ω_m) — an identity of the calibration, not an emergence
```

Ω_m ≡ 8πGρ_m,0/(3H₀²) by definition, and the modified Friedmann H₀² = 8πGρ_m,0/(3C₀) forces **C₀ = Ω_m identically, for any coherence form** — Session #100 itself calls this calibration "a tautology." So Ω_Λ = 1 − Ω_m and flatness are built into the ρ_DE = ρ_m(1−C)/C split, not derived from it; the earlier "upgraded from CONSTRAINED to DERIVED" framing overclaimed. γ is the cosmic sector's only free parameter. (Sharpens the standing "internal-consistency reproduction" caveat; source: 2026-08-10 site-lane finding, `Research/proposals/dark_energy_sector_exists_and_forbids_desi_quadrant_20260810.md`; verified `simulations/publisher_20260811_w_eff_erratum_check.py`.)

**Physical Interpretation:**
Dark matter AND dark energy are both coherence effects—unified through C(a). At galactic scales, low coherence enhances gravity ("dark matter"). At cosmic scales, the coherence floor creates an effective vacuum energy ("dark energy").

**Erratum on the source sessions (2026-08-10; independently verified 2026-08-11):** Sessions #100/#101 state `w_eff = −1 + (1/3)·d(ln ρ_DE)/d(ln a)`; the continuity equation gives `w = −1 − (1/3)·d(ln ρ_DE)/d(ln a)` (the published form returns w = −2 for matter), and the tabulated w(z) column follows neither formula — it drops the leading −1. Corrected values: w(0) = −1.24 at γ = 2 (not > 0 as Session #100 reports), monotone with w → −1 as a → ∞ and w → −2γ as z → ∞. Session #101's "category error" conclusion (that the galactic C form cannot serve at cosmic scales) is an artifact of this sign error: at γ = 1/2 the galactic form with C₀ = Ω_m is *identically* Ω_m(z) — the exact ΛCDM background, w ≡ −1. Structural consequence, now **class-level** [corrected 2026-08-12; the 2026-08-11 form of this sentence conditioned the result on the absence of a covariant derivation of the 00-component — that derivation was executed the same day, site lane, and *hardens* the conclusion]: the substituted model obeys **sign(w₀ + 1) = sign(wₐ)** — no phantom crossing, which is the quadrant DESI DR2 prefers in all four of its dataset combinations (0 of 16 γ values reach it) — and the two minimal covariant completions of Appendix D §D.3's equation close the escape the substitution left open. **Completion A** (the equation as written, no new fields): the Bianchi identity forces ρ/C ∝ a⁻³, so the background is *exactly Einstein–de Sitter* for every γ and calibration — the dark-energy sector vanishes identically, and the vacuum floor ρ/C ≥ ρ_crit/γ leaves the equation no FRW solution beyond a ≈ 1.037 under Session #100's own calibration. **Completion B** (C as a Brans-Dicke scalar pinned to C(ρ̄(a))): the literal sign lock *dies* — the w = −1 attractor is destroyed, the γ = 1/2 Λ-degeneracy is broken (no member of the completed family is ΛCDM), and w₀ < −1 with wₐ > 0 members exist — but the DESI quadrant (w₀ > −1, wₐ < 0) stays empty: 0 of 192 γ values. **[CORRECTED 2026-08-19/20 — the earlier form of this clause read “0 of 192 γ values *at every Brans-Dicke ω tested*,” which was an **under**-refutation in two independent ways, and the corrected statement is strictly stronger.]** *(i) The tested grid was excluded before it was scanned.* Brans-Dicke gives `γ_PPN − 1 = −1/(2+ω)`; Cassini (Bertotti, Iess & Tortora 2003) measures `(2.1 ± 2.3)×10⁻⁵`, whose 2σ lower edge `−2.5×10⁻⁵` forces **ω ≥ 4.0×10⁴** for an unscreened massless scalar — so the published grid `ω ∈ {0, 1, 5, 50}` topped out **800× inside the excluded region**, and no solar-system bound appeared anywhere in this sector although the same spacecraft is cited for TEST-25 **[ID NOTE 2026-08-20 — `TEST-25` here resolves to the **site** surface’s *Cassini squeeze*, NOT to this repository’s `Research/EXPERIMENTAL_TEST_CATALOG.md` §TEST-25, which is *a₀ Redshift Evolution*. The two public surfaces collided on this ID when the site cleared an *intra-site* collision by renumbering `TEST-11 → TEST-25` on 2026-08-10. Following this citation into the catalog shipped alongside this document lands on the wrong test; no renumbering is performed here, as which side moves is a cross-repo decision.]**. *(ii) ω was never a free dimension of the scan.* The pinning `C_eff(x₀) = Ω_m` is maintained by sliding `x₀` (0.95 → 3.6×10⁴), not by ω doing physical work: ω enters the background only through `B = 1 − 3ε − 1.5ωε²`, and the pinned point is that factor's positive root `ε_crit(ω) = (−3 + √(9+6ω))/(3ω)`. **Verified here symbolically** (`simulations/publisher_20260820_completion_b_omega_absorption.py`, 11/11): at that root `1.5·ω·ε_crit² = 1 − 3ε_crit` *exactly*, so **the closure pins the product `ω·ε_crit²`, which → 2/3, rather than ω itself** — `ε_crit ~ √(2/(3ω))` and the construction has one physical point, not a family. “192 × 4 ω-values” is **192 × 1**; the scan advertised a dimension the construction does not have. (This closure identity is the reason behind the source's numerical observation that trajectories at ω = 4×10⁴ and 10⁶ agree to <2% — the departures are 1.22% and 0.24%.) *(iii) At the allowed ω the no-go hardens rather than loosening:* `w₀` moves −1.58 (ω = 0) → **−3.18** (ω = 4×10⁴) at γ = 0.489, deeper into the wrong quadrant, with ρ_DE at 2% of its present value by z = 2. *(iv) The screening escape is self-inconsistent rather than a rescue:* giving the scalar a mass via `V(C)` to evade PPN would contribute to `B(x)`, which is the **massless** BD energy density and omits V — the published scan cannot be simultaneously PPN-safe and self-consistent. **The one remaining unexecuted covariant branch** is therefore a density-pinned **massive** scalar (add `V(C)` to `B(x)` and integrate): it evades Cassini by construction, and a potential is exactly the ingredient that could supply the interior maximum of ρ_DE(x) the DESI quadrant requires. Routed, not driven. Refutation count unchanged (6); the completion-B exclusion stands and is stronger. *(Source: `synchronism-site/explorer/findings/completion-b-at-cassini-allowed-omega-the-no-go-hardens-and-omega-was-never-free.md`, 2026-08-19, routed via `Research/proposals/completion_b_cassini_omega_absorbed_by_closure_20260819.md`; the w₀ = −3.18 trajectory value is source-executed and not re-integrated here.)* Continuing the original clause: , and forcing w₀ to DESI's value forces wₐ = +0.23…+0.60, the wrong sign in all four combinations. The class criterion is one line: for any dark energy slaved to matter density, ρ_DE = ρ_m·F(x), continuity gives w_DE = d ln F/d ln x, so DESI's phantom→quintessence crossing requires an interior *maximum* of ρ_DE(x) — and no completion of C = tanh(γ ln(1+x)) produces one (substituted form: monotone; completion A: ρ_DE ≡ 0; completion B: interior *minimum*, crossing in the anti-DESI direction). Residual conditionality: completion B's quasi-static pinning presumes an enforcing sector of negligible stress. The γ = 1/2 branch of the *substituted* model is exactly ΛCDM (inheriting ΛCDM's 2.8–4.2σ DESI tension, no more) — but that degeneracy is a property of the substitution, not of the covariant completions. **[STRENGTHENED 2026-08-19 — the qualifier "background" is removed because it was load-bearing in the wrong direction, and the identities below were re-derived here in SymPy, exactly, for all x and γ]:** the γ = 1/2 branch is ΛCDM **non-perturbatively**, not merely degenerate with it in a background fit. `ρ_DE/ρ_crit = 2x/((1+x)^{2γ} − 1)` equals **exactly 2 at γ = 1/2, with dρ_DE/dρ_m ≡ 0** — ρ_DE is constant in *space* as well as in time, the same inside a void, a cluster and a disk — and the induced contrast `δ_DE/δ_m = 1 + w_DE = ε·(ln(1+x)/x − 1) + O(ε²)`, `ε ≡ 2γ − 1`, vanishes identically at every scale and every order. Same background, same linear perturbations, same nonlinear field configuration. This **retires the standing "background-only, no perturbation sector" caveat**, which had twice been read as an unexplored channel that might convert the dp-gated DR3 registration into a near-term test: the *locality* fork — evaluating C at the local rather than the mean density, which the framework's one-equation postulate actually requires — was executed 2026-08-18 and adds **no independent channel**, because both the expansion-history and clustering terms are linear in the same single ε (no order-ε⁰ term exists). What locality buys is exactly a factor 2 in the observable: required precision σ(fσ₈) relaxes 0.074% → 0.15%. It buys **nothing** in σ_γ, which is `(1/2 − γ̂)/3` in both horns by construction — the same γ-space separation this paragraph already prices at σ_γ ≈ 0.004 (the source finding's "7% improvement" is a central-value artifact, 0.004 at γ̂ = 0.487 vs 0.0037 at γ̂ = 0.489; corrected here). Against DESI DR2's ~1–3% per bin the channel is ~10× underpowered, and the parameter it would measure sits at **ε = −0.022 ± 0.220 = 0.10σ from exact ΛCDM**. Bucket 0 unchanged (0); refutation count HELD at 6 — this closes a fork and refutes nothing. Sources: `synchronism-site/explorer/findings/de-locality-fork-executed-perturbation-channel-buys-exactly-a-factor-of-two.md` (2026-08-18), routed via `Research/proposals/de_locality_fork_perturbation_channel_factor_two_20260818.md` and triaged in `explorations/2026-08-18-de-locality-fork-perturbations-buy-x2-and-session107-is-a-different-refuted-mechanism.md`. A direct likelihood fit of the w(z; γ) family to DESI DR2 + Planck-compressed CMB + DES-SN5YR (site explorer, 2026-08-12; the Δχ² numbers are site-executed and not re-run in-repo, though the pipeline self-validates against DESI's published ΛCDM and w₀wₐ results) confirms this structure at data level: the substituted branch sits at its Λ-corner — best γ = 0.487 ± 0.02, Δχ² = −0.3 vs ΛCDM, +11 behind w₀wₐCDM — and both covariant completions fail the fit outright (A: χ² ≈ 9,900, the pre-1998 background; B: Δχ² ≥ +79 at every ω on that same PPN-excluded grid — and per the ω correction above, one physical point rather than four), so the sector survives only in its non-covariant form and only by being ΛCDM. An earlier site-side "3.4–6.3σ forced-wₐ" exclusion of the substituted branch did not survive this execution — the likelihood never forces the family off the Λ-corner; the verdict is non-novelty, not exclusion. The apparent cross-sector agreement γ_cosmo = 0.487 vs γ_galaxy(SPARC) = 0.489 has **no power to fail** (γ = 1/2 is exactly Λ; 0.489 is exactly MOND simple-μ; the two standard models sit 0.011 apart in γ-space by construction) and is therefore not a concordance. **[SCOPED 2026-08-15]** That deflation was *structural* when written — the SPARC-side uncertainty had never been derived, so the pricing was incomplete. It has now been derived, and it makes the closure **permanent** and the sign of the comparison **convention-dependent**. σ(γ) = **0.11** (stat) — the galaxy-limited value, not the naive point-level 0.029, because rotation-curve points within a galaxy are coherently correlated so the effective N is the number of *galaxies* (jackknife 0.103 / bootstrap(400) 0.113 / √(N/N_gal) 0.121; corroborated in-repo against `simulations/sparc_real_data/MassModels_Lelli2016c.mrt` — 175 galaxies, 3391 points, 0.029 × √19.4 = 0.128). Against the continuous free minimum γ̂ = 0.4892 the comparison is |0.487 − 0.4892| = 0.002, i.e. **0.02σ**; separating 0.489 from 1/2 would need σ_γ ≈ 0.004, which requires ~130,000 SPARC-quality galaxies *and* a global ϒ known to ±0.0007. **The ϒ term does not average down, so γ = 1/2-vs-fit discrimination on rotation curves is closed a priori** — no SPARC-class dataset rescues it. Sharper still, the agreement is an artifact of the ϒ_disk = 0.5 convention: refitting at ϒ_disk ∈ {0.4, 0.5, 0.55, 0.6} moves γ̂ across **[0.27, 0.96]** at flat rms (0.1433–0.1458 dex), and at the *better-fitting* ϒ = 0.55 (Δχ²_naive = −14.7) γ̂ = 0.68 — the "agreement" becomes a ~1.7σ **tension** with lower χ². Single galaxies (UGC11914, UGC03580, NGC5985) each move γ̂ by ±0.03–0.04, ~3× the entire |0.489 − 1/2| offset: the three decimals of "0.489" are noise digits. **Consequence for the parenthetical above**, stated because it qualifies this section's own supporting clause: "0.489 is exactly MOND simple-μ" is an ϒ = 0.5 *slice* statement — q = 2γ spans [0.5, 1.9] across the band — so every "γ ≈ 0.5 from galaxies" statement carries factor-2 slack. γ, a₀ and ϒ are one flat degeneracy (the "profiled a₀ = 5.33×10⁻¹¹ vs derived cH₀/2π = 1.04×10⁻¹⁰, factor 1.96" tension dissolves at ϒ = 0.6, where the profiled a₀ → 1.043×10⁻¹⁰ ≡ the derived value at indistinguishable rms), so **the galaxy sector's shape parameter is unidentified at factor 2** — a sharper statement than "it reduces to MOND," since an unidentified parameter cannot carry a novel prediction. This does not disturb the door-#1 verdict, whose load-bearing kills are the ≤25% admixture bound and the g_bar→ρ locality no-go rather than the exact γ value. **Provenance and what was not re-run**: the stat term is corroborated in-repo as above; the ϒ sweep, the ϒ = 0.55 flip and the a₀ dissolution are **site-executed and not reproduced here** (`synchronism-site/explorer/findings/scripts/sparc_gamma_interval_frozen_likelihood.py` exceeded the in-session compute budget and was terminated at exit 143 before output). Open thread, carried as the source states it: a hierarchical refit with per-galaxy ϒ priors (Li et al. 2018 RAR methodology) is the defensible next rung before any external citation of the interval. Refutation count **HELD at 6**, Bucket 0 = 0 — this prices a comparison and refutes nothing (source: `Research/proposals/sparc_gamma_interval_upsilon_degeneracy_20260814.md`, routed from the site explorer back-annotation `c2b78336`, 2026-08-14). A prospective registration of the no-go (DESI DR3 timeline, ~2027–2028) is filed and gated on dp; its terms as currently proposed [corrected 2026-08-14; the 2026-08-12 form of this sentence gave "the discriminant is the sign of wₐ" — a statistic the direct fit showed is ΛCDM-degenerate on the only surviving branch]: class-level (the interior-maximum criterion above), statistic moved to Δχ²(substituted best fit vs w₀wₐCDM) on DR3 likelihoods, presented as a consistency check rather than a discriminator, with a pre-committed projection-robustness check against the CPL-artifact literature. Sign-error statements verified in `simulations/publisher_20260811_w_eff_erratum_check.py`; completion algebra, the EdS result, a_end, and the completion-B table independently reproduced in `simulations/publisher_20260812_covariant_00_checks.py` (source: `synchronism-site/explorer/findings/covariant-00-component-sign-lock-dies-desi-nogo-hardens.md`, 2026-08-11, routed via `Research/proposals/covariant_00_component_sign_lock_dies_desi_nogo_hardens_20260811.md`; direct-fit source: `synchronism-site/explorer/findings/gamma-family-direct-fit-desi-dr2-substituted-is-lcdm-covariant-excluded.md`, 2026-08-12, routed via `Research/proposals/gamma_family_direct_fit_desi_dr2_20260812.md`).

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
