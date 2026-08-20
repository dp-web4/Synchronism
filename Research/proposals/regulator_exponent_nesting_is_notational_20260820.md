# The nesting of the galaxy sector inside MOND is imposed by a notational convention, not by the model

**Back-annotation from synchronism-site explorer, 2026-08-20.**
Source: `synchronism-site/explorer/findings/regulator-exponent-the-nesting-in-mond-is-a-notational-convention.md`
Scripts: `synchronism-site/explorer/findings/scripts/regulator_exponent_n_real_sparc.py`,
`regulator_exponent_outer_slope.py`, `regulator_exponent_n_crossval.py`.

**Refutation count UNCHANGED at 6.** Nothing newly refuted, nothing newly confirmed. What
changes is the *status* of an a-priori argument, plus one factual correction, plus one
registrable executed negative result.

## 1. The claim

`C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))` contains an unnamed index. Writing the density ratio to a
general power,

    C_p(x) = tanh(γ·ln(1 + x^p)),      the archive's equation is p ≡ 1,

the deep limit is `C → γ·x^p` **exactly** (verified to 8 significant figures at x = 10⁻¹⁰,
independent of γ). In the reading where `x = g/a₀` — the reading under which the archive's own
identity `C ≡ μ` holds (`γ=1/2, p=1 ⇒ C = x/(x+2) = μ_simple(x/2)`, exact to 1.1×10⁻¹⁶) —

    g = (g_bar·a₀^p/γ)^{1/(p+1)},   V² ∝ r^{(p−1)/(p+1)},   M ∝ V^{2(p+1)},   s_deep = 1/(p+1).

`p = 1` ⇔ asymptotically flat rotation curves ⇔ BTFR slope 4 ⇔ MOND's deep limit ⇔, by
Milgrom (2009, ApJ 698, 1630), **spacetime scale invariance**. MOND *derives* the index from a
symmetry. Synchronism *inherits* it from a regulator chosen for finiteness at ρ = 0.

## 2. Why this matters to the archive's own conclusions

The site's a-priori argument (`/for-researchers`) runs: the bounded boost `B ≤ 1/Ω_m` is *"the
framework's only structural difference from MOND"*; a ceiling is a restriction; therefore the
galaxy sector is `MOND ∩ {B ≤ 3.17}`, a strict submodel, which "cannot win."

The premise is false. There is a second structural difference and it has the **opposite sign**:
a free `p` is an *extension*, so the family **contains** MOND at `p = 1`. A supermodel can be
selected over its nested special case. The tie-or-lose conclusion is conditional on `p ≡ 1`.

Two further archive results are conditional on the same index:
- `d²C/dx² ≈ γp(p−1)x^{p−2}` ⇒ concavity-everywhere in ρ (and hence "the S-curve is a log-axis
  artifact") holds **only for p ≤ 1**; a real inflection appears at p = 1.5 (x = 0.545) and
  p = 2 (x = 0.819).
- The TEST-09 statement that the bounded-boost deep limit is BTFR slope → 2, "the opposite end
  of the ladder from MOND's 4", is a statement about `p ≡ 1` plus the asserted `B_max`.

## 3. Executed on real SPARC (Lelli+2016; 2700 pts, 149 galaxies; Q<3, i>30°, e_V/V<0.10)

**Power gate first**, because Desmond, Bartlett & Ferreira 2023 (MNRAS 521, 1817) showed that
with a *free* functional form SPARC cannot recover the deep slope even from MOND mocks.
Injection–recovery inside the rigid 3-parameter family: `p` recovered unbiased (worst bias
0.038) with scatter ≈ 0.05 at every injected `p ∈ {0.70, 0.85, 1.00, 1.15, 1.30}`,
`dp̂/dp_true = 1.039`. **The gate passes — but the power comes from the family's rigidity, not
from dynamic range** (only 0.4 % of points lie below a₀/100), and must be labelled as a prior.

**Result.** In-sample `p̂ = 0.762`, `Δ(2 lnL) = 15.0` vs `p ≡ 1`. But the two best-fit models
differ by **at most 0.038 dex** across the SPARC-covered range against `σ_int = 0.122 dex` —
the significance is 2700 correlated points each carrying ≲0.3σ, the same structure as the
ΔBIC = +184 the archive corrected to ≈+7 at N_eff = 175.

**Adjudication out-of-sample** (10-fold CV over galaxies, 3 shuffles, 149 galaxies):

    p free vs p ≡ 1 :  Δ lnL/pt = +0.00194 ± 0.00592  =  0.33σ
                       84/149 galaxies better, binomial p = 0.140, Wilcoxon p = 0.234

**Not selected. `p = 1` stands.** Independently corroborated: arXiv:2608.08945 (Aug 2026)
establishes a 0.106 dex per-galaxy nuisance floor for one-parameter structural corrections to
the RAR in SPARC; our model separation is 0.038 dex, i.e. under the floor.

**Three robustness checks, all clean.** (i) `p̂` moves only 0.091 = 0.62 bootstrap σ across
three qualitatively different outer sigmoids (tanh / `u/√(1+u²)` / `1−e^{−u}`) while γ̂ moves
×36 and â₀ ×213 — the argument exponent is separately identifiable from the outer-sigmoid
freedom that the 2026-05-02 sweep found unconstrained. (ii) Neither end constrains `p` alone:
deep-only (`g_bar < a₀/3`, N=1568) gives 1.186 ± 0.697, high-only (N=627) gives 0.885 ± 0.851;
only the join reaches 0.146. **The precision comes from the shape the family imposes between
the ends, not from either end's data** — label it as a prior. (iii) No evidence against a
universal `p`: per-galaxy fits (125 galaxies) give IQR [0.591, 0.916] with correlations
−0.064 / +0.103 against acceleration regime / mass proxy. Permutation null piles at the
parameter bound (median 8.0), confirming the estimator responds to real structure — but because
the null is at a boundary its "6.6σ" is **not** a significance and must not be quoted as one.

**This null is a real null, not an unrunnable one** — the archive's recurring failure mode
(TEST-02 hung both ways, TEST-04a on the wrong observable, the γ concordance with no power to
fail, C doubly unanchored). Here the estimator was shown to have power *before* the
measurement, and returned a null.

## 4. The by-product that bears on the archive's parameter accounting

Same likelihood, same data, nested comparison:

| freeing | in-sample Δ(2 lnL) | out-of-sample Δ lnL/pt |
|---|---|---|
| **γ** (MOND simple-μ → the archive's equation) | **0.39** | **−0.00159 ± 0.00119 (−1.34σ)** |
| **p** (archive's equation → the extension) | 15.01 | +0.00194 ± 0.00592 (+0.33σ) |

Freezing γ at MOND's exact ½ **predicts held-out galaxies better** than fitting it. The
archive's headline fitted parameter is not merely powerless — out of sample it is mildly
harmful, which is the signature of overfitting a degenerate direction. This is independent
empirical support for the "Reparametrization-until-proven" prior and for the standing result
that γ = 1/2 is the exact MOND point.

Pipeline validation: at `p ≡ 1` this independent implementation returns **γ = 0.4834**,
reproducing the archive's γ_SPARC = 0.489 to 1 %.

## 5. Systematics, and a registered prediction that came back partly wrong

Registered before running: *"ϒ rescaling translates the RAR horizontally and cannot change a
slope, so p̂ should be ϒ-invariant where γ̂ was destroyed."* **Partly wrong** — ϒ scales only the
stellar component, so it changes the star/gas mix and hence the *shape* of g_bar, and it
re-selects which radial windows pass the cuts. Measured spans across ϒ_disk ∈ [0.30, 0.80]:
γ ≈ 140 % of its mean (the 2026-08-14 result), **p only 23 %**. The defensible form is
comparative: **p is ≈ 6× more ϒ-stable than γ**, not ϒ-invariant.

The independent kinematic estimator (`s_V = d lnV/d lnr`, no ϒ in the statistic, windows
required deep *and* mass-converged; 42 of 132 galaxies survive) gives `p` in [0.60, 0.91]
across estimators and [0.599, 0.889] across 15 window definitions — **direction robust
(`p < 1` in all five estimators), significance not** (0.8σ to 5.1σ). The naive
inverse-variance mean is driven by three galaxies carrying 57 % of the weight; the median does
**not** exclude `p = 1`. Quote the trimmed mean or the median.

## 6. Scope limit that must travel with this result

The `p` direction is alive **only** in the `g`-variable reading. In the archive's literal
variable `x = ρ/ρ_crit`, `C_p(ρ)` is a local algebraic function of ρ, and the 2026-08-02
form-free bound (`σ(log B | ρ) = 0.1613 dex` vs `σ(log B | g_bar) = 0.1178`) is an infimum
over *all* functions of ρ — no `p` escapes it. So this is a one-parameter generalization **of**
MOND's deep limit, not an escape from MOND's variable. It remains a direction MOND cannot
take, because there `p = 1` is a theorem.

## 7. Recommended archive actions

1. Correct wherever the `+1` is described as *excluding* power-law behavior as ρ → 0. It
   **creates** the power law and pins its index to 1.
2. Replace *"the framework's only structural difference from MOND"* with *"the only structural
   **restriction**"*, and record the extension direction alongside it.
3. Register the executed negative result with a TEST-ID and an out-of-sample registered
   statistic (per-galaxy Δ lnL under galaxy-level K-fold), including the pre-stated power
   σ(p̂) ≈ 0.05. Without an ID it will be invisible to ledger audits.
4. Credit line required: the theoretical content — that scale invariance fixes the deep
   exponent — is Milgrom (2009). This deliverable is the *measurement* and the *null*.
5. Do not stack the four `p < 1` directions (RAR family fit 0.762; kinematic outer slope
   ≈0.76; the archive's own observed BTFR 3.75 ± 0.10 ⇒ p = 0.875 ± 0.050; arXiv:2006.06700's
   independent-data "weak preference for a steeper low-acceleration slope"). Three are
   SPARC-correlated; only the last is independent, and it is explicitly weak. The out-of-sample
   adjudication is the one that counts, and it is null.
