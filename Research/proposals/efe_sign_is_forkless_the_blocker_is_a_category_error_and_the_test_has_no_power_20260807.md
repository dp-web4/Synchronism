# EFE = 0 is fork-free; the registration blocker is a category error; the test has no power

**Date**: 2026-08-07
**Origin**: synchronism-site explorer session, executing
`explorer/topics/efe-sign-convention-dependence.md` (P0, maintainer, same day).
**Supersedes** the operative reading of
`force_law_fork_blocks_efe_registration_and_makes_count_convention_dependent_20260807.md`
on the EFE axis only. That proposal's *count* argument (TEST-09/10 are one root) is untouched.
**Full derivation and execution**: `synchronism-site/explorer/findings/efe-sign-is-not-convention-dependent-the-blocker-is-a-category-error.md`

---

## 1. The sign does not fork

Embed a system in a uniform external field `g_ext`, baryons fixed. A uniform field exerts no
tide, so ρ is unchanged by construction. Then:

| reading | law | response to `g_ext` |
|---|---|---|
| amplitude | `v² = v_b² + (V_flat·C)²` | **0** |
| division | `g_obs = g_bar/C(ρ)` | **0** |
| multiplication | `g_obs = C(ρ)·g_bar` | **0** |

Executed at `e = 0, 0.05, 0.5, 5.0`: identical to every printed digit in all three columns.
**Zero has no sign. The force-law fork does not reach the EFE.**

## 2. The fork that does decide it is the ARGUMENT fork, and it is untracked

| C's argument | class | field equation | EFE |
|---|---|---|---|
| `ρ` | source | `∇·[C(ρ)∇Φ] = 4πGρ`, **linear in Φ** | 0 exactly |
| `g_bar` | source (Poisson image of ρ) | linear in Φ | 0 exactly |
| `\|∇Φ\|` | **solution** | nonlinear — this *is* AQUAL with ν = 1/C | MOND's EFE |

The program has been tracking *how C multiplies* and has never separately tracked *what C
eats*. The EFE is a pure function of the second axis and a null function of the first. No
branch of the 3×3 space yields a nonzero EFE of framework-specific sign — you get 0, or you
get MOND's.

**Archive correction.** `Session215_EFE_Predictions.md` is the original EFE = 0 derivation and
it runs on an acceleration-keyed `C(a)`, concluding *"C(a) depends ONLY on local acceleration a
→ No external field effect."* **That is a non-sequitur**: MOND's μ also depends only on the
local acceleration, and that is precisely how the EFE arises, because the acceleration at a
point inside an embedded system includes the external field. The conclusion is true under the
density-keyed C now in use; **the argument on file is invalid and was never redone.** This is
the numbers-outliving-their-computation failure mode applied to a *structural* claim.

## 3. The 2026-08-04 "opposite sign" is a different channel

It reproduces, and it is about `ρ_ext`, not `g_ext`:

- division: ambient density ↑ ⇒ velocity **deficit**
- multiplication / amplitude: ambient density ↑ ⇒ velocity **surplus**

| channel | keyed on | forks by force law? | status |
|---|---|---|---|
| EFE | `g_ext` (external field) | **no** — always 0 | claimed, unregistered |
| ambient-density effect | `ρ_ext` (external mass density) | **yes** | **TEST-05: registered, executed, r² = 0.0001, kill fired** |

The convention fork blocks the ambient-density channel — which was already registered and
already refuted. It cannot also be what blocks the EFE channel. **The stated blocker is a
category error.**

## 4. Executed: the test has no power (N = 141)

Chae+2020 Table 2 (erratum-corrected arXiv v2) × the registered TEST-08 per-galaxy run.
Matched filter `Δ_pred = log₁₀(−s+√(s²+1))`, `s = e_env/(2√(g_bar/a₀))`, regressed against the
measured RAR offset with ambient density controlled. Pre-declared three-branch verdict rule.

Raw: **β_E = +2.118 ± 0.740**, 95% CI [+0.67, +3.57] — CI excludes 0, EFE-directional, 2.9σ.

**That reading is wrong.** `Δ_pred` is 71.5% the galaxy's own mean acceleration (r = +0.846 with
⟨x₀⟩) and only 13.9% environment. 20,000 permutations of `e_env` at fixed ⟨x₀⟩ and density:

```
permuted null (EFE = 0):  β_E = +1.351 ± 0.390       <- NOT zero
observed vs null:         z = +1.97,  two-sided p = 0.051
```

Against the correct nulls:

| hypothesis | β_E | distance | inside 95% CI |
|---|---|---|---|
| EFE = 0 (confound only) | +1.351 | 1.04σ | yes |
| MOND+EFE | +2.351 | 0.31σ | yes |
| **separation** | | **1.35σ** | |

**Verdict: no power.** Signal `sd(Δ_pred) = 0.0140 dex`; noise `σ_resid = 0.1221 dex`; ratio
0.115. Chae's own low-acceleration cut (N = 106) gives 1.05σ — no better. Using all 175 SPARC
galaxies leaves a factor ~2 shortfall.

**The deficit is dynamic range, not N.** SPARC is a field sample: `e_env ∈ [0.011, 0.057]`.
3σ needs ~2× that range, reachable at `e_env ≈ 0.1–0.3` — group/cluster members and satellites
at small host-centric distance. **That is Session 184 Scenario B and Session 215's
satellite-vs-field design, on file since 2025-12, never executed.** The design was right; the
sample was never assembled.

## 5. Even powered, it is refutation-only

`EFE = 0` is the Strong Equivalence Principle. **GR predicts it; ΛCDM predicts it.** The
ledger's stated alternative is the composite *MOND+EFE+ΛCDM*. A confirmation of EFE = 0 is
therefore shared, not selecting:

| outcome | refutes | selects Synchronism? |
|---|---|---|
| EFE detected | Synchronism **and** ΛCDM | no |
| EFE = 0 confirmed | MOND+EFE only | **no** |

TEST-12 can cost a seventh refutation and can never earn a first confirmation. Registering it
would leave "0 executed tests could select Synchronism" exactly where it is. The framing of
EFE = 0 as *"the framework's only structurally discriminating prediction"* is an over-claim on
the selecting side.

## 6. Methodology contribution

> **Declare the null by permutation, not by convention.** When a test statistic is a function of
> more than one input, generate the null by permuting *the input the hypothesis is about* while
> holding the others fixed. "Coefficient = 0" is an assumption about the null.

Prior instances on this ledger were *fabricated* nulls (invented p-values, retro-fitted
thresholds). This one would have been an honestly computed statistic against an honestly stated
but structurally wrong null — harder to catch, and it would have survived a citation walk,
because every number in it is real and reproducible.

## 7. Scope discipline

- **Refutation count stays at 6.** Nothing here is a refutation; §4 is a *prevented* one.
- The 2026-08-04 retraction of "EFE = 0 refuted by Chae+2020" (baseline off 2–4 dex ⇒
  not-evaluable) is untouched and was not re-derived.
- The count argument in the 2026-08-07 force-law proposal (TEST-09/10 = one division-law root)
  is unaffected — that is the `f_DM = 1−C` axis, not the EFE axis.
- **Do not register TEST-12 until an environment-selected sample exists.** Registering an
  unpowered test is what produced TEST-04a.
