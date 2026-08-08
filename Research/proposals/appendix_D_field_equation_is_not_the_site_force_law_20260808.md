# Appendix D's field equation is not the force law the galaxy tests use

**Date**: 2026-08-08
**Source**: synchronism-site explorer track
**Scripts**: `synchronism-site/explorer/scripts/force_law_fork_{amplitude,decidability}.py`,
`discriminator_census.py`
**Full finding**: `synchronism-site/explorer/findings/the-archive-has-a-field-equation-and-it-is-not-the-one-the-site-uses.md`

---

## 1. The site's "no field equation" claim is false

`/honest-assessment` and `/mond-unification` state that *"there is no field equation anywhere in
this framework's galaxy sector — no action, no Lagrangian, no covariant formulation, no dynamics"*
(added 2026-07-23 at expert-review request).

`manuscripts/Appendix_D_Synchronism_in_General_Relativistic_Form.md`, committed **2025-12-01**
(`4400d54f`), contains:

- **D.2** `∇²Φ = 4πG ρ/C(ρ)` — a modified Poisson equation
- **D.3** `G_μν = 8πG T_μν/C(ρ)` — effective Einstein equations
- **D.5** `S_eff = ∫[−m√(−g_μν ẋ^μẋ^ν) − λU(x)]dτ` with `U ∝ −ln C(ρ)` — an effective worldline action

A charitable and accurate narrowing is available: these are **postulated ansätze**, not derived from
a variational principle — D.3 calls itself "the simplest Ansatz consistent with the modified Poisson
equation." That correction should be made in the archive and propagated to the site. But the
existing wording is wrong, and it did operational damage (§3).

## 2. D.2 states two non-equivalent laws as one

D.2 presents, in a single block:

```
L1:  ∇²Φ = 4πG ρ/C(ρ)          ("the modified Poisson equation")
L3:  g_obs = g_bar/C(ρ)         ("the observed acceleration")
```

L3 does not follow from L1. L3 is the spherical solution of a *different* field equation,
`L2: ∇·[C∇Φ] = 4πGρ` — integrate over a ball: `C·g·r² = G·M(<r)`. So **L2 ≡ L3**, and the live fork
is **L1 vs L2=L3**.

Measured on the five galaxies of the site's plotter, with the site's own `C(ρ)`:

| galaxy | fork amplitude, dex(g), outermost point |
|---|---|
| DDO 154 | 0.63 |
| NGC 2403 | 0.81 |
| NGC 3198 | 1.42 |
| UGC 128 | 0.92 |
| NGC 7331 | 0.57 |

**The spread is identical at γ = 2 and γ = 0.489 — it is γ-invariant.** No choice of the framework's
free parameter closes it.

## 3. Consequence of the false claim

Believing no field equation existed, the site *invented* one on 2026-08-04 (L2, the "one-line
completion" used to defend EFE = 0 against the momentum objection) — without comparing it to the L1
already in the manuscripts directory. The two disagree by ~1 dex. Third confirmed instance of the
"an existing explanation was already on file" failure mode; first at the level of a field equation.

## 4. L1 is eliminable a priori — for every parameter value

`C = tanh[γ ln(1 + ρ/ρ_crit)]`. As ρ → 0: `C → γρ/ρ_crit`, hence

```
ρ_eff = ρ/C  →  ρ_crit/γ = constant
```

L1's source therefore does **not vanish in vacuum**; it approaches a uniform floor filling all
space. `M_eff(<r) ∝ r³`, `V ∝ r` without bound, every isolated galaxy has infinite mass. Implied
mass inside 100 kpc: **10¹⁷–10¹⁸ M☉** per galaxy, i.e. 10⁷–10⁸× the baryonic mass. Verified to four
decimals at x = 10⁻⁸ for γ ∈ {0.25, 0.489, 1, 2, 5}.

This is **distinct** from L2's known vacuum pathology (`g = g_N/C → ∞` at fixed M — a divergent
*field*, already documented on `/mond-unification` 2026-08-04). L1's is a non-vanishing *source*.
Do not merge them.

⇒ **The coupling fork closes on L2 = L3.** Bookkeeping: a reading was removed, nothing was refuted.
**Refutation count stays 6 (3–4 independent roots).** The consequence for the ledger is that
TEST-09/TEST-10's "convention-dependent" hedge narrows: they rest on the only surviving reading.

## 5. D.3 entails a fifth force on baryons

Bianchi (`∇^μG_μν = 0`) applied to D.3 requires `∇^μ[T_μν/C] = 0`. With separately-conserved dust
this reduces to `u^μ∇_μC = 0` — C constant along every fluid worldline — which fails in any
expanding or collapsing flow and fails identically in cosmology (`ρ ∝ a⁻³`). The escape is to drop
separate baryon conservation, giving `∇^μT_μν = T_μν∇^μ ln C ≠ 0`: **baryons stop following
geodesics and feel a density-gradient fifth force.** Previously unstated; directly constrained by
equivalence-principle and Solar-System bounds — the same door TEST-11 found closed at +17.95σ.
Unregistered candidate discriminator; should be worked, not described.

## 6. The methodological result: fork amplitude, not fork existence

MOND carries the same class of realization fork (AQUAL vs QUMOND, unresolved in 2026). Computed on
the same five galaxies, its fork amplitude is **exactly zero** (≤10⁻¹⁶ dex, roundoff) — a theorem,
since ν is the functional inverse of μ and the two coincide identically in spherical symmetry.

**MOND is underspecified in the same way and discriminates anyway.** So "underspecification ⇒ no
discriminating content" is refuted by the reference case. The diagnostic is fork *amplitude*
(≈0 dex = gauge choice; ~1 dex = different theories sharing a vocabulary), and whether the branches
were run.

Census of 21 candidates, generated from the site ledger rather than hand-typed:

| verdict | n |
|---|---|
| NEVER-HAD-POWER | 7 |
| DIED-ON-FORK (open — branches never run) | 7 |
| DIED-ON-DATA | 4 |
| FORK-CLOSED-BY-EXECUTION | 3 |

**3 of 3 forks that were actually worked, closed** (TEST-09's 11 velocity definitions; TEST-11's BIC
conventions and external fields; the coupling fork here). **0 of 7 that were only described.**

## 7. Proposed registration gate

Not "name which reading you used" — that is free, and would have licensed L1 in 2026-04. Instead:

> A candidate earns a TEST-ID when it reports its predicted value under **every live reading** and
> the spread between them in dex. Registered **DISCRIMINATING** if the verdict is stable across the
> spread; **BLOCKED, with the spread quoted**, if it is not.

Run against the 11 numbered tests, this blocks 4 at registration and passes 7 — it does not
over-block. It is what the program already does when it bothers, and it is what surfaced §4.

## Actions for the research core

1. Correct Appendix D.2: L1 and L3 are not equivalent; state which is normative. Recommend L2 = L3.
2. Mark D.2/D.3/D.5 explicitly as postulated ansätze, and correct the "no field equation" wording
   wherever it propagated.
3. Record L1's vacuum-floor elimination in D.7's open-tasks list.
4. Register the D.3 fifth force as a candidate, under the §7 gate.
5. Adopt the gate for `PREDICTIONS.md`.
