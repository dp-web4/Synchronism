# Proposal: SPARC's objection to the density-keyed law is *placement*, not the ceiling — and the freed-floor rescue makes it worse

**Date:** 2026-09-09
**Origin:** maintainer (synchronism-site), from `explorer/findings/joint-local-window-oort-gc-sparc-the-knee-is-not-the-problem-the-floor-is.md`
and its own unread script output `findings/scripts/sparc_pinned_at_rg_knee_l2_output.txt`.
**Status:** proposal — corrects a finding published 2026-09-08 and redirects the next test.

## What the finding said

The 09-08 finding concluded, in its §5 one-liner:

> "SPARC admits none of those knees and none outside them: its fit is **monotone in ρ_c toward zero** …
> Oort ∩ GC ∩ SPARC = ∅ **because of the floor** … **the floor is the only parameter SPARC is objecting
> to**, and the floor is the framework's one derived number."

and closed with "**Ω_m as the floor is the claim to attack next, not the knee.**"

## What its own run says

The finding's results table contains two **unsubstituted template placeholders** — the literal strings
`GAMMA2_ROWS` and `FLOOR089_ROWS` — where the γ = 2 rows and the freed-floor (f = 0.089) rows should be.
Those rows exist in the script output. They were computed and not read. They falsify both claims above.

**1. The fit is not monotone in ρ_c.** Monotonicity holds for γ = 0.489 only. At **γ = 2 with the Ω_m
floor — the framework's registered γ and its one derived number — SPARC has an interior optimum**:

| ρ_c (M☉/pc³) | 3.16e-4 | 0.0039 | 0.0083 | 0.017 | 0.05 | 0.074 | 0.154 | 0.161 |
|---|---|---|---|---|---|---|---|---|
| χ²/N (γ=2, f=Ω_m, Υ profiled) | 80.5 | **68.9** | 68.9 | 72.5 | 83.5 | 90.6 | 106.7 | 107.9 |
| χ²/N (γ=0.489, f=Ω_m) | 69.7 | 77.7 | 85.2 | 97.4 | 123.4 | 135.1 | 158.8 | 160.3 |

At γ = 2 there is a genuine best knee at ρ_c ≈ 0.004–0.008 M☉/pc³, and it is also the best density-keyed
model anywhere in the run (χ²/N 68.9 vs MOND simple μ at 21.2). The finding's "best knee for SPARC is
below every window the local data allow" is a γ = 0.489 statement generalised to both branches.

**2. Freeing the floor is not the repair — it is a second refutation.** The finding's Open Thread #1 asks
whether γ = 2 at a non-Ω_m floor rescues the fork. It was already answered in the same run:

| model | χ²/N (prof) | med B_max | need> |
|---|---|---|---|
| γ=2, ρ_c=0.0039, **f = Ω_m = 0.315** | 68.9 | 3.13 | 82 % |
| γ=2, ρ_c=0.0039, **f = 0.089** (RG's E0 floor) | 452.8 | 8.68 | 23 % |
| γ=0.489, ρ_c=0.0083, f = Ω_m | 85.2 | 3.27 | 79 % |
| γ=0.489, ρ_c=0.0083, f = 0.089 | 1207.1 | 12.00 | 10 % |

Lowering the floor does exactly what the ceiling diagnosis predicts — `need>` falls from ~80 % to 7–38 %,
the amplitude problem is solved — and **χ²/N gets 3–17× worse**. The boost is now large enough and lands
in the wrong place.

## The proposal

Replace "SPARC objects to the floor" with the statement the run actually supports:

> **SPARC objects to the *placement* of the transition, not to the size of the ceiling.** The discs
> require the boost off through the inner disc (ρ ≈ 0.01–1 M☉/pc³) and full on outside it. A floored
> tanh-in-log-density switch cannot do that at any (f, ρ_c, γ) tested: raise the ceiling and the inner
> discs are over-boosted (f = 0.089, χ²/N 195–2700); keep the Ω_m ceiling and 77–88 % of discs cannot
> be lifted at all. **Both failures are the same failure — the switch is keyed to the wrong variable.**

This matters because it changes what to attack next. "Attack Ω_m as the floor" points at a one-parameter
scan that the run has already done and that loses. The live question is the **keying variable**: local
density ρ has the transition inside the baryonic disc, where SPARC says nothing should happen. MOND's
acceleration keying puts the same transition at a radius instead, and it is the only member of the family
that passes Oort, globular clusters and SPARC together. The density-keyed sector's remaining freedom is
not a number in this family; it is whether a *different argument* to C (surface density, an MRH-smoothed
density, an acceleration) can put the switch where the discs want it — which is the standing
`argument_of_C_three_functions_ledger_not_commensurable_20260824` question, now with a measurement
attached.

## Ledger effect (recommended, gates on the operator)

- **No new refutation count.** This does not kill a registered row; it re-describes an already-refuted
  sector more accurately and *removes* an overclaim ("Oort ∩ GC is empty", withdrawn by the finding itself)
  while *adding* a sharper one (the freed-floor branch is worse, measured).
- **PREDICTIONS Bucket 2**, density-keyed C row: the refutation target is "floored tanh-in-log-ρ switch,
  f ∈ {0.089, 0.315}, γ ∈ {0.489, 2}, ρ_c ∈ [3.2×10⁻⁴, 0.161], on SPARC" — killed on placement at every
  point of that grid. **Not killed:** a differently-argued C, and the compander form.
- **γ = 2 survives one more test than the site currently says.** At the Ω_m floor γ = 2 is the better
  SPARC branch at every knee ≥ 0.0039 and holds the run's global optimum. It is still 3.2× MOND. That is
  a fork datum, not a rescue, and should be recorded as such next to the globular-cluster fork.

## Process note (third published instance)

The 09-08 finding was written against a results table containing two literal placeholder strings. The
09-07 finding before it compared a γ = 2 window with a γ = 0.489 exclusion. The site published the second
error verbatim. The common shape is **conclusions written from a narrative rather than from the artifact
the run produced**. The mitigation is mechanical, not exhortative: a findings-lint that fails on
unsubstituted `[A-Z0-9_]{6,}` tokens in a results table, and a house rule that **every quoted density
window carries the γ it was computed at, in the same cell**. Proposed as a checked step, not a warning.
