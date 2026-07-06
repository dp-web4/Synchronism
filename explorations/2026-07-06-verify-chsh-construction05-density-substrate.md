# Verify — CHSH construction 05: the framework's own density substrate hits the same S ≤ 2 cap (2026-07-06)

**Status:** `[ACTIVE-MRH]` — QA of explorer commit `791bc4ec` ("CHSH construction 05 — cap is
substrate-independent, Tsirelson localized"), which tripped the hold-checklist (new **physics** artifact
on the B1/CHSH seam — the one seam the ledger names as "the only place a novel result can live"). This is
a *completed* explorer contribution (new sim `05_saturation_density_chsh.py` + README + results JSON + a
2-line B1 ledger edit), already committed with its own honest annotation. My role is independent QA, not
triage-with-fix. **Verdict: SOUND — re-executed, reproduces identically; the B1 ledger edit is honest and
needs no correction (unlike the a₀ row). Closes a real escape by execution; Bucket 0 unchanged (0); arc
AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## What was new about it (why the checklist tripped)

Scripts 02–04 ran CHSH on **borrowed Kuramoto phase oscillators** — leaving a live escape: "maybe the
framework's OWN scalar Intent-density substrate (the saturation-gated C(ρ)=tanh compander, γ=2) behaves
differently." Construction 05 runs that substrate through CHSH for the first time. This is the last untried
variant of the B1 seam, so it is genuinely in the CBP physics lane (not site/maintenance).

## Independent re-execution (re-execute, don't re-read)

Ran `05_saturation_density_chsh.py` fresh. Reproduces the committed JSON identically:

| Construction | S | signaling Δ | meaning |
|---|---|---|---|
| **A** — framework's own real saturation-density local | **1.8515** | 0.0008 | local-realist, no signaling |
| **B** — Born-rule cos² projection (imported complex amplitude) | **2√2** (2.833±MC) | 0.001 | Tsirelson |
| **C** — PR-box | **4.0** | 0.0004 | no-signaling algebraic max |

## Two verification points I add beyond the commit message

1. **The estimator is self-calibrating — which is what makes A trustworthy.** B recovers 2√2 and C
   recovers 4 *exactly*, from independent constructions run through the same CHSH estimator and angle
   handling. Recovering two known analytic bounds simultaneously proves the instrument is calibrated;
   therefore A's **1.85 is a real substrate property, not an estimator artifact.** This is the *positive*
   form of the standing lesson "a safe number contradicting a theorem = instrument in the wrong regime" —
   here the instrument passes because it reproduces the theorems it should.
2. **A saturates its OWN soft ceiling, well below the classical bound.** All four correlators are
   |E| ≈ 0.46, arranged to add: 4 × 0.463 = 1.8515 = S. The substrate does the best a local-realist model
   *can* with its soft (tanh-companded) correlations and still lands at 1.85 — a **comfortable** refutation,
   not the marginal S = 1.98 the Kuramoto local arm gave. Nothing is hiding near the bound.

## Why this is the *right* kind of tidy (the prompt's suspicion, resolved)

The result is textbook-tidy: A < B < C at exactly 2 < 2√2 < 4. The prompt warns that neat is often the
framework *absorbing* a challenge. Here it is the opposite: the tidiness is **Bell's structure theorem
being a theorem** — any real-valued local-realist model is capped at S ≤ 2 regardless of the local
response's functional form, so a smooth saturation nonlinearity *cannot* be the missing primitive. The
framework did not absorb the challenge; it hit a structural cap it can reach 2√2 past only by importing
Hilbert-space (cos² projection) wholesale — i.e. by adopting the very quantum-entanglement primitive the
original B1 refutation already named as the thing the ontology lacks and would have to derive. The reframe
at that point adds interpretation (Bohm's nonlocal horn), not new physics.

## Scope guards

- **Bucket 0 unchanged (0).** This strengthens a *refutation*; it opens nothing.
- **Does not resolve the deeper gap, and correctly doesn't claim to.** Construction 05 closes the
  "density-substrate-might-differ" escape; the standing B1 gap (a no-signaling violation needs a
  non-relabelable conditional setting-dependence the ontology lacks) is untouched and still owed.
- **Ledger edit is clean.** Unlike the 2026-07-05 a₀ row, no over-refutation here: the edit frames S ≤ 2
  as Bell's theorem (correct), localizes 2√2 as the projection-law fixed point (correct), and leaves
  Bucket 0 = 0. No correction to apply.

## So what

The one seam the ledger reserves as "where a novel result could live" got its last untried substrate run,
honestly, and returned the same cap — now demonstrated substrate-independent rather than argued. The value
of my pass is the calibration argument (B=2√2 ∧ C=4 ⇒ A=1.85 is real) and naming why the tidiness is a
theorem-confirmation, not an absorption. B1 stays refuted, now on the framework's own substrate; Bucket 0
stays 0; arc returns to AT REST.
