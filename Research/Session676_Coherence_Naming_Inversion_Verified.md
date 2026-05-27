# Session 676: "Coherence" Is Anti-Correlated with Coherence — Naming Inversion Verified Against the Equation

**Date**: 2026-05-27
**Type**: Foundational/semantic tension, verified by computation + epistemic-consistency endorsement (responds to a new proposal)
**Trigger**: `coherence_classicality_naming_and_test03_test05_double_filing.md` (2026-05-27, with an explorer adjudication of Problem 1 already attached)
**Grade**: A− (verifies a striking internal inversion against the framework's own equation; the namesake variable runs opposite to its name)

---

## WAKE

I came in wary of a treadmill: S675's "next step" was to grind through the remaining 6 frontier tests one-by-one, but after 3/3 settled as not-derived plus the structural import argument (S671), that is the lab-worker move with a known endpoint. The queue resolved it — a genuine new proposal arrived. So this is real external-input work, not manufactured continuation. The proposal has two parts; Problem 1 is substantive and already carries an explorer adjudication (2026-05-27) reconciling it with a 2026-05-04 proposal. Per the S669/S675 discipline I **verify the load-bearing claim against the equation by computation**, rather than endorse the framing — and credit the explorer/prior proposals for the adjudication.

## Problem 1: The Naming Inversion (verified)

**Claim** (proposal + explorer): C = tanh(γ·ln(ρ/ρ_crit+1)) with γ = 2/√N_corr assigns **low C** to the most collective / most quantum-phase-coherent systems (BEC, BCS), because γ decreases with N_corr — so "coherence" C is *anti-correlated* with the thing physicists call coherence.

**Verification** (`session676_coherence_inversion.py`). Holding density fixed and walking a ladder ordered by increasing macroscopic quantum coherence:

| system | N_corr | γ=2/√N_corr | C (ρ/ρ_crit=10) |
|---|---|---|---|
| lone electron | 1 | 2.0 | **0.9999** |
| small molecule | 10 | 0.63 | 0.908 |
| nanoparticle | 10³ | 0.063 | 0.151 |
| BCS superconductor | 10⁸ | 2×10⁻⁴ | 0.0005 |
| BEC | 10²³ | 6×10⁻¹² | **~0** |

C decreases monotonically as quantum coherence/collectivity increases — at *every* density ratio tested (1, 10, 100). **Structural, not a tuning artifact**: γ=2/√N_corr is strictly decreasing in N_corr, and C is strictly increasing in γ (for ρ>ρ_crit), so C is strictly decreasing in N_corr. The most quantum-coherent macroscopic systems (largest N_corr) are pinned at C≈0 by the equation.

**Three things the names invoke, vs what C does:**
1. **Quantum phase coherence** (ODLRO) is *maximized* in BEC/BCS — exactly where C≈0. C is anti-correlated with quantum coherence.
2. **Synchronization** — a marching band and a BEC are both maximally collective/synchronized → both large N_corr → both low C. **C is anti-correlated with synchronization, the framework's namesake.** A framework called *Synchronism* has a central variable that is smallest for the most synchronized systems.
3. **The framework's own defining metaphor** — landing page/glossary define high coherence by "electrons in a superconductor, in lockstep." But a superconductor (large N_corr → tiny γ) is pinned at C≈0 by the equation. **The exemplar used to define high coherence is a low-C system by the framework's own formula.** (The explorer notes `terms.ts` even contradicts itself: γ entry says "γ≪1 = quantum," N_corr entry calls a crystal at γ≈10⁻¹² "classical.")

**What C actually is** (per the equation, and consistent with S613 = no decoherence parameter, S652 = no dC/dt, the compander/γ-dual-role/three-C findings): a density-driven saturation index — high for *dense, weakly-correlated* matter, low for *sparse or strongly-correlated* matter — carrying no ℏ, temperature, action, or decoherence rate, hence **zero quantum-vs-classical content.** This is why the explorer's "rename to classicality / decoherence-fraction" options are *wrong*: they would convert a confusing label into an *asserted falsehood* (that a lone electron is "maximally classical," a BEC "maximally quantum"). The honest fix is a **scope statement**, not a vocabulary swap (explorer Options I/II): either restrict C to the single-particle/galaxy regime where it was calibrated, or keep multi-scale and state explicitly that C does not measure quantum-vs-classical character.

**Why this is a foundational tension, not a site nitpick (Tension #5 — what is the framework protecting?):** the protected, fact-looking assumption is that **C measures coherence/synchronization.** It does not — by its own equation it is *anti-correlated* with both. The name, the namesake of the whole framework, and the defining metaphor all point one way; the equation points the other. This is the lossy multi-axis→one-scalar projection (compander class, γ dual-role, three-C — the whole arc) made visible at the front door, in the central word.

## Problem 2: TEST-03/05 Double-Filing (endorsed, connected to the archive)

The proposal flags that the site files one environment-dependence result as both "FAILED" (R²=0.14 < 20% kill criterion) and implicitly "discovery" (p=5×10⁻⁶, no kill badge), letting a reader reach opposite conclusions. **The resolution is correct and is a clean significance-vs-effect-size point:** p=5×10⁻⁶ and R²=0.14 are *one* signal — statistically significant but explaining 14% of variance, below the pre-registered 20% effect-size kill threshold. A tiny effect can be arbitrarily significant with enough data; the kill criterion was wisely set on effect size, so **R²=0.14 is a failure regardless of p-value.** Both cards must give that verdict.

This is the *site's* TEST-03/05; in the **archive** this is the environment-dependence result I executed in **S637** (RAR σ_int environmental slope ≈ 0.00016 dex, ~120× below the SPARC scatter floor) and **S654** (environment tests MOND+EFE-degenerate). So the archive already reached "fails by effect size"; the site's double-filing is a presentation inconsistency, not a new physics question. (Note the recurring **site↔archive numbering discrepancy** I flagged in S674/S675: site TEST-03/05 = TFR-scatter/RAR-environment ≠ archive TEST-03/05 = compact-elliptical/CMB-cold-spot. The numbering should be reconciled.)

## Self-Check (SESSION_PRIMER STOP list)

- **Compute, don't assert** (S669/S675): the inversion is computed across a ladder at three densities and shown monotone + structural, not argued.
- **Credit honestly**: the Problem-1 adjudication is the explorer's (2026-05-27) reconciling the 2026-05-04 γ-collision proposal; my contribution is the numerical verification against the equation, the synchronization/namesake point, and the connection to S613/S652/compander arc.
- **Consensus-attractor watch** (Tension #2): invoking ODLRO/BEC-as-maximally-coherent is *not* an imported external standard used to score the framework — it is the framework's *own* defining exemplar (superconductor lockstep), so this is an internal-consistency check, not a comfort reach.
- **Not the treadmill**: responded to genuine new input; did not grind another frontier-test provenance check.

## Files

- `Research/Session676_Coherence_Naming_Inversion_Verified.md` (this document)
- `simulations/session676_coherence_inversion.py` (C across the coherence-ordered ladder; monotone-decreasing verification)

## So What?

The framework is named *Synchronism* and its central variable is named *coherence*, defined by a lockstep/superconductor metaphor. Its own equation makes that variable **anti-correlated with both synchronization and quantum coherence**: the lone electron scores C≈1, the BEC scores C≈0. This is not a terminology quibble — it is the clearest single demonstration of what the whole audit arc has found, surfaced in the framework's own name: C is a density-saturation index dressed in the vocabulary of coherence/synchronization it does not measure. The right response is a scope statement (it is a galaxy-regime density map), not a rename (renaming to "classicality" would assert a falsehood, since C has no quantum/classical content at all — S613). And Problem 2 confirms the recurring discipline: a small-but-significant environmental effect (R²=0.14) is a *failure* by the pre-registered effect-size criterion, exactly as S637 already found by execution.

Cumulative: 40 audit/governance (S676 = coherence-naming inversion verified + TEST-03/05 double-file resolution) + 1 executed refutation (S661) + 1 post-hoc disfavoring kill-triggered (TEST-04a, S672) + novel-survivor 0 + 2 foundational-tension proofs (S665/S666) + 1 synthesis (S667) + 1 executed null (S669) + 1 method-specificity test (S670) + 1 frame resolution (S671). Frontier still 3/9 settled (not-derived), 6/9 unverified.
