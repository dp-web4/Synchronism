# Proposal: the unfalsifiable-badge class — where refutations hide on an honesty-branded site

**Date**: 2026-07-14
**Track**: Explorer (synchronism-site)
**Status**: proposal — one mechanical rule, two executed instances
**Evidence**: `synchronism-site/explorer/findings/2026-07-14-btfr-bounded-boost-refutation.md`

---

## The finding that forced this

`/tier-1-existing` badged **TEST-09 (BTFR)** as **"MOND-Shared"** with the note: *"A positive result
is consistent with Synchronism AND standard MOND — it cannot discriminate between them."* It had
never been run.

Executed on real SPARC (2026-07-14): the framework's own bounded law predicts BTFR slope
**n = 3.35 ± 0.07** against an observed **3.75 ± 0.10** — **3.3σ**, and its *registered* kill
criterion (slope deviation > 0.3) **fires at 0.41**. MOND passes at 0.6σ. It is not MOND-shared; it
is the sharpest discriminator in the Tier-1 ledger and the framework loses it.

**TEST-10** carries the same badge, was also never run, and is *structurally* impossible: it
predicts apparent DM fraction → 100%, while the framework's boost ceiling (1/Ω_m = 3.17) caps that
quantity at **68.5%**. No data required.

**2 of 2 MOND-Shared badges audited are concealed refutations.**

---

## The mechanism, stated as a rule

> **A badge that asserts a *tie* is unfalsifiable, and therefore terminates tests without executing
> them.**

Every other badge in the taxonomy either **commits to an outcome** (Failed, Kill Criterion
Triggered, Reparametrization, Refuted) or **admits ignorance** (Untested, Speculative). Both classes
are honest: the first can be wrong and checked; the second says "we don't know."

"MOND-Shared" does neither. It says *we tie*. A tie:
- requires no execution,
- produces no number,
- cannot be contradicted by data that was never collected,
- and **sounds like a result** — indeed sounds *modest*, which is why it survives an honesty audit.

That last property is the load-bearing one. On a site whose entire epistemic immune system is tuned
to catch **over-claiming**, a badge that reads as a concession is invisible. But "we tie with MOND
on the BTFR" is not a concession — **it is a concealed win claim**, because tying with MOND on the
BTFR would be a genuine success, and the framework does not tie.

---

## Why the existing audit rules didn't catch it

This is a third error direction, and the two rules proposed on 2026-07-14 (morning, maintainer) do
not cover it:

| Failure mode | Where it lives | Caught by |
|---|---|---|
| Over-refutation in verdict statistics | TEST-03, TEST-04a | citation-walk; criterion-verdict substitution rule |
| Over-promotion in explanatory prose | critical exponents, γ/B collision | analogy citation-walk rule |
| **Over-claiming *inside a verdict*** | **TEST-09, TEST-10 (MOND-Shared)** | **nothing — until now** |

The 2026-07-09 reflexivity finding established that **self-referential statistics break far more
often than physics statistics** (SELF 4/4 vs physics ~4/27): the program is *harder on itself than
on its theory*. TEST-09 is the counterexample that completes the picture. Where the program is
**soft** on itself, it is soft in exactly the place the reflexivity predictor says to stop
looking — a verdict that already sounds self-critical.

Both biases share one root: **the audit checks the direction of a claim's sentiment, not whether
the claim was ever tested.** A pessimistic-sounding untested claim passes; an optimistic-sounding
tested claim gets walked. The correct axis is not sentiment. It is **execution**.

---

## Proposed rule (mechanical, not judgment)

**Any badge asserting a tie, a degeneracy, or a shared prediction must carry the same execution
burden as one asserting a kill.** Specifically, a "shared with [rival theory]" badge is admissible
only if:

1. The two theories' predictions for the observable are **computed** — not argued from the
   theories' verbal descriptions — and
2. The computed predictions are shown to **agree to within the data's discriminating power**.

If either is missing, the correct badge is **Untested**, not "shared."

This is cheap to enforce and would have caught both instances immediately: neither TEST-09 nor
TEST-10 had a computed Synchronism prediction at all. The badge was inferred from a *verbal*
resemblance to MOND ("regime-dependent slope is a textbook MOND signature") — and the verbal
resemblance is false, because MOND's boost diverges and Synchronism's saturates.

### Generalization beyond this program

The rule is transferable to any research program that benchmarks itself against an incumbent theory.
The failure is not "we cheated"; it is that **degeneracy claims are the only claims that are cheaper
to assert than to test**, so they accumulate silently wherever a program is graded on honesty rather
than on execution. A ledger of "we tie here, we tie there" can look like scrupulous modesty while
being the least-audited surface in the program.

---

## Immediate consequences for this program

- **TEST-09 and TEST-10 move to Bucket 2 (refuted).** Done — PREDICTIONS.md updated 2026-07-14.
- **Every remaining MOND-Shared badge must be audited.** The general criterion is now known:
  *any observable sensitive to the asymptotic boost discriminates*, because boundedness is (per the
  site's own `/galaxy-rotation`) "the only form whose prediction differs from MOND."
- **The boost-ceiling refutation must be linked to the tests it refutes.** `/honest-assessment`
  already calls it "the strongest direct refutation in the framework's own internal audit." It sat
  two pages away from a badge asserting non-discrimination on the observable it most directly kills.
  The refutation was never *false* — it was never **connected**.

---

## The deeper structural result (worth a paper, not just a ledger row)

Stated without reference to Synchronism at all:

> **A modified-gravity theory whose boost function is bounded cannot reproduce the BTFR.**
> Asymptotically flat rotation curves require `g_obs ∝ g_bar^(1/2)` exactly — the unique exponent
> for which `V² = g_obs·r` is r-independent. Any bounded boost is asymptotically a constant
> rescaling of G, i.e. Newtonian, giving BTFR slope n = 2. The BTFR's observed n ≈ 3.85 therefore
> **excludes the entire class of bounded-boost modifications**, independently of parameters.

Synchronism is one instance. The result is the transferable content, and it is *cleaner than the
locality no-go* (the program's current #1 transferable null) because it needs **no data at all** to
state — boundedness alone forces n = 2, and the data only says by how much that misses.

This gives the proposed preprint a **second, independent structural blade**:
- **Locality no-go** — local ρ cannot key a non-local g_bar threshold (a quantified Milgrom
  non-locality instance).
- **Boundedness no-go** — a bounded boost cannot produce a BTFR (new, 2026-07-14).

Both are one-line statements with worked refutations, and both kill *classes*, not just this
framework. Recommend folding into `locality-nogo-standalone-writeup` and retitling to reflect two
blades. **Decision gates on dp.**
