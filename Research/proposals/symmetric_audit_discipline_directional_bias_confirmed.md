# Symmetric Audit Discipline: The Site's Review Process Has a Confirmed Directional Bias

**Date**: 2026-07-09
**Source**: Site maintainer session, synthesizing 2026-07-08 explorer citation-walk +
2026-07-09 four-persona visitor pass (WAKE-phase proposal, before site fixes)
**Status**: Proposal — a methodology/process finding, not a physics finding

## The pattern, now confirmed across two independent tracks on two consecutive days

The explorer track's 2026-07-08 citation-walk found the TEST-03 "kill criterion triggered"
claim was manufactured (five separate provenance errors, all in the direction of making the
framework look *worse* than the data support). The visitor track's 2026-07-09 four-persona
pass, run the next morning with no knowledge of the explorer's findings, independently
re-derived the same TEST-03/05 conflation from scratch by direct arithmetic (t≈48.7 at
N=14,585 implies p~10⁻⁵⁰⁰, not the quoted 5×10⁻⁶) — and found three *more* instances in the
same session: a Σ₀ arithmetic error (110/12% reported, 119-123/4-0.5% actual — always
understating the match), a Bullet Cluster refutation aimed at a claim the framework doesn't
make (a strawman that makes the failure look sloppier than the real, stronger argument would),
and an A2ACW null denominator inflated 70-550× on the one page written for citers (the site's
only *methodology* overclaim, which is the mirror image of the physics under-claims).

Tally as of 2026-07-09: **6 independent provenance-break instances found, 6 over-refute the
framework's physics, 0 over-claim the physics.** The one over-claim found (A2ACW denominator)
is on the *methodology* track, not physics — meaning the physics-audit machinery is
systematically asymmetric, and the methodology-audit machinery doesn't exist at all.

## Why this is a research-direction finding, not a site-content finding

This is not "the site has some bugs." It is a property of *how corrections get generated and
propagated* across the visitor → maintainer → explorer loop: the loop's implicit incentive
structure rewards finding overclaims (the site's stated mission — "researcher, not lab worker,
question the frame") and has no symmetric mechanism that rewards finding *under*-claims. A
program that scores itself only on one axis will drift on that axis, even — especially — when
the drift direction is "more honest than the truth warrants." An audit trail that always finds
the framework guilty is not more rigorous than one that always finds it innocent; it's the
same failure mode with the sign flipped, and it is currently invisible to a review process
built entirely around catching the other sign.

## What's already been done (site-side, 2026-07-09 maintainer session)

Fixed on the live site: TEST-03/05 conflation (9 carrier pages), Σ₀ arithmetic (4 pages), the
Bullet Cluster strawman (2 pages), the A2ACW denominator inflation (3 pages), a mis-attributed
LIV citation (also fixed in this repo's `PREDICTIONS.md`), and the Claim 3 badge reasoning
(now grounded in the mass-cancellation locality no-go instead of the retracted statistic).

## What's still open — the process fix, not the content fix

1. **Adopt a standing rule**: any "kill criterion fired" or "prediction refuted" claim must
   walk to (a) the registered threshold text, (b) its registration date vs. measurement date,
   (c) the variable actually measured — *before* it ships, not after a visitor catches it.
   (Already proposed once, in `test03_site_kill_manufactured_and_test08_unrun.md`; this
   proposal extends it from "apply to failure claims" to "apply symmetrically to methodology
   claims too.")
2. **A2ACW self-audit**: the methodology track has never had its own claims (denominators,
   sensitivity without specificity, no control arm) checked with the same rigor the physics
   track receives. This proposal recommends the explorer track run one dedicated citation-walk
   on `/a2acw`, `/research-philosophy`, and `/for-researchers` artifact 2 specifically —
   treating the methodology's own claims as a test subject, not just the audit apparatus.
3. **Track the tally itself.** If the 6/6-over-refute pattern continues (i.e., stays
   asymmetric) after this correction round, that asymmetry is itself worth a session:
   something about how corrections get generated is structurally one-directional, and fixing
   today's six instances won't fix the generator.

## Decision for dp

Is this worth a standing process change (e.g., "every maintainer session's fix list must
include at least one re-derivation of a *negative* verdict, not just positive-claim checks"),
or is six instances over roughly six weeks within the noise of an otherwise-working
self-correcting loop? The loop caught all six eventually — the question is whether "eventually,
via an alert visitor persona" is good enough, or whether the asymmetry should be designed out.
