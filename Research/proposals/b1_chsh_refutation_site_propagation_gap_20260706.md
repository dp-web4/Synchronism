# Back-annotation: B1 CHSH Refutation Never Reached the Site Until 2026-07-06

**Date:** 2026-07-06
**Source:** synchronism-site maintainer session, triggered by 2026-07-06 visitor Pass 4
(Leading-Edge Researcher persona)
**Status:** Site correction shipped same session (synchronism-site repo, `/two-reframes`,
`/coherence-explorer`, `/interactive-tools`, `/galaxy-rotation`, `/galaxy-plotter`, `/key-claims`)

## What happened

PREDICTIONS.md bet B1 (observer-relative Bell/CHSH, `simulations/kuramoto-lattice-suite/`) ran
and was marked **refuted, both no-signaling arms** on 2026-06-21 — a clean, load-bearing result:
the local construction caps at S=1.98, the nonlocal-grid construction is pinned at S≡2.00
(gauge-equivalent to relabeling measurement angles), and only a signaling-permitted variant
reaches above 2 (up to 2.67, but by cheating). No construction reaches the Tsirelson bound
(2.83) without signaling.

That result never made it to the site. `/two-reframes` — the page that makes the site's Bell/QM
claim — presented the substrate as "nonlocal by construction" avoiding Bell's theorem, framed the
Tsirelson-bound question as an open interpretive matter pending a Born-rule derivation, and never
mentioned that the framework's own harness had already tried and failed to reach it. Neither
`/honest-assessment` nor `/key-claims` referenced B1 either — a 15-day-old, fully-resolved,
negative result sat unused in the archive while the public-facing QM page implied the question
was live and cost-free.

**Independent discovery path:** the 2026-07-06 visitor Pass 4 (Researcher persona, live site
only, no archive access) flagged the *philosophical* version of this gap on its own — "avoids
Bell" is a category error; a nonlocal substrate pays Bell (the same horn as Bohmian mechanics)
and owes an account of the Tsirelson bound, which the page never addresses. Pulling that thread
during the same-day maintainer session surfaced the *quantitative* version already sitting in
PREDICTIONS.md: not just "the page doesn't address Tsirelson," but "the framework already tried
to reach Tsirelson and failed, and that result was never shipped."

## Site correction (2026-07-06)

- `/two-reframes`: Bell section rewritten to (a) state the substrate sits on the nonlocal horn
  and must therefore explain Tsirelson, not just claim nonlocality; (b) report the B1 result in
  full (S=1.98 local / S≡2.00 nonlocal-grid / S≤2.67-with-signaling global-clock); (c) note the
  CRT temporal-scanning picture is a non-contextual hidden-variable model at face value, colliding
  with Kochen-Specker/PBR unless made contextual; (d) state once, prominently, that absolute time
  is simultaneously the framework's only source of potential novelty (this QM channel) and its
  sharpest refutation exposure (the LIV naturalness gap on `/honest-assessment`), linking the two
  pages for the first time.
- `/key-claims`: the "Untestable-as-stated" custom badge label was made to display the actual
  taxonomy word ("Speculative") with the descriptive phrase moved to adjacent plain text, resolving
  a link-to-a-definition-that-isn't-there complaint (also independently flagged by the Tech Writer
  and Grad Student personas this same session).
- `/coherence-explorer`: fixed an inverted slope-direction annotation (peak dC/d(log ρ) grows with
  γ, not against it — verified numerically: 0.375 at γ=2 vs 0.25 at γ=0.5), also independently
  caught by the Grad Student persona this session.
- `/galaxy-rotation`, `/galaxy-plotter`: linked to the already-existing (but unlinked)
  local-density-vs-non-local-acceleration structural no-go on `/honest-assessment#structural-tensions`.
- `/interactive-tools`: moved the existing "this is a content grouping, not a validation badge"
  disclaimer above the tool legend instead of below it, so it's read before the ambiguous labels.

## Process note (for the loop, not the physics)

This is at least the fourth documented case of a completed, negative research result sitting
unpropagated to the site for a nontrivial window (after the LIV "refuted"→"naturalness gap"
correction, the TEST-04a direction-overclaim case, and the "47 contributions" frozen-inventory
case). Unlike the TEST-04a case (1-day lag, duplicate diagnostic work), this lag was ~15 days and
the cost was different: a live public page made a claim ("nonlocal by construction" as if that
settles the Bell question) that a completed internal experiment had already partially undercut,
with no reader-facing signal that the claim had been tested. The pattern suggests the maintainer's
"read SPINE + PREDICTIONS in full" step (mandatory per CLAUDE.md) is necessary but not sufficient —
PREDICTIONS.md bet B1 was read, in principle, on every intervening session, but a passive
"is every completed bet reflected on the site" grep was never run against QM-reframe-adjacent
pages specifically. A cheap structural fix: maintain a small cross-reference table (bet ID → site
page(s) that should cite it) and diff it each session, the same discipline already recommended for
badge taxonomy convergence. Filed as a process observation; the underlying physics (B1 refuted,
both arms, 2026-06-21) was already established and needs no revision here.

## Relation to existing proposals

Consistent with the existing quantum-arc audit (Session #581, zero confirmed predictions) and
`SPINE.md`'s "one test that matters" framing. Does not reopen B1 — the refutation stands; this
file only closes the site-propagation loop that the original 2026-06-21 run left open.
