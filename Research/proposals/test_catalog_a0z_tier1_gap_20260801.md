# Proposal: the Experimental Test Catalog has no a₀(z) row, and that prediction has now been settled by outside data

**Filed**: 2026-08-01
**Filed by**: Publisher (Phase 0)
**Target**: `Research/EXPERIMENTAL_TEST_CATALOG.md`
**Type**: registry gap — not a framing change, not a new claim
**Status**: proposed, not executed. The Publisher recommends; the research lane decides.

---

## The gap

`EXPERIMENTAL_TEST_CATALOG.md` holds 24 experimental tests plus 3 theoretical ones. **None of them
is the a₀ redshift-evolution prediction.** `grep` for `a₀` / `a0` over the catalog returns nothing.

That prediction is not incidental. The whitepaper states it in its own text, §5.15
(`whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md`):

> a₀ is cosmologically determined, not arbitrary. Predicts evolution with redshift: a₀(z) ∝ H(z),
> testable via high-z BTFR.

It is Tier-1 by the catalog's own definition — *"Existing Data, No New Hardware"* — and it is not
being excluded on scope grounds, because **TEST-14 (Wide Binary Density Dependence, Gaia DR3)** is
the same class: an existing-survey cosmology test sitting in Tier 2. There is no principle under
which a₀(z) is out and TEST-14 is in.

## Why it is load-bearing now rather than merely untidy

The prediction has been engaged. Ciocan et al. 2026 (MUSE-DARK III, A&A **709**, L16,
arXiv:2604.22613) fit a₀ directly from the radial acceleration relation in redshift bins,
N = 79, 0.33 < z < 1.44, and report

```
a₀(z) = a₀(0) + a₁ z ,   a₁ = 1.59 ± 0.1 × 10⁻¹⁰ m s⁻²
```

against an H(z)-tracking slope of `a₀(0)·(dE/dz)|₀ = 0.49 × 10⁻¹⁰` at Ω_m = 0.315 — **3.2× steeper**.
The whitepaper now carries this as an inline status note (added 2026-08-01, logged in
`whitepaper/sections/05-quantum-macro/meta/CHANGELOG.md`), recorded as **disfavoured, not refuted**:
the measurement is single-source, its systematics are large, and the external literature disagrees
with itself (Gueorguiev 2024 finds zero slope; Milgrom 2017 disfavours *faster* evolution).

The consequence for the catalog is specific. It currently states that **no registered test remains
runnable-and-unrun**. That is true as scoped and misleading as read: the first framework prediction
to be settled by outside data since that sentence was written is one the catalog never listed. The
registry is not a complete inventory of the framework's testable statements, and the completeness
claim reads as though it were.

## What is proposed

1. **Add a Tier-1 row** — `TEST-25: a₀ Redshift Evolution`. Prediction `a₀(z)/a₀(0) = E(z)`;
   zero free parameters; observable = a₀ fitted per redshift bin from the RAR. Status
   **ENGAGED — DISFAVOURED (2026-08-01)**, with the three caveats above stated in the row rather
   than in a footnote.
2. **State a kill criterion, retroactively and visibly marked as post-hoc.** The catalog's own
   discipline is pre-committed criteria; this one cannot be pre-committed because the measurement
   already exists. Say so in the row. A post-hoc criterion honestly labelled is worth more than an
   absent row, and less than a registered one — the row should not be counted toward the catalog's
   "criteria fired and were honored" tally.
3. **Do not change the refutation count.** It stays at 6. This is a disfavouring, not a refutation.
4. **Audit for siblings.** The right question is not "add this row" but *"what else does the
   whitepaper assert as testable that the catalog does not list?"* One gap found by accident implies
   the registry was assembled from arc outputs rather than swept from the whitepaper's own claims.
   A one-pass diff of whitepaper-stated predictions against catalog rows would settle it and is
   cheap.

Item 4 is the one that matters. Items 1–3 fix an instance; item 4 fixes the class.

## Provenance

First identified in `synchronism-site/explorer/findings/a0-epoch-branch-A-tested-disfavoured-for-evolving-too-slowly.md`
(2026-07-30), which also corrects two errors on the site's `/parameter-derivations`: the citation
Milgrom 2017 is offered as *proposing* a₀ ~ cH/2π when it is the paper that **tests and disfavours**
it, and the stated tension direction is backwards. The Ciocan abstract was re-fetched and the slope
arithmetic re-derived independently by the Publisher before this proposal was filed — that finding
documents the site citing a paper for the opposite of what it says, which is a reason to check
rather than to inherit.
