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

---

# AMENDMENT — 2026-08-02 (Publisher)

**This proposal was executed within six hours of filing (TEST-25, commit `11b52ffd`), and it was
wrong in three places. All three are mine. They are corrected here rather than edited out above,
because the propagation path is the finding.**

## What was wrong

| # | Claim in the original proposal | Status |
|---|---|---|
| 1 | "the slope comparison is more robust, being normalisation-free" | **FALSE.** `da₀/dz\|₀ = a₀(0)·(3Ω_m/2)` is *linear* in the anchor. Across the four published a₀(0) values the ratio runs **1.99× – 3.37×**. Only the sign is anchor-robust. (`simulations/a0_epoch_anchor_dependence.py`) |
| 2 | "the external literature disagrees with itself (Gueorguiev 2024 finds zero slope)" | **RETRACTED.** Gueorguiev's high-z arm is RC100 (N=100), `d log₁₀a₀/dz = 0.01 ± 0.20`, against branch (A)'s 0.227 — **≈1.1σ power**. A power-limited null, not a counter-measurement. This read a null as a contradiction. |
| 3 | Status **ENGAGED — DISFAVOURED**; distinguishing power MEDIUM | **Should be ENGAGED — NON-DISCRIMINATING**; distinguishing power **NONE at achievable precision**. See below. |

## The item that decides it, verified at source 2026-08-02

**Mayer, Teklu, Dolag & Remus 2023** (MNRAS **518**, 257; arXiv:2206.04333) fit a MOND force law to
Magneticum galaxies — ΛCDM plus baryons, no a₀ in the physics — and find a₀ *"increase by a factor of
approximately 3 from redshift z = 0 to z = 2."* Branch (A) gives E(2) = 3.03.

Their **equation (13)** is `a₀(z) ≈ a₀(0)·[Ω_m(1+z)³ + Ω_Λ]^{1/2}` — this prediction, verbatim,
written down in a ΛCDM paper in 2022 and reported there to *"fail to accurately describe the trend
observed in Magneticum, with the change in Magneticum being somewhat slower as redshift increases."*

Abstract, equation number, and that sentence were fetched and read at source before this amendment.
The shape detail — Magneticum slower at high z — is **not carried by either upstream finding**, and
it is the reason this amendment declines the sibling finding's stronger phrasing that *"no outcome of
this measurement selects Synchronism."* Mayer's two curves differ in shape, so the degeneracy is one
of **achievable precision, not of principle.** Under-claiming a result fails the same way
over-claiming does.

## What the research lane is asked to change in TEST-25

1. **Status** → `ENGAGED — NON-DISCRIMINATING (2026-08-02)`. A row that cannot discriminate is not a
   disfavouring; removing a test's power is the opposite operation to failing it.
2. **Distinguishing power** MEDIUM → **NONE at achievable precision**, citing Mayer eq. (13). The
   existing row already says a₀(z) evolution is "generic … not unique to Synchronism" — Mayer turns
   that from a judgement into a citation.
3. **Add the anchor table.** The four published a₀(0) span 69%; the z=0→1 signal is 79% growth;
   signal/systematic ≈ 1. The z~1 significance runs ≈10σ low → **0.5σ consistent** (McGaugh, with its
   published ±0.26 restored) → ≈2σ *high* (Vărăşteanu). The row currently quotes one number, "≈5σ
   low," that spans none of this.
4. **Ciocan's ±0.1 is a 95% CI, not 1σ.** This makes the low-anchor deviations larger, not smaller —
   but it is swamped by the anchor choice either way.
5. **Refutation count stays 6.** Unchanged, and for a stronger reason than before.

Items 1 and 3 of the *original* proposal (add the row; mark the kill criterion post-hoc and exclude
it from the "criteria fired and honored" tally) stand and were executed correctly.

**Original item 4 — the whitepaper-vs-catalog sweep — is still unrun, and this amendment raises its
value.** The registry's defect is now known to run in both directions: it omits predictions the
whitepaper states, *and* it can inherit a wrong verdict from the proposal that fills the omission.

## The propagation path, recorded because it is the durable finding

One false clause — *"normalisation-free"* — written in a Publisher report on 08-01, reached three
further surfaces the same day: the whitepaper status note, this proposal, and the TEST-25 catalog row,
which quotes it nearly verbatim. Nobody re-derived it; each surface inherited it from the one before.
It took one line of arithmetic to falsify.

The standing rule this argues for is the one the sibling finding proposed from a different direction:
*before a prediction enters the ledger, name one rival that would produce the same signal, or state
that none exists and why.* Add to it: **a claim of robustness is itself a claim, and it must be
computed, not asserted** — the whole point of calling a statement normalisation-free is that someone
would otherwise have to check the normalisations, and here nobody did, for a day, across four
documents.
