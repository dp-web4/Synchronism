# Triage — a₀(z)∝H(z) engaged by data and disfavoured (new TEST-25 catalog row) + DESI DR2 prospectivity timing flag (2026-08-01)

**Status:** `[ACTIVE-MRH]` — gate-fired on two site/publisher proposals. **Verdicts: (1) the a₀(z)∝H(z)
prediction is genuinely engaged by outside data for the first time and DISFAVOURED at 3.2× (single-source);
verified the in-repo pieces, added TEST-25 to the catalog with honest post-hoc scoping. (2) The DESI DR2
prospectivity claim has a real, open timing-verification risk I cannot close in-repo; flagged it caveated on the
ledger. Bucket 0 unchanged (0); refutation count unchanged (6); arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Why this is a real gate-fire, not boost-ceiling churn

Yesterday I parked the boost-ceiling headline as over-worked. This is the opposite case: a**genuinely new
prediction engaged by new data for the first time** (a₀(z)∝H(z), a MOND-scale evolution claim the framework
carries in whitepaper §5.15 with zero free parameters). It clears the "productive vs churn" bar cleanly.

## (1) a₀(z)∝H(z) — engaged, disfavoured; TEST-25 added

Verified the in-repo, checkable pieces:
- **The prediction is real and in-repo:** whitepaper §5.15 line 72 — "a₀ is cosmologically determined … Predicts
  evolution with redshift: a₀(z) ∝ H(z), testable via high-z BTFR." (The whitepaper also already carries a
  thorough, well-caveated 08-01 status note on the disfavouring — this is not new to the whitepaper, only to the
  *catalog*.)
- **The arithmetic checks (by hand):** for a₀(z)∝E(z), da₀/dz|₀ = a₀(0)·dE/dz|₀ = a₀(0)·(3Ω_m/2) = 0.491×10⁻¹⁰
  (Ω_m=0.315, a₀(0)=cH₀/2π≈1.04×10⁻¹⁰) — matches the quoted 0.49. Ratio to the measured a₁=1.59×10⁻¹⁰ is
  **3.24× (steeper)** — matches the quoted 3.2×. So the framework predicts a₀ evolving *too slowly*.
- **What I could NOT verify:** the Ciocan et al. 2026 measurement itself (external paper, not in-repo) — so
  a₁=1.59 travels as the paper's reported value, caveated.

**The honest status is ENGAGED — DISFAVOURED, not refuted**, and I inscribed it that way: single-source; the
external literature is *mutually inconsistent* (Gueorguiev 2024 finds zero slope, Milgrom 2017 disfavours faster
evolution — a₀(z)∝H(z) is squeezed from both sides by measurements that disagree with each other); Ciocan's
systematics are large. The **kill criterion is POST-HOC** (the measurement predates any pre-commitment) — marked
as such in the row, and explicitly *not counted* toward the catalog's "criteria fired and honored" tally. **The
refutation count is unchanged (6).**

**The deeper, honest framing (which I preserved from the whitepaper note):** the result is *worse* for the
alternative anchor 2πa₀≈c²(Λ/3)^{1/2} (which predicts no evolution at all) — so **what is embarrassed is the
general programme of deriving a₀ from cosmology, not Synchronism specifically.** This is the dynamical extension
of the a₀=cH₀/2π Bucket-3 reparametrization (my 07-05 entry): the reparametrization now has a *dynamical*
prediction that has been engaged and disfavoured.

**The catalog gap, and the class flag.** The catalog (24 experimental + 3 theoretical tests) had no a₀(z) row
despite the whitepaper stating the prediction — so the catalog's "no registered test remains runnable-and-unrun /
complete inventory" reading was misleading. Added TEST-25 (Tier 1) to close the instance. **The class point
matters more (proposal item 4):** the catalog appears assembled from *arc outputs*, not swept from the
whitepaper's own testable claims — one gap found by accident implies others. I flagged a one-pass
whitepaper-vs-catalog diff as recommended and unrun; I did **not** execute the full sweep this session (it is a
larger task, and executing it well means reading every whitepaper "predicts/testable" claim — a scoped job worth
its own gate-fire, not a bolt-on).

## (2) DESI DR2 prospectivity — a real timing risk I cannot close

The DESI DR2/TEST-04a registration (the one I routed to dp, adopted 07-17) is called "the program's first
genuinely prospective test." A site/maintainer integrity check flags that a **DESI DR2 talk (PIRSA:26040071)
was given April 2026, three months before the 07-17 adoption** (and the fσ₈≤0.46 threshold was git-committed
07-01/07-14 — also after the talk). "Prospective" is a claim about what was *knowable*: if that talk showed a
threshold-relevant fσ₈(z≈0.5) figure, the prospectivity claim is compromised — the same failure mode the program
already caught once on this exact test (TEST-04a's original post-hoc σ₈ calibration).

**I cannot resolve this** — it requires watching/reading the PIRSA recording, which is not in-repo — so per the
don't-assert-what-you-can't-re-execute discipline I flagged it on the ledger as an **OPEN timing-verification
item**, caveating the "genuinely prospective" claim rather than either affirming or retracting it. The registered
statistic and verdict branches are unaffected; only the *prospectivity* claim is at risk, and it stays flagged
until someone with recording access confirms the talk showed no threshold-relevant number.

## Disposition

- **Research/EXPERIMENTAL_TEST_CATALOG.md** — TEST-25 added (a₀ redshift evolution, ENGAGED-DISFAVOURED, post-hoc
  kill visibly marked, count unchanged, sibling-audit flagged).
- **PREDICTIONS.md** — DESI prospective-registration entry caveated with the open timing-verification item.
- **Whitepaper §5.15** already carries the disfavouring note (maintainer lane) — verified consistent, not edited.
- **Flagged / not executed:** the whitepaper-vs-catalog sibling diff; the PIRSA-talk content check (needs
  external access).
- **Bucket 0 unchanged (0); refutation count 6; arc AT REST.**

## So what

For the first time in the arc, a framework prediction was engaged by *new outside data* rather than re-analysed
in-repo — and the honest verdict is a disfavouring, not a win: a₀(z)∝H(z) evolves too slowly by 3.2× against the
one measurement, which is itself single-source and contradicted by two others. The right reading is the one the
whitepaper already reached — this embarrasses *deriving a₀ from cosmology in general*, not Synchronism uniquely,
and it is a disfavouring, not a refutation. The catalog gap it exposed is the more durable finding: a registry
built from arc outputs will silently omit the framework's own stated predictions, and the completeness claim then
reads as more than it is. Meanwhile the DESI prospectivity claim I helped inscribe now carries a timing risk I
can flag but not close. Bucket 0 stays 0.
