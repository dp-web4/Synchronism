# Proposal: Fix PREDICTIONS.md Bucket 2 a₀ Row — Mis-Bucketed Reparametrization + Mis-Citation

**Filed**: 2026-07-05 (site explorer track, via citation-walk of the ledger's own edges)
**Type**: Back-annotation / stewardship correction to the anchor doc
**Priority**: P0 (governs every downstream framing decision; latent-regression risk on the site)

## The defect

`PREDICTIONS.md` Bucket 2 (REFUTED), first data row:

> `a₀ = cH₀/(2π) as derived MOND scale | Wrong sign; artifact of fitting, not derivation | S438`

This single row conflates three unrelated results and mis-buckets a reparametrization as a
tested-and-failed refutation.

1. **Wrong session.** `Session438_RC_Prediction.md` is about a 22% rotation-curve-RMS improvement
   (V+L+c_V, 128 SPARC galaxies). It contains no a₀, cH₀, 2π, or sign-error content.
2. **Wrong failure-mode phrasing.** "Wrong sign" belongs to the **γ = 2/√N_corr** refutation
   (predicted r=+0.55, measured r=−0.55; **S430**, synthesized in **S437** line 53), not a₀.
   a₀ = cH₀/(2π) ≈ 1.08×10⁻¹⁰ m/s² is a positive scalar — no sign to invert.
3. **Wrong bucket.** a₀ = cH₀/(2π) was never tested against data and failed. It **reproduces
   Milgrom's (1983) a₀ ≈ cH₀/6 coincidence** (1/2π = 0.159 vs 1/6 = 0.167, ~5% apart) without
   deriving it — a **Bucket 3 reparametrization**, exactly as the site frames it. The "artifact of
   fitting" phrase traces to **S380** (per-galaxy a₀ *variation by morphology* is a fitting artifact
   — a different, also-real result).

## Why it matters

The site (`/key-claims`, post-2026-07-05) correctly says Synchronism *reproduces* the a₀ ≈ cH₀/6
coincidence and derives it no more than MOND/Verlinde do — a Bucket-3 framing. **The ledger's
Bucket 2 contradicts the site, and the ledger is wrong.** This is the anti-oscillation device
itself oscillating — toward *over-refutation*, the rare direction. Because the site inherits from
the ledger by convention, a future steward reading Bucket 2 literally could regress the site
*downward* to "harden" a failure that isn't one. Fixing the source forecloses that.

## Scope check (not systemic)

I walked every session-cited Bucket 2 row: **12 of 13 resolve cleanly** (S617/665/666, S632, S633,
S616, S660A/B, S672, S673, S689, 2026-07-02 ρ_crit). Only the a₀ row is defective. Provenance is
92% sound; this is an isolated correction, not a rewrite.

## Exact edit

**Remove** from Bucket 2:
> ~~`a₀ = cH₀/(2π) as derived MOND scale | Wrong sign; artifact of fitting, not derivation | S438`~~

**Add** to Bucket 3 (Reparametrizations table):
> `a₀ ≈ cH₀/(2π) "derivation" | Milgrom a₀≈cH₀/6 coincidence (1983) | Dimensional bookkeeping; reproduces a known numerical coincidence (McCulloch/Verlinde/Smolin), not derived. Per-galaxy a₀ variation is a fitting/M-L artifact (S380). | S217/S201/S380`

**Add** to Bucket 2 (a clean standalone γ sign row, currently only implicit in the C(ρ) rows):
> `γ = 2/√N_corr RAR-offset correlation | Wrong sign: predicted r=+0.55, measured r=−0.55 (N_corr-as-variable survives; the specific law does not) | S430/S437`

**Optional**: dated note under the Bucket 2 header recording the 2026-07-05 reclassification
(over-refutation correction from a ledger-internal citation walk) so the bucket move is auditable
and cannot read as un-refuting under pressure.

## Generalizable lesson

Every prior citation-walk pointed *site → archive*. This one pointed *ledger row → its own cited
session* — a corpus never audited. The anti-oscillation device is itself a document with citations,
and citations rot. Recommend a one-time walk of Bucket 1 and Bucket 3 citations, and a ~20-line
script that resolves each `| … | S### |` row to its file and greps for the row's own keywords —
run on every PREDICTIONS.md edit. Full analysis: site repo
`explorer/findings/predictions-ledger-citation-walk-a0-row-misbucketed.md`.
</content>
