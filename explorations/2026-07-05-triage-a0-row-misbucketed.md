# Triage — the ledger's own a₀ row over-refuted: mis-cited + mis-bucketed; fixed (2026-07-05)

**Status:** `[ACTIVE-MRH]` — triage of `Research/proposals/predictions_ledger_a0_row_misbucketed_20260705.md` (site-explorer citation-walk of the ledger's own edges, 2026-07-05), tripped the hold-checklist. This is a correction to **my own anchor doc** (PREDICTIONS.md). **Verdict: ACCEPT — verified the three defects independently and applied the fix. The Bucket-2 a₀ row was over-refuting (the rare direction): it mis-cited S438, mis-attached a "wrong sign" refutation that belongs to γ, and mis-bucketed a reparametrization as a tested-and-failed refutation. Discipline-4 respected: nothing is un-refuted — the genuine refutation (γ wrong sign) is preserved on its own correct row. Bucket 0 unchanged (0); arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Independent verification (all three defects confirmed)

The defective row was: `a₀ = cH₀/(2π) as derived MOND scale | Wrong sign; artifact of fitting, not derivation | S438`.

1. **Mis-cited (verified).** Session #438 is "Predicting Individual Rotation Curves" (2026-02-06,
   RC-RMS improvement on SPARC) — grep confirms it contains **no** a₀ / cH₀ / 2π / sign content. The
   S438 citation is simply wrong.
2. **Wrong-sign mis-attached (verified).** a₀ = cH₀/(2π) ≈ 1.08×10⁻¹⁰ m/s² is a **positive scalar —
   no sign to invert.** The genuine "wrong sign" result (predicted r=+0.55, measured r=−0.55) is the
   **γ = 2/√N_corr RAR-offset correlation** — grep confirms it lives in **S430** (Theory Revisit) and
   **S437** (Grand Synthesis). It was mis-stapled onto the a₀ row.
3. **Mis-bucketed (verified by arithmetic).** a₀ = cH₀/(2π) **reproduces Milgrom's a₀ ≈ cH₀/6
   coincidence** (1/2π = 0.15915 vs 1/6 = 0.16667, ~4.5% apart) without deriving it — a Bucket-3
   reparametrization (same status MOND/Verlinde/McCulloch/Smolin have), not a Bucket-2
   tested-and-failed refutation.

## The fix applied to PREDICTIONS.md

- **Bucket 2:** replaced the defective a₀ row with a clean, correct γ row — `γ = 2/√N_corr RAR-offset
  correlation | Wrong sign: predicted r=+0.55, measured r=−0.55 | S430/S437` — re-attaching the real
  refutation to the claim it actually refutes.
- **Bucket 3:** added `a₀ ≈ cH₀/(2π) "derivation" | Milgrom a₀≈cH₀/6 coincidence | reproduces a known
  coincidence, not derived; per-galaxy a₀ variation is a fitting/M-L artifact`, with a dated note that
  it moved from Bucket 2 as an over-refutation correction. (Origin-session provenance S217/S201/S380 is
  cited per the proposal's walk; I verified the load-bearing facts — S438-wrong, γ-sign-in-S430/S437,
  coincidence-arithmetic — but not those specific origin IDs, and flagged that in the row.)

**Discipline-4 held:** this is not un-refuting a real result. The γ "wrong sign" refutation is
*preserved* (moved to its own accurate row); only the a₀ *reparametrization* moves to where it belongs
(Bucket 3, matching what the site already says). The ledger is now more accurate, not softer.

## The pattern (worth noting): the ledger over-refutes, mildly and twice now

This is the **second over-refutation correction in ~10 days**, both on the framework's own strongest-
sounding negatives: my dim-4 c_μν "sharpest falsification" was walked back to a standard naturalness
problem (2026-06-30), and now the a₀ row moves refuted→reparametrization. The anti-oscillation device
is real and valuable, but it has a mild tilt toward **over-failing** — the rare direction the
recalibration (dp, symmetric standards) flagged and that my own "over-failing is as seductive as
over-claiming" lesson names. A decisive-sounding "refuted / wrong sign" is satisfying and slips in
without the scrutiny an over-claim would draw. The proposal's generalizable recommendation is sound:
**the anti-oscillation ledger is itself a document with citations, and citations rot** — a one-time
walk of Bucket 1/3 citations + a small `| … | S### |`-resolver run on every PREDICTIONS.md edit would
catch this class. (Tooling recommendation is maintainer/operator lane; flagged.)

## Disposition

- **Ledger fixed** (Bucket 2 γ row corrected; a₀ moved to Bucket 3). Scope check accepted: the
  proposal walked all 13 Bucket-2 rows; 12 resolve cleanly, only a₀ was defective — isolated, not a
  rewrite. **Bucket 0 unchanged (0).**
- **Maintainer/operator lane:** the citation-resolver tooling recommendation; any site propagation.
- **Arc AT REST.**

## So what

The signal lifted the hold and pointed the citation-walk *inward* — at the ledger's own edges, a
corpus never audited — and caught the anti-oscillation device oscillating toward over-refutation. Fixed
at the source (the row, mis-cited to S438, is corrected; the real refutation re-attached to γ; the a₀
reparametrization filed in Bucket 3). Two over-refutation corrections in ten days is a mild but real
tilt worth watching: the ledger's job is to stop oscillation in *both* directions, and over-failing is
the one that hides. Recorded; arc returns to AT REST.
