# Verify Peres–Mermin (contextuality) + triage TEST-03 "manufactured kill" + sharpen the over-refutation meta-claim (2026-07-08)

**Status:** `[ACTIVE-MRH]` — checklist tripped by explorer commit `24049d42` (Peres–Mermin suite 06 + two
proposals). Three things handled: (1) QA the contextuality artifact; (2) triage the "TEST-03 kill is
manufactured / 6-over-refutation pattern" proposal; (3) sharpen its meta-claim before it becomes a new
unscrutinized slogan. **Verdicts: (1) Peres–Mermin SOUND (re-executed, textbook KS state-independent
contextuality) — added a distinct contextuality leg to B1. (2) TEST-03 kill is real but ENTIRELY SITE-SIDE —
the ledger never carried it (grep-confirmed clean); load-bearing primary verified; routed to maintainer/dp.
(3) "6/6 over-refutations, 0 over-claims" is directionally right but OVER-STATED — the honest form is an
audit-posture ASYMMETRY, and it lives on the site/whitepaper compilation layer, not the ledger. Bucket 0
unchanged (0); arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## (1) Peres–Mermin contextuality — SOUND (re-executed)

`06_peres_mermin_contextuality.py` reproduces identically: 6 contexts (3 rows, 3 cols) with products
+I,+I,+I,+I,+I,**−I**; exhaustive search over 2⁹ assignments → **0/512 satisfy all six, best 5/6**; parity
argument confirms it (each of the 9 observables sits in exactly 2 contexts, so ∏ context-products = (∏
values)² = +1, but the constraints require −1 → contradiction). The QM structure is **derived from the Pauli
algebra** (commutation within each context checked), not asserted. This is textbook state-independent
Kochen–Specker contextuality.

**Meaning (the CRT-scanning proposal, confirmed by execution):** the framework's CRT temporal-scanning
picture — a fast deterministic cycle through *definite* per-location values, sampled by the observer — is a
**non-contextual hidden-variable model**, exactly the class KS excludes for Hilbert dim ≥ 3. So "CRT vs KS"
was never an open problem pending formalization; it is theorem-determined. The remaining question is only
whether a *contextual* variant exists — and contextualizing the scan means the sampled value depends on the
global measurement context, not just sampling timing, which **sacrifices the very content that made the
analogy compelling** ("just timing"). Added a precise B1 annotation: the non-classicality gap now has **two
distinct theorem-level legs — Bell/Tsirelson (CHSH-05) and KS-contextuality (this)**. Kept them distinct
(Bell ⊂ contextuality, but they are different theorems — no over-unification, per the 07-07a lesson). Bucket
0 untouched; this strengthens a refutation, opens nothing.

## (2) TEST-03 "manufactured kill" — real, but site-side; ledger clean

The proposal alleges the site's TEST-03 verdict ("Failed — R²=0.14 < pre-registered kill 0.20") is
manufactured via five defects. **My-lane check first: does the ledger carry it?** `grep -iE
"TEST-03|TEST-08|R²=0.14|pre-registered kill" PREDICTIONS.md STATUS.md` → **empty.** The ledger does not
carry the TEST-03 kill — same situation as last session's B3 "0.64 also rejected" fabrication (site drift,
clean ledger). So there is **no ledger fix to make** and, importantly, **no ledger softening happening** —
the both-directions worry (a run of un-refutations over-softening PREDICTIONS.md) does not bite, because the
over-refutation lives on the site, which I don't own.

**Verified the load-bearing primary anyway** (to gauge reliability): `EXPERIMENTAL_TEST_CATALOG.md` TEST-08
(lines 126–138) registers **Expected: >20%** and **Falsification: r² < 0.09**. So 0.20 was the *success*
target, not a kill bar; R²=0.138 **passes** the registered threshold (0.138 > 0.09); and TEST-08 is about
*environment* (cluster/field/void), whereas 0.138 was S377's Hubble-type *morphology* term — a different
variable. The proposal's claim #1 is exactly right. (Claims 2–5 — catalog postdates measurement, N=14,585
belongs to S591's BTFR, p/N/R² mutually impossible, "Session 616" misread — are internally coherent and
consistent with the verified #1; I did not separately re-derive each, since the disposition is site-lane
regardless.)

**Routing:** site "manufactured kill on ~12 pages" = **maintainer P0** (site is do-not-modify for me).
Running TEST-08 as registered (SPARC×environment catalogs, one afternoon) and adopting the standing
kill-criterion rule = **dp decisions** (the proposal itself files them as "Decision for dp"). No CBP
physics-ledger action.

## (3) Sharpening the "6/6 over-refutations, 0 over-claims" meta-claim (the real catch)

The proposal's transferable thesis — "on an honesty-branded corpus, self-refutations are the least-
scrutinized claims, so manufactured failures pass; compilation layers drift toward the brand" — is a real and
useful observation. But the tally **"6 over-refutations, 0 over-claims" is itself mildly over-claimed**, and
saying so matters so it doesn't become the next unscrutinized slogan:

- **There IS a documented over-claim, on the strength side:** 07-07a, the publisher's "two-scale locality
  no-go" over-unification (bundling CHSH-05 + door-#1 as "one locality deficit"). That was an over-*claim* of
  a strength, which I flagged. So "0 over-claims" is false unless you scope it to *refutation-provenance
  walks only*.
- **The honest form is an ASYMMETRY, not a one-sided count:** an audit posture that hunts overclaims
  under-scrutinizes failure-claims → systematic **over-refutation on the failure side**, *and* over-claiming
  still occurs wherever the audit isn't looking — which on a failure-hunting honesty-corpus is the strength/
  packaging layer (the two-scale synthesis). Over-claims don't vanish; they **migrate to the unwatched
  surface.** This is a strictly sharper version of the proposal's own "drift toward the brand."
- **The drift lives on the SITE/whitepaper compilation layer, not the ledger.** The proposal's 6 exemplars
  are site-side (TEST-03, S63 "0.64" fab, CDM inversion, TEST-04a, wide-binary, ΔBIC). The ledger's *own*
  defects were a different, subtler class — citation-rot (dim-4 over-claim of a no-go; a₀ mis-bucket; B3
  salience-mislabel) — which the ledger-edge walks fixed at source. **PREDICTIONS.md is the comparatively
  clean control surface; brand-drift over-refutation accumulates downstream on the compilation layers.** Good
  news for the anti-oscillation device; the action item is on the surfaces that inherit from it.
- **Discipline consequence (mine):** the existence of a strength-side over-claim means I must NOT adopt "the
  corpus only ever over-fails" — that belief would license reflexive softening, the exact failure mode. The
  both-directions check stays mandatory; here it resolved cleanly only because the ledger was already clean.

## Disposition

- **PREDICTIONS.md B1** — added the contextuality (KS) leg (verified strengthening; Bucket 0 unchanged).
- **No other ledger change** — TEST-03 is not in the ledger; nothing to un-refute in my lane.
- **Maintainer P0:** site TEST-03 kill (~12 pages) + adopt the kill-criterion-walk rule on the site.
- **dp decisions:** run TEST-08 as registered; adopt the kill-criterion rule as a standing ledger norm.
- **Bucket 0 unchanged (0); arc AT REST.**

## So what

Two quantum-foundations legs now close by *theorem+execution* on the framework's own construction — Bell
(CHSH-05) and KS-contextuality (Peres–Mermin) — which is the strongest form the non-classicality cage takes:
not "we couldn't find a violation" but "no real-valued non-contextual local model can, and here it is on your
own square." Meanwhile the TEST-03 walk is genuine but site-side, and its most valuable output — the
over-refutation pattern — is best stated as an audit-posture *asymmetry* (over-claims migrate to the unwatched
surface) rather than a "0 over-claims" tally that my own 07-07a two-scale catch already contradicts. The
ledger stays clean and honest in both directions; the drift to fix is downstream. Bucket 0 still 0; arc AT REST.
