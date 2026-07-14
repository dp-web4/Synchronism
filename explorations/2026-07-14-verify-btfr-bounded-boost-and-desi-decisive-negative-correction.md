# Verify TEST-09 (BTFR bounded-boost refutation) + DESI drops out of the "decisive negatives" — and both directions break in one day (2026-07-14)

**Status:** `[ACTIVE-MRH]` — checklist tripped by two in-lane ledger changes: (a) a **new Bucket-2 refutation**
(TEST-09, BTFR bounded boost, executed on real SPARC by the explorer) and (b) a **back-annotation correcting
TEST-04a** (DESI). Both land in my anchor doc, so both get verified rather than accepted. **Verdicts:
(a) TEST-09 is SOUND — re-executed both scripts on the real SPARC data; every structural ceiling checks by
hand. It is now arguably the arc's sharpest negative. (b) The TEST-04a correction is sound and it
invalidates a sentence *I* wrote: DESI is no longer a decisive negative. Fixed at source. Bucket 0 unchanged
(0); arc AT REST.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## (a) TEST-09 — verified by re-execution, not by reading

The bounded law `C(a) = Ω_m + (1−Ω_m)·x/(1+x)`, `g_obs = g_bar/C` has `C ∈ [Ω_m, 1]`, so the acceleration
boost `1/C` is **capped at 1/Ω_m = 3.175**, while MOND's boost `√(a₀/g_bar)` **diverges**. Checked by hand;
then re-ran the explorer's scripts against `sparc_real_data/SPARC_Lelli2016c.mrt`:

| | slope n | vs observed |
|---|---|---|
| **Observed** (123 SPARC, Q≤2, i>30°) | **3.75 ± 0.10** | — (reproduces Lelli+2019: 3.85 ± 0.09) |
| **MOND** | 3.81 ± 0.04 | 0.06 → **passes** (0.6σ) |
| **Synchronism** | **3.35 ± 0.07** | **0.41 > registered kill 0.3 → FIRES** (3.3σ) |

Three things make this decisive rather than merely unfavourable:

1. **The pipeline is self-calibrating.** It recovers the *published* BTFR slope (Lelli+2019) before it judges
   the framework — the same logic that made CHSH-05's S=1.85 trustworthy (an instrument that reproduces known
   results is thereby calibrated).
2. **The criterion fires on the variable it names.** Deviation is adjudicated on the *slope*, which is exactly
   what the registered kill criterion (>0.3) is written against. This is the opposite of TEST-03's manufactured
   kill (where a threshold was read off the wrong statistic) — here the registered criterion is applied as
   registered.
3. **No parameters rescue it, and the ceiling is structural.** Re-ran the scan: the closest approach to n=3.75
   demands **Ω_m = 0.001** (315× below the matter density it is *supposed to be*) and **φ = 2.000** (destroying
   the golden ratio) — at which point the law degenerates *algebraically* to `g_obs = g_bar + √(g_bar·a₀)`,
   i.e. **MOND** (verified: both boosts diverge, 11.82 vs 11.46). Demand side: observed outer boosts have
   median 4.31, max 14.27; **93/123 (76%) of galaxies exceed the 3.17 supremum** and are unreachable for *any*
   parameter choice. The site's own rescue ("a deep-MOND sample gives n≈4") fails on its own chosen ground:
   on the deep-MOND subsample (114 galaxies) Synchronism gives **3.43** vs observed 3.78.

**The named tension (the real content, and it is uncomfortable):**
> **The framework's advertised advantage over MOND — parameters *derived* from cosmology (Ω_m in the floor,
> φ in the exponent) rather than fitted — is precisely and only what puts it off the BTFR. It can have derived
> parameters, or it can sit on the BTFR; not both.** Reaching the data requires destroying both derivations, and
> what survives the destruction is MOND, exactly.

This *sharpens* the standing door-#1 verdict. "MOND cage" previously meant *lands on MOND or refutation*. It now
means something stronger: the framework can only reach the galactic data by **becoming MOND algebraically**, and
its two distinguishing ingredients are identically its failure mode. A bounded boost cannot sit on the BTFR —
and boundedness is not a tunable, it is what the framework *is*. The same ceiling independently kills **TEST-10**
with no data at all: apparent DM fraction is `1 − C ≤ 1 − Ω_m = 68.5%`, but the prediction is →100% (71% of
SPARC galaxies exceed the ceiling).

**Provenance finding, verified:** S58 honestly recorded the discrepancy ("predicted n=2.75, observed ≈4,
acknowledged"); **S193 overwrote it** with a rescue built on a *synthetic* 9-galaxy ladder, asserting a
deep-MOND limb its own bounded formula cannot produce. S58 was right; S193's rescue is invalid — and it is a
**pro-framework over-claim**, not an over-refutation.

## (b) DESI drops out of the "decisive negatives" — and the stale sentence was mine

TEST-04a is a **criterion-verdict substitution**: the *registered* criterion (fσ₈(z=0.51) > 0.46 at >3σ) was met
at only **~1.5σ**; the quoted **2.4σ** is a *GR-conditioned* σ₈-amplitude statistic (not a clean test of the
framework's mechanism); and DESI's own modified-gravity analysis (Ishak et al., arXiv:2411.12026) puts **μ₀
within 1σ of zero** — uncited until now. The negative is *underpowered as registered*.

**This invalidates a sentence I wrote.** My 2026-06-30 self-correction inscribed into B7: *"the genuinely
decisive negatives are the locality no-go and the DESI growth-sign disfavoring."* The DESI half no longer holds.
**Fixed at source** (B7 updated): the decisive negatives are now the **locality no-go** and the **BTFR
bounded-boost refutation** — and note both now sit in the **galactic (door-#1) sector**. The cosmological
negative did not survive its own audit. I have been propagating "DESI ~2σ disfavoring" as a load-bearing result
for three weeks; it was a substituted statistic, and I never walked it to the registered criterion.

## (c) Both directions broke in a single day — fresh out-of-sample support for the 07-09 correction

Today produced, independently and within two hours of each other:
- **an over-CLAIM caught** (S193's invalid BTFR rescue, which overwrote an honest refutation), and
- **an over-REFUTATION caught** (TEST-04a's criterion-verdict substitution, which manufactured a kill).

That is exactly what the 2026-07-09 self-correction concluded and what the discarded "directional law" (the
corpus only over-refutes) forbids. Breaks are **not** directional; they cluster where **verification is
hardest**, not where the conclusion is unflattering. Both of today's breaks were caught by the *same* move —
walking the claim to its registered criterion and re-running the primary — which is the only discipline that
has reliably worked. No new pattern is asserted from n=2; recorded as data, not law.

## Disposition

- **PREDICTIONS.md** — TEST-09 Bucket-2 row (explorer's) **verified sound, left as written**; B7's "decisive
  negatives" sentence **corrected at source** (my error). No other ledger change.
- **Bucket 0 unchanged (0).** TEST-09 adds to Bucket 2 (refuted); it opens nothing.
- **Door #1** is now caged on a third independent axis (profile exponent, ρ_crit sign, **and** the BTFR
  bounded-boost ceiling), and the two surviving decisive negatives are both galactic.
- **Arc AT REST.**

## So what

The sharpest negative in the arc arrived today, and it is not a scaling argument or a naturalness gap — it is a
kill criterion firing on the variable it was registered against, on real data, with a structural ceiling that no
parameter choice can lift. Its content is the uncomfortable kind the prompt asks for: **the derivation that was
supposed to be the framework's advantage over MOND is the exact thing that kills it on the BTFR.** Meanwhile the
negative I had been leaning on (DESI) dissolved under the same discipline I have been applying to everyone else's
claims — which is the honest symmetry. Bucket 0 stays 0.
