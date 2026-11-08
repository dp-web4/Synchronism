# Literature statements fail where computations hold; TEST-10's "≲ 14" was a sample maximum

**Date:** 2026-11-08 (site maintainer, from simulated-visitor graduate and researcher personas and site explorer 2026-11-07)
**Status:** proposal. No bucket moves, count 6, range 4–8, Bucket 0 = 0. Rulings gate on the project lead (section 4).

## 1. The a₀ prior-art credit was wrong, and so was part of the correction to it

The ledger's Bucket 3 row reads "Milgrom's a₀ ≈ cH₀/6 coincidence (1983) … (1/2π = 0.159 vs 1/6 = 0.167, ~5%; cf. McCulloch / Verlinde /
Smolin)". The site went further: three frameworks "derive the same relation with the same geometric factor", and McCulloch "derived
a₀ = cH₀/(2π)". A researcher persona read the papers. I re-read them through a separate fetch before editing:

| Source | What it prints | × cH₀ |
|---|---|---|
| Milgrom 2008 (arXiv:0801.3133 §VI), of his 1983a | "observed that 2πa₀≈cH₀" | 0.16 |
| Sanders & McGaugh 2002 (astro-ph/0204521) | abstract "~cH₀/6"; §2.1 Milgrom noticed a₀≈cH₀ "to within a factor of 5 or 6" | 0.17 |
| Famaey & McGaugh 2012 (arXiv:1112.3960 §4.1) | "a₀≈cH₀/2π", no derivation, no attribution | 0.16 |
| McCulloch 2007 (astro-ph/0612599 Eq. 12) | 0.99×cH, "about six times larger" than MOND's | 0.99 |
| Verlinde 2016 (arXiv:1611.02269 Eq. 1.7) | a_M = a₀/6 in d = 4, a₀ ≡ cH₀ | 0.17 |
| Smolin 2017 (arXiv:1704.00780 Eq. 3) | a₀ ≈ a_Λ/8.3, a_Λ = c²√Λ | ~0.17 |

- **"Same geometric factor" is false for all three**, and for McCulloch by a factor of 6. The shared thing is the dimensional origin.
- **The persona's own correction was also partly wrong.** It said Verlinde's is "the only /6 I located". Sanders & McGaugh print /6 in
  their abstract, and Milgrom 2008 attributes **the 2π form** to Milgrom 1983a. So cH₀/2π, which the archive treats as its version, is the form
  the MOND literature itself prints. That strengthens Bucket 3.
- **Milgrom 1983 itself was not read** (not on arXiv; ADS refused the fetch). "Milgrom wrote cH₀/6" has no locator anywhere in this archive.
  Read that row as "Milgrom 1983's a₀ ~ cH₀ coincidence".

Edit to the row gates on the project lead (10-16 (i)). The site carries the table at `/mond-unification#a0-prefactors`.

**Also from the re-read:** the a₀ text is the fifth compiled-vs-primary drift since July (DESI's QSO row, Hulse–Taylor 1.0013,
Ciocan's CI level, the chemistry table). /wide-binaries had Saad & Ting at γ = 1.12 ± 0.25 and Chae at "4.9σ"; the current abstracts give
γ = 1.00 (+0.24/−0.19) and no σ, and Chae's says the strong anomaly is dominated by two pairs (γ ≈ 1.2 without them).

## 2. TEST-10's "B_max ≲ 14" is not a bound

This is the only row that refutes the framework's own content on outside data. Its badge read "B_max ≲ 14 excluded at every convention". That
14 is the rounded face-value B of one disc in 153, with no significance. Pre-registered (site `0be9704`, before running) on the explorer's 10-29
Statistic B, identity control exact (14.2 / 9.5 / 35.1):

| | P | Result |
|---|---|---|
| Identity | P1 | held |
| Largest cap with a disc at χ² ≥ 9, gas on, is 7–10 | P2 | **failed: 14**, through UGC 06667 |
| Binding disc is UGC 06667 or NGC 2915 | P3 | held |
| Family-wise bound (χ² ≥ 18.4) is 6–8 | P4 | **failed: 9**, same disc |
| Doubling the gas prior moves the 3σ cap by ≤ 1.5 | P5 | held (no move) |

**Why P2 and P4 failed, post-hoc:** UGC 06667's excess at high caps sits in its inner points, where SPARC lists V_gas = 0.0. The gas nuisance cannot act
on a zero, which is what the 10-29 table's unexplained "gas × ∞" was (UGC 00731 is the same). Restricted to the **outer half of each curve**, closer to
TEST-10's own statistic (post-hoc, `maintainer/scripts/test10_bmax_3sigma_diag_posthoc.py`):

| Cap | discs ≥ 3σ | discs clearing χ² ≥ 18.4 (153-disc correction) |
|---|---|---|
| 3.175 (1/Ω_m, the framework as written) | 28 | 10 |
| 5.389 | 4 | 1 |
| 6.39 (Ω_m/Ω_b) | 2 (NGC 2915, NGC 3741) | **0** |
| exact edges | **≈ 7.5** (NGC 2915) | **≈ 5.7** (NGC 2915) |

**Reading.**
- The framework as written is not close. The 3.17 cap fails overwhelmingly.
- "TEST-10 fires under every candidate cap" stays true under the rule as registered (no significance term), so the count and range do not move.
- But at the most permissive cap the rotation-curve case is **per-disc, not trials-corrected**. The max-over-discs test is conservative: a population statistic could
  be much stronger and has not been run. Lensing is the independent leg, and it is census-conditional.
- This qualifies the 10-30 (iii) "3 discs stay ≥ 3σ": two on the outer curve, the third through gas-free inner points.
- Limits: no asymmetric drift, beam smearing or non-circular motions; χ² treated as χ²-distributed despite one-sided residuals.

## 3. Explorer 11-07 is now the site's statement of TEST-25

Every in-house Cassini σ is now labelled M/L-pinned, and one status sentence renders on three pages. The robust form is the Υ_d each family needs in order to pass: MOND δ = 6 ≈ 0.54
(no bulges; ~1.8σ on SPARC), compander ≈ 1.2–1.4 (36–42% of discs over-maximal; 2.9σ at γ = 3); γ ≤ 1 and δ = 1 pass at no Υ_d ≤ 1.5. Headlines lead with the published 1.9–15σ. **Open flip:** AQUAL.
Desmond+ say it raises the tension, and if MOND δ = 6 goes above 2σ the route closes for both families.

The Desmond abstract also says the quadrupole "implies, to high precision, no deviation from Newtonian gravity in nearby wide binaries", and Park+ 2026 bound the MOND boost at the Sun at 2% (95%). This was
not on /wide-binaries and bears directly on reading DR4 (2 December 2026). It is there now.

## 4. Frame, for the project lead

**(a) What is a correction note for, if no test can move the count?** A researcher persona asked this directly, and it restates 10-04 (stopping) with a new observation from a technical-writer persona: every
fix added words, so the loop slowly degrades the pages it repairs. Today I folded dated corrections into per-page revision logs and removed provenance-only notes (visible words on 9 long pages −1,213, despite
additions). That is reversible and is not the same as a retirement rule. **Recommended:** adopt the rule that a correction is folded into the page's log in the session after it ships, and that a page has a
word budget.

**(b) The most useful number this program could publish is its own error rate.** Same persona: sample 50 verdict sentences at random, have them re-derived independently, publish the fraction that hold, split
by literature vs computed. Today's evidence suggests literature statements score worse. It is internally runnable, its outcome is unknown, and an outsider could cite it. It is a claim about method and cannot
move Bucket 0. Explorer topic seeded. **Whether it replaces further physics corrections as the maintainer's default gates on the project lead.**

**(c) Still open:** the post-cutoff arm (routed since 09-10), the count rule, and adoption of the DR4 reading key, which needs a decision before 2 December.
