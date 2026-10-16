# Proposal: the ledger's tables lag its own correction blocks; TEST-09 under three rules; the direct-tail floor

*Site maintainer, 2026-10-16. Sources: simulated-visitor pass 2026-10-16 (four personas), explorer finding 2026-10-15.
No bucket moves; count 6; Bucket 0 = 0.*

## 1. The pattern the visitor named applies to PREDICTIONS.md too

The visitor's cross-persona synthesis found that on the site "the bottoms of pages are more honest than the tops."
Badges, fine print and revision notes had caught up, but lead sentences, headings and sibling pages had not. Four of the
top six fixes were this one defect:

- /key-claims Claim 1 still said the reframe "generates different [predictions]" above sub-cards all badged
  Reparametrization.
- /dark-matter's body said a₀ "emerges from cosmology" above its own "1983 coincidence" box.
- The a₀ card's first line said a₀ "follows from the coherence function" under a Dimensional-Analysis badge.
- /cdm-discrimination's net verdict said "consistent with CDM" above its own suspension paragraph.

All four are fixed today, and lint rules now stop each wording from coming back.

PREDICTIONS.md has the same structure, inverted. Corrections arrive as dated 📌 blocks at the top of Bucket 1. The
bucket **tables** further down are where a reader looks for the current verdict, and those rows are not edited:

- **Bucket 2, row "C(ρ) ⇒ MOND (γ free)":** "γ=2-pinned refuted at ΔBIC = +184 on SPARC". Two later blocks qualify it:
  - 09-08: the +184 is C_g at pinned γ, not C_ρ, so the row label names the wrong function.
  - 09-29: +184 is N_eff-dependent. Galaxy-level ΔBIC is ≈ +11 (1.6–2.2σ), "at the threshold, not decisive".
- **Bucket 2, row "Locality is the discriminating axis":** "C(ρ) γ=2 ΔBIC=+184", with the same two issues.
- **Table rows overall:** 0 rows mention galaxy level or N_eff. The 09-29 block is not reflected anywhere in the tables.
- **Footer:** still reads "Last updated: 2026-06-21". About forty dated blocks have been added since.

A reader who goes straight to the table (the ledger's own instruction: "pin the prose to the board") reads the
pre-09-29 verdict.

**Recommendation (gates on the project lead):** when a block corrects a table row, edit the row in place with a short
dated bracket and a pointer to the block. The block stays as the record. Optionally, add a ledger lint like the site's
`maintainer/tools/site_lint.py`, a list of retired row phrasings checked before commit. I have not edited the rows,
because the tables are the ledger's canonical surface and the edit style is yours to set. This proposal is the request.

## 2. TEST-09 under three rules (number supplied for the counting-rule question)

A researcher persona noted that TEST-09's registered rule (|Δn| > 0.3) has no error term. From
`synchronism-site/maintainer/scripts/which_C_carries_the_floor_output.txt`:

- Observed: n_obs = 3.75 ± 0.10.
- Predicted, three candidate ceilings: 3.35 ± 0.07, 3.46 ± 0.06 and 3.49 ± 0.06.

As σ-distances, with errors combined in quadrature: **3.4σ, 2.6σ, 2.2σ**. Both slopes come from the same galaxies, so
quadrature is approximate. Positively correlated errors would make these larger.

| Rule | 1/Ω_m | (Ω_m−Ω_b)/Ω_b | Ω_m/Ω_b | TEST-09 in the count? |
|---|---|---|---|---|
| Registered: point estimate > 0.3 | fires | no (0.30) | no | convention-dependent |
| Threshold excluded at 2σ (10-13) | no (~0.9σ past) | no | no | out |
| Zero tolerance, 2σ | fires | fires | fires | in, convention-free |
| Zero tolerance, 3σ | fires | no | no | convention-dependent |

So the TEST-09 entry in "6" moves under every candidate rule, in both directions. A zero-tolerance 2σ rule is the only
one that makes it convention-free. But that rule also applies to rows where a tolerance was registered on purpose. This
adds a row to the 10-12/10-13 counting-rule question and does not settle it. The site card now shows the σ-distances
beside the registered verdict.

## 3. Direct-tail floor (explorer 2026-10-15, propagated to the site today)

The compander's power-law return to Newton gives an EFE-free anomalous pull at the planets. On EPM2011 supplementary
perihelion precessions it excludes γ < 0.87 at 3σ:

- Saturn binds. γ = 0.489 predicts −9,576 mas/cy against −0.32 ± 0.47.
- INPOP15a's conservative interval gives 0.78.

That empties the point-level SPARC intersection without Q₂ (ΔBIC_profile +59). It does not empty the galaxy-level
intersection: the 95% interval there is [0.42, 0.95], which leaves an EFE-free corner at γ ≈ 0.87–0.95 for about 1.2σ.
For "EFE = 0" readings, the next adjudicator is the SPARC EFE detection (Chae et al. 2020, ApJ 904, 51, >4σ), not
Cassini. This is inherited at γ = ½ (simple μ) and does not change the count.

Ledger note: TEST-25's reading × test matrix gains a column. The "EFE = 0" escape now has a named Solar-System floor
(0.87) and a named untested galaxy-side test.

## For the project lead
1. Should table rows be edited in place when a block corrects them (item 1)? If yes, the two +184 rows are the first
   instances.
2. TEST-09's σ-distance table (item 2) joins the counting-rule decision.
3. Carried: the "fired" rule, the γ statement for FUNDAMENTALS/STATUS, the post-cutoff A2ACW arm, contacting the RG
   authors.
