# Proposal: "Citable" should mean "passed a pre-send check". Today's three checks: LLR kernel and bound edition, Cassini calibration, Ciocan's redshift

*Site maintainer, 2026-11-28. From the site visitor log 2026-11-28 (graduate and researcher personas, both simulated, LLM) and
the site explorer's 11-27 pre-send check of the RG note. No bucket moves. Count 6 (range 4–8). Bucket 0 = 0.*

## What happened
The explorer's 11-27 pre-send check of the note to the Refracted Gravity authors found that **every number held, and the note
still wasn't ready to send**. It asked the authors something their papers already answer. It used an LLR bound edition 9.2×
looser than the current one. And it did not name the escape a referee would raise first (kernel shape). The next morning a
researcher persona read the five "results about other people's theories" on /for-researchers and found the same kinds of gap in
the items the site already calls citable:

1. **LLR Ġ/G no-go for ambient-density-keyed G** (the "if you read one thing" item). Its 1 pc rate was "set by α Cen and
   Barnard's Star", at 1.34 and 1.83 pc. That is a Gaussian-tail number: under a 1 pc top-hat kernel neither star is inside.
   **Executed today** (`synchronism-site/maintainer/scripts/llr_density_keyed_kernel_and_bound_edition.py`, PREREG `13639ec`
   before the script; it reuses the explorer's 11-18/11-27 inputs; 3 of 4 predictions held, and P2 failed on the high side):

   | ℓ | Gaussian (HM18) | top-hat, gas only (HM18) | top-hat (B21) |
   |---|---|---|---|
   | 1 pc | 9.6×10⁻⁹ | 3.6×10⁻⁷ – 1.2×10⁻⁵ | 3.9×10⁻⁸ – 1.3×10⁻⁶ |
   | 3 pc | 1.3×10⁻⁷ | 4×10⁻⁷ – 1.3×10⁻⁵ | 4.4×10⁻⁸ – 1.4×10⁻⁶ |
   | 10 pc | 5.2×10⁻⁶ | 1.2×10⁻⁴ – 4.1×10⁻³ | 1.3×10⁻⁵ – 4.4×10⁻⁴ |

   The site's "at every ℓ at least about four orders below O(1)" was true for a smooth kernel only. **For a sharp-edged kernel at
   10 pc it is ~2.4 orders (2018 bound) or ~3.4 (Biskupek 2021).** The verdict stands, because no cell comes within 100× of O(1). The
   10 pc gas rate is a cloud-edge model and is the weakest input. Star crossings during the 48 yr span are visible steps (8–23 %
   chance at 10 pc).
2. **Cassini squeezes every compander.** With Υ free the in-house pipeline reproduces Desmond, Hees & Famaey 2024's mitigated
   1.9σ (explorer 11-07: 1.8σ). It does not reproduce their threshold: they find Cassini demands δ ≳ 2.5, while in-house the δ
   family first passes at δ ≈ 4–6 depending on the prior width. The compander's σ values are read off a pipeline that is harsher
   than the published one on the published family. Queued to the explorer
   (`reproduce-dhf24-cassini-delta-threshold-before-the-compander-item-is-citable.md`).
3. **a₀(z), Ciocan et al. 2026.** Checked against the paper today (arXiv:2604.22613 §3.1, Eq. 2). The 2.38 (+0.12/−0.10, 95 % CI)
   is **one fit to the whole sample, 0.33 < z < 1.44**, labelled "z ∼ 1", and no median z is given. The paper's own line, a₀ = 1.00 +
   1.59z (Eq. 4, both fitted, 95 % CI), reaches 2.38 at z = 0.87, where E = 1.66, not E(1) = 1.79. At z = 0.87 the level rows
   move as follows: the framework's own cH₀/2π goes from 9.8σ to 12.5σ, Ciocan's intercept from 9.4σ to 11.9σ, McGaugh from +0.5σ
   to +0.9σ, and MIGHTEE-HI from −2.3σ to −1.4σ. **No verdict moves** (the level comparison was already ruled non-discriminating;
   the slope statistic doesn't use E at a single z). Script: `synchronism-site/maintainer/scripts/ciocan_a0_effective_z_rows.py`.
4. **The grid's p⁴ term.** The convention was unstated: cell = 1/M_Pl with the unreduced M_Pl = 1.22×10¹⁹ GeV (from
   `grid_cubic_p4_sme_and_thresholds.py`). Because η ∝ (h·M_Pl)², a reduced-Planck cell multiplies η by 8π, while the verdicts sit at
   factor-of-two precision. **With the spacing free, every bound is a lower limit on the spacing, not a test.** That belongs beside
   "commit a fermion action" (11-27): the commitment that would make B7 a bet is a stencil **and** a spacing.

## Proposal
1. **Adopt a pre-send checklist as the meaning of "citable"** on /for-researchers and in any outreach. Before an item is called
   transferable it needs four things:
   - (a) the newest bound edition, named;
   - (b) every free convention stated, with the factor it moves the verdict by (kernel, Planck mass, CI level);
   - (c) the published instrument reproduced on the published case (a positive control at the real configuration);
   - (d) the escape a referee would raise first, named and quantified.
   The RG note now has all four. LLR has (a), (b) and (d) as of today. Cassini lacks (c). a₀(z) is a proposal, not a result.
2. **Record the checklist state per item on /for-researchers** (one line each). That costs nothing and tells an outside reader
   which items have been treated like the RG note.
3. **For dp:** the RG note v2 is send-ready
   (`draft_note_to_refracted_gravity_authors_llr_bounds_the_smoothing_length_20261127.md`, "Pre-send check" section). DR4 is
   on 2 December. If no reading key is adopted by then, the maintainer records the reading per named variant, as stated 11-27.

## Why this is a research-direction point, not a site fix
The program's external value now rests on a handful of framework-independent constraints (SPINE: "unproven for lack of
instruments"; the site: "no dataset on the horizon can confirm a novel prediction"). Those constraints are the only output an
outside physicist can test. Their failure mode is not wrong arithmetic: ~30 hand checks by two personas today all held. Their
failure mode is unstated conventions and uncalibrated instruments. The program has built lints for wording and a findings-lint for
registration. It has no check for "would a referee in this subfield accept this as stated". The RG note shows what that check
finds and what it costs: one explorer session per item.
