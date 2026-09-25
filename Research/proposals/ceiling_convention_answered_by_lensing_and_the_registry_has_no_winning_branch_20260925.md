# The ceiling-convention question has a data answer, and the registry has no branch that moves Bucket 0

**Date:** 2026-09-25 · **From:** site maintainer (back-annotating visitor log 2026-09-25 and explorer 2026-09-02 §5)
**Gates on dp:** items 1 and 3. Item 2 is a record correction and is inscribed in PREDICTIONS.md today.
**Count:** unchanged (6). **Bucket 0:** 0.

---

## 1. The "convention-dependent" asterisk on TEST-09/TEST-10 is answered by weak lensing

Since 2026-07-29 (TEST-10) and 2026-09-18 (TEST-09), the site and the ledger have said that both galaxy
kills fire at B_max = 1/Ω_m but not at every candidate ceiling, so "whether they stay in the count of 6
gates on dp". Three candidate ceilings are in play: 1/Ω_m = 3.175, (Ω_m−Ω_b)/Ω_b = 5.39, Ω_m/Ω_b = 6.39.

Galaxy–galaxy weak lensing reaches accelerations rotation curves cannot. The KiDS-1000 isolated-lens RAR
(Brouwer et al. 2021, A&A 650, A113) follows the extrapolated MOND branch from g_bar = 5×10⁻¹² down to
≈10⁻¹⁵ m/s². The site explorer set the ceiling against it on **2026-09-02** (finding
`synchronism-site/explorer/findings/the-last-escape-is-mond-induced-and-the-column-was-chi2-not-rho-c.md` §5,
P0 for the maintainer). That P0 sat on one site card for 23 days, while every public scoreboard carried
"convention-dependent, pending dp". A visitor researcher persona re-derived the argument from citations
today.

Today's recomputation (`synchronism-site/maintainer/scripts/lensing_ceiling_every_convention.py` +
`_output.txt`; predictions written into the header before running; not blind, since the f = 1 and f = 5
values at 10⁻¹⁴ were hand-estimated while reading the log) changes one thing about the explorer's version.
Hidden gas lowers the boost the data *require* by the full factor f, B_req = ν_obs/f, and not by √f
(the √f in 09-02 §5 is the MOND re-prediction view, which does not apply to a ceiling).

| g_bar (m/s²) | B_req, f=1 | B_req, f=5, −0.3 dex | ÷ 3.175 | ÷ 6.389 | ÷ 20.3 (1/Ω_b) |
|---|---|---|---|---|---|
| 10⁻¹³ | 35 | 3.5 | 1.11 | 0.55 | 0.17 |
| 10⁻¹⁴ | 110 | 11.0 | 3.47 | **1.73** | 0.54 |
| 10⁻¹⁵ | 347 | 34.8 | 10.95 | 5.44 | **1.71** |

f = 5 puts every cosmic baryon inside the lensing radius (abundance matching gives M★/M_h ≈ 0.02–0.03 for
these lenses). −0.3 dex is a downward allowance on g_obs relative to the MOND branch. All three registered
predictions held: every Ω_m-based cap is exceeded by ≥ 5× at g_bar ≤ 10⁻¹³ with no allowance, and 6.39 is still
exceeded at 10⁻¹⁴ with both allowances.

**Reading.** The *registered* TEST-09 kill is convention-dependent (a statement about SPARC's slope). The
*ceiling it tests* is not: every Ω_m-based cap is excluded by lensing with generous allowances. The kill
rests on the two lowest bins, 0.4–4 Mpc from the lens, where the isolation criterion and neighbouring haloes
are the live caveats. The Brouwer data were not re-read here (no network in the sandbox). "Tracks the MOND
branch" is the paper's own statement, relayed by the explorer.

**For dp:** retire "convention-dependent" as the reason the count is under review. The remaining reason is
independence (#1 and #2 are one inequality; ≤ 5 independent), which the site's new one-table ledger
(`/honest-assessment#refutation-ledger`) now shows directly.

### 1a. The frame point underneath: every cap number is CDM-conditional

All three conventions are ratios built from Ω_m, and Ω_m counts cold dark matter. In the no-CDM reading
(Ω_m → Ω_b), the site's own convention 1/Ω_m becomes 1/Ω_b = 20.3. The other two degenerate to 0 and 1.
That is the only cap a no-DM framework could motivate, and it survives lensing only at 10⁻¹⁴ with every
allowance granted. It fails at 10⁻¹⁵ (row 3). In the CDM reading, the explorer's 09-22 result applies:
a real halo plus the boost overshoots by +0.28 dex. So the convention question was never independent of the
dark-matter fork (proposal `de_sector_cannot_carry_the_dark_matter_cosmology_needs_cdm_20260922.md`). It is
the same fork seen from the galaxy side. **Recommendation:** fold the convention question into the DM-fork
decision instead of keeping it as a separate gate.

## 2. Record corrections (inscribed today, no bucket moves)

**(a) TEST-09's "deep limit n → 2 (verified numerically: 2.01)" is a fixed-radius statement.** A saturated
boost is Keplerian, and a Keplerian law has no radius-independent BTFR slope. For measurement radii R ∝ M^α
it gives n = 2/(1−α): n = 2 at a common radius, ≈ 3 for real disc radii (α ≈ 0.3–0.4), and **n = 4 at a fixed
g_bar** (α = ½). What is radius-independent is the declining curve shape. The verdict does not move: the
executed 3.35 used each galaxy's real radii. (Visitor graduate-physics persona, 2026-09-25. The n = 4 case
is the maintainer's extension.)

**(b) B4's "homogeneous case confirmed (r = 0.994)" is a correlation read as proportionality.**
`Research/Compatibility_Synthon_Experiment.md` Experiment A: across compatibility 0.2 → 1.0, p_crit goes
0.0320 → 0.0185. An affine fit gives p_crit ≈ 0.0151 + 0.0034/C (residuals ≤ 0.0014). The intercept is **82%**
of p_crit at C = 1. The stated law p_crit ∝ 1/⟨C⟩ predicts 0.0925 at C = 0.2, and 0.0320 was observed
(2.9× off). The source itself says "not a perfect 5× inverse (predicted)". The ledger then carried "confirmed".
What the data support is monotone decrease of p_crit with compatibility, affine in 1/C with a dominant
intercept. The proportionality is refuted on its own five points. B4 stays in Bucket 1 (the heterogeneous
law is still untested), and its odds cell now carries the correction.

## 3. The registry has no branch that moves Bucket 0 — and B4 is the one bet that could

Visitor Pass 4 asked: "Is there any registered test with an outcome that would favour Synchronism over
MOND + ΛCDM?" Checked branch by branch:

- **TEST-04a (DR2):** branch A (fσ₈ ≤ 0.46) is booked "Bucket 0 stays 0" by its own registration.
- **TEST-26 (proposed):** thawing/crossing → kill; ΛCDM → tie. Freezing-preferred at γ measurably ≠ ½ *would*
  select the mean-density family over ΛCDM, so the persona's "no winning branch" is slightly too strong.
  But that family is Cardassian-class (Freese & Lewis 2002), so even that branch is prior art: no Bucket-0 move.
- **TEST-02:** self-eliminating-or-tie.
- **TEST-01…24 unrun:** "none still unrun can select the framework as postulated" (site, standing).

So **no test in the TEST-01…26 namespace has a branch that moves Bucket 0.** Discipline 3's loan (bar b,
"the frame must commit to being pushed where current models fail") currently has no repayment instrument
in the registry. That is not a refutation. It is a statement about where effort goes: every registered
outcome is kill, tie or prior-art. The last four weeks of maintainer work have been refutation-side
audits (controls, conventions, corrections), and that work has been productive and honest. It cannot
change the one number the ledger exists to track.

**Where a winning branch still exists:** Bucket 1 bet **B4** (heterogeneous compatibility law). It is
runnable in-house with no instruments, it has a genuine win branch, and after correction (b) it is sharper.
Prior art makes a specific rival prediction: for heterogeneous coupling networks the synchronization
threshold goes as 1/λ_max of the coupling matrix (Restrepo, Ott & Hunt 2005, PRE 71, 036151), not 1/mean.
The two coincide for homogeneous or regular structures and separate for block-diagonal and heavy-tailed
compatibility, which are exactly the "≥ 5 structure types" B4's refutation criterion names. So the
heterogeneous run is a three-way test: 1/⟨C⟩ (B4 as stated), 1/λ_max (prior art), or the affine form
(correction (b)). Seeded to the explorer as
`synchronism-site/explorer/topics/b4-heterogeneous-compatibility-vs-lambda-max-prior-art.md`.

Caveat on novelty: even a 1/⟨C⟩ win would be a coupled-systems result, not a physics prediction, so it
would move Bucket 0 only if dp scopes Bucket 0 to include the applied axis. That scoping question is item 3
for dp. Absent it, the honest line for the site is: *"No registered test has a branch that would confirm a
novel physics prediction."* Say that plainly, next to "0 confirmed", as the second clause of discipline 1.

## For dp

1. Retire the ceiling-convention gate (§1). Fold what remains into the DM-fork decision (§1a).
2. (Inscribed) TEST-09 n → 2 precision; B4 proportionality correction.
3. Does Bucket 0 admit an applied-axis confirmation (B4), or is it physics-only? If physics-only, the site should
   say "no registered test can move Bucket 0". If not, B4 is the program's only live winning branch, and it is
   cheap.
