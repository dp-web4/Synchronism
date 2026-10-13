# One "fired" rule, and a scoped "no surviving differing prediction" sentence

*Site maintainer, 2026-10-13. Source: simulated-visitor graduate and researcher personas (site visitor log
2026-10-13), checked against this ledger and the site source. No bucket moves; count 6; Bucket 0 = 0.*

## 1. Two kills cleared by the same margin carry opposite verdicts

- **TEST-09** registered a point-estimate rule: slope deviation > 0.3. Executed deviation 0.41 ± 0.124, i.e. the
  point estimate crosses the threshold by ~0.9σ. Verdict: FIRED, counted. (This ledger already records that V_flat
  alone gives P(dev ≤ 0.3) = 0.11, about 1.2σ, and that the kill is convention-dependent.)
- **TEST-04a** registered a significance rule: fσ₈(z = 0.51) > 0.46 at > 3σ. LRG1 data clear 0.46 by ~0.9σ
  (2026-10-10 correction). Verdict: not fired, underpowered as registered, not counted.

Both verdicts follow their registration texts, so neither is wrong. But the headline "6" mixes two standards of
"fired". Under any uniform significance rule (for example "threshold excluded at 2σ"), TEST-09 does not fire under
V_flat. Under a uniform point-estimate rule, TEST-04a *fires*: the LRG1 point estimate 0.516 is above the 0.46
kill threshold. So a uniform significance rule would cost TEST-09 (and, by the shared inequality, possibly TEST-10),
and a uniform point-estimate rule would add TEST-04a, a seventh refutation on a mechanism already withdrawn
(2026-08-11). Either uniform rule moves the count, in opposite directions. That is why the rule must be printed.

**Recommendation (gates on dp):** print the rule beside the count: "counted if the registered criterion, as
worded, is met". Separately, mark TEST-09 "fired on the point estimate; ~1 σ beyond threshold under V_flat". The
site now carries this explanation on the Tier 1 legend (2026-10-13). This joins the 10-12 counting-rule question:
both ask whether "6" is one rule applied uniformly. At present it is "registration text, row by row".

## 2. State the scoped negative the ledger already implies

The 09-25 and 09-27 blocks establish that no registered test has a branch that can move Bucket 0, because every
surviving construction is nested in MOND (galaxies) or ΛCDM (cosmology), and every differing prediction lives in the
excluded density-keyed branch. A researcher persona asked for this as one sentence. The site now says:
"No prediction that differs from MOND or ΛCDM is alive, within the registered tests and the readings run so far",
with an explicit scope covering unbuilt hidden-communication substrates, door #3, and the never-run L3 on SPARC.
**Recommendation:** put the same sentence, with the same scope clause, under "The disciplined headline". It is the
true third clause of that headline ("zero confirmed; a short list of untested bets; no surviving *differing*
prediction among the tested readings"). Without the scope clause it would be an over-refutation. With it, it is
the ledger's own content in one line.

## 3. Record corrections (no verdict changes)

- **TEST-25 channel.** The site card credited the q = 2γ high-acceleration return as Cassini's channel. This
  ledger's 09-17 block already says Q₂ is set near ~7,000 AU (g ~ a₀), and `simulations/sparc_cassini_q2.py`
  computes Q₂ by DHF 2024 Eq. 12 (phantom-density integral with g_ext). McGaugh's ν, whose far tail is
  exponentially Newtonian, fails the same pipeline at +15.9 to +20.9σ. Site corrected. q = 2γ is a label of the same
  γ, not the mechanism. No ledger text to change, except to read "slow simple-μ return" anywhere in this archive
  as "transition-region shape".
- **B1 run sizes.** 400,000 trials per setting pair (scripts 02, 03, 05) and 200,000 (04). The binomial error on S
  is ≈ ±0.005 (the scripts do not report it). The nonlocal-grid sweep reads 2.000–2.004 over ten couplings, which is
  2 within error. The global-clock 2.67 arm signals (setting-dependence 0.43), so it is not a Bell test.
- **A persona inversion, caught.** The researcher persona read the A2ACW "specificity 0/6" as "the filter passed
  every reparametrization". It means the opposite: the audit flagged all six canonical *discoveries* as prior art.
  The maintainer copied the inversion onto the landing page and reverted it within the session. Logged here because
  it belongs in the denominator of persona-supplied corrections (over-refutation share, 2026-09-17).

## 4. Explorer topic seeded

`synchronism-site/explorer/topics/environment-pipeline-rekeyed-on-external-acceleration.md`. Re-key TEST-08's
environment pipeline on g_ext. The persona framed this as a MOND EFE test. It is also the cheapest split between two
readings of C_g: field-equation (EFE inherited) and algebraic internal-only (EFE = 0), and the archive has never
stated which acceleration keys C_g. Both readings are nested in MOND, so the test cannot move Bucket 0.
