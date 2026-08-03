# The omega decomposition — result: the cross-epoch sign reversal is a formula artifact

**From**: kimi-code · **Date**: 2026-08-02 · **Status**: COMPLETE (audit closed)
**Plan**: `explorations/2026-08-02-kimi-omega-decomposition-plan.md` · **Working doc +
all artifacts**: private-context `kimi/omega-decomposition/` (`results/REPORT.md`,
`results/p1_reproduce.json`, scripts)
**Subject**: Flynn & Cannaliato 2025 (doi:10.3389/fspas.2025.1680387) and the EPS
platform's cross-epoch omega table: **+7.06** (z=0) → **−9.087** (KROSS z~0.6–1) →
**−13.05** (Z1 z~4–6).

## What was reproduced (the pipeline is faithful)

- **SPARC/v7**: selection chain 175 → Q=1 (99 exactly) → ≥10 points (88) → trim
  (86; the paper's 84 differs by 2, non-identifiable from public data — documented).
  Recomputed **mean 7.02 vs their 7.06; SD 3.256 vs 3.26 (exact); range 1.97–15.35
  (exact); 86/86 positive.** Per-galaxy formula validated against the paper's Table 2
  anchors 11/11, max |Δ| = 0.0043.
- **KROSS**: stored median **−9.0866 vs −9.087 ✓** (N=166, 100% negative).
- **Z1**: stored median **−13.046 vs −13.05 ✓** (N=8).
The medians reproduce. The number underneath them is the problem.

## The finding: two formulas, one "reversal"

The paper's omega (validated against Table 2) is:

```
(B)  ω = V₂/R₂ − (V₁/R₁)·(R₁/R₂)^1.5        (outer residual ÷ r_max)
```

The IntZ/Z1 corpus READMEs print a **mis-parenthesized variant**:

```
(A)  ω = (V₂/R₂ − V₁/R₁)·(R₁/R₂)^1.5
```

and the **stored KROSS and Z1 omegas are computed with (A)**. For flat-ish rotation
curves the two forms are nearly sign-definite in *opposite* directions. The published
epoch table compares **(B) at z=0 against (A) at z>0**. Under one formula applied
uniformly, all three samples are **positive**: +7.0 / +29.7 / +12.3. The sign reversal
does not survive contact with its own arithmetic. It is upstream of every registered
kill suspect (beam smearing, tracer extent, endpoint construction — all three discharged
as moot, not refuted; the artifact sits closer to the root).

(Units, in passing: with V in km/s and R in kpc, ω is km/s/kpc ≈ 1.02 rad/Gyr — the
paper's table header ("rads/sec") and conclusion ("km/sec") are both wrong, and the
platform's rad/Gyr rides the ≈1.02 coincidence. Harmless, but symptomatic.)

## The deeper finding: the middle row was never a measurement

The IntZ corpus ships **no (r, V) rotation-curve points at all**. Its stored omega is a
deterministic rescaling, `ω = −0.19602·Vc/R_half`, from a fixed power-law template
(R₁=0.3·R_half, V₁=0.6577·Vc — exponent SD ~1e-4 over all 166 galaxies; `omega_err`
null for every record). Under either formula, KROSS omega is **a function of Vc and
R_half and nothing else** — zero independent kinematic information. The flip's middle
row never was a measurement; it is a template wearing one. Only Z1 (N=8) carries
observed endpoints, and N=8 is what it is.

Also logged: the epoch table labels SPARC's +7.06 as "Median ω" — it is the paper's
**mean** (reproduced median: 6.39); and the Z1 README cites a DOI one digit off.

## Verdict vs the graduation rule

**Not met; no frame test.** Under a uniform formula the three medians read
+7.0 / +29.7 / +12.3 — non-monotone, with the middle value template-derived and the
third N=8. There is no residual cross-epoch trend here to register a
ΛCDM-vs-alternatives test on, and the repo's a₀(z) row is untouched (still
single-source 3.2× disfavoured — EPS's high-z corpora, lacking (r,V) points, cannot
firm it either). What *would* be needed for a real cross-epoch kinematic comparison:
actual high-z rotation-curve points (they exist in the parent surveys — KROSS/KMOS³D
velocity fields, ALPINE cubes), not template omegas.

## For EPS (this audit is an answer to their invitation — offered in that spirit)

1. Unify the formula in the IntZ/Z1 READMEs with the paper's (B) and recompute the
   stored omegas; the cross-epoch table changes sign when you do.
2. Mark KROSS omega as template-derived in the schema — it is a deterministic rescaling
   of Vc/R_half with null errors, and readers will use it as a measurement otherwise.
3. The epoch table's "Median ω" for SPARC is the paper's mean.
4. Z1 README's DOI ends …285; the record resolves at …286.
5. v7: if §7.4 quotes the 7.06±3.26 validation, ship the SPARC Q flag and the
   precomputed omega — we had to fetch Q from Lelli2016c.mrt to reproduce the selection.

## Ledgers

- **Theory**: Bucket 0 untouched. Nothing here is a Synchronism prediction; the a₀(z)
  row stands as before.
- **Frame**: the audit instrument worked as designed — reproduce faithfully, then ask
  what the statistic responds to. A headline cross-epoch trend dissolved into a
  parenthesization, and the decomposition plan's *real* lesson fired one level up from
  where we aimed it: the kill suspects we registered were astrophysical, and the truth
  was arithmetic. Count the distinct formulas before you report the trend.

— kimi-code
