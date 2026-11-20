# The dSph cap exclusion used Wolf's projected-radius coefficient on a 3D radius; the reinvention claim needs an ablation, not only a post-cutoff arm

*Site maintainer, 2026-11-19. From the 2026-11-19 visitor log (simulated graduate-physics and researcher personas) and explorer
2026-11-18. No bucket moves; count 6; range 4–8; Bucket 0 = 0. Items marked **gates on the project lead** need a decision.*

## 1. Correction: the dwarf-spheroidal exclusion of the 6.39 cap is 4/3 too strong in mass (executed, not pre-registered)

A researcher persona asked which mass estimator stood behind "Draco ≈ 100, Crater II ≈ 60" and guessed the site's choice
was the conservative one. It is the most aggressive of the three I checked, by a slip.

- `synchronism-site/explorer/findings/scripts/pressure_supported_boost_ceiling.py` (2026-08-23) says "Wolf+2010:
  M_1/2 = 4 σ² r_1/2 / G" and feeds it the **3D** half-light radius. Wolf et al. 2010's 4 goes with the **projected** R_e;
  with the 3D radius the coefficient is 3. So every B_req in the 08-23 table, and the 11-15 6.39 test that reused its
  inputs unchanged, overstate the dynamical mass by 4/3.
- Re-read under both published estimators (`synchronism-site/maintainer/scripts/cap6p39_vs_dsph_estimator.py`, inputs
  unchanged; identity row reproduces 101.5 / 60.2 and 4.9σ):

| | site 11-15 (4·r3D) | Wolf 2010 (3·r3D) | Walker 2009 (2.5·R_e) |
|---|---|---|---|
| Draco B_req (M/L 2) | 101.5 | **76** | 48 |
| Crater II B_req (M/L 2) | 60.2 | **45** | 28 |
| Draco z vs 6.39 at M/L 4 | 4.9 | **4.5** | 3.7 |
| Crater II z vs 6.39 at M/L 4 | 4.9 | **4.2** | 2.9 |
| Cap excluded up to (Draco, M/L 4, Šidák 3.51 / local 3σ) | 14.6 / 18.5 | **11.0 / 13.9** | 6.9 / 8.7 |
| Escape M/L (Draco / Crater II) | 9.2 / 7.0 | **6.9 / 5.3** | 4.3 / 3.3 |
| σ tidal inflation needed (Draco / Crater II) | 1.51 / 1.32 | **1.31 / 1.15** | 1.04 / 0.91 |

- **What survives:** under the standard estimator both dwarfs still exclude 6.39 past the six-system threshold, so the 11-15
  verdict ("the most permissive cap is excluded by matter dynamics, no light rule") stands. **What changes:** "excluded up to
  ≈ 14.7" becomes ≈ 11; the Crater II leg needs only a 15% tidal inflation of σ to escape, and Crater II is a known tidally
  disturbed system, so it is the weaker leg; under Walker's estimator Crater II drops out and Draco excludes 6.39 only just.
  The exclusion is estimator-conditional at its edge. Site corrected on six surfaces.
- **Pattern (third of its kind):** a docstring naming the right paper made the formula look checked. The coefficient was taken
  from the paper and applied to a different radius definition. Same family as 2026-10-19 "a correction is a claim" and the
  catalogue-column lesson: read the variable the coefficient belongs to, not just the citation.
- **For the project lead:** the 10-04 (ii) and 11-15 (ii) statements about the no-CDM cap 1/Ω_b = 20.3 used the same inputs
  and are 4/3 too strong in mass; not recomputed here (routed to the explorer).

## 2. Explorer 2026-11-18 drained: the LLR register bound on a measured gradient

The site's `/for-researchers#gdot-llr` assumed a density e-fold over ℓ. The explorer measured d ln ρ_s/dt from the Gaia DR3
census to 30 pc: 2.3×10⁻⁵ yr⁻¹ at ℓ = 1 pc (α Cen, Barnard's Star), 4–8×10⁻⁸ at 10 pc (√N averaging). The 2σ bound is
|d ln G_eff/d ln ρ_s| ≲ 1×10⁻⁸ (1 pc), 1.5×10⁻⁷ (3 pc), 3–6×10⁻⁶ (10 pc), ≈ 5×10⁻⁵ beyond ~30 pc. The 11-18 (ii) "1×10⁻⁸–1×10⁻⁷"
was right at 1 pc and 25–50× too strong at 10 pc. The 68/74 count is unaffected. Prior art for both ingredients
(Berengut & Flambaum 2010; Olive & Pospelov 2008) is now credited; the combination was not found in six sources.

## 3. Where the density law's +2843 comes from: a tilt, not the outskirts (post-hoc)

A graduate persona: with one global h the law is Σ-keyed, and g_bar ≈ 2πGΣ in a thin disc, so it should track the
acceleration law in inner discs. Decomposed on the exact 2,438 points
(`synchronism-site/maintainer/scripts/dbic2843_where_it_comes_from.py`; one of my two written predictions held):
- by quarter of each curve, inner → outer: 18 / 26 / 26 / 31 %;
- by g_bar: 33 % from g_bar > 10⁻¹⁰ and 39 % from g_bar < 10⁻¹¹;
- mean residual of the density law: **−0.21 dex** (over-predicts g_obs) in bright inner discs, **+0.21 dex** in the deep regime;
  the acceleration law is within 0.01 dex at both ends.

So the loss is a tilt. At the fitted γ = 0.046, C ≈ 0.046 ln(ρ/ρ_c) varies too little across SPARC to give a boost near 1 in bright
centres and near 10 in outskirts at once. Σ-keying does not track g-keying even where the disc is thin and bright. Also
recorded: the 122-vs-166 sample difference is the loader's V_flat > 0 requirement (SPARC Table 1).

## 4. Cassini and the compander: the width floor is 0.16 dex from the MOND function that passes

The site said the compander "cannot reach δ = 6" because its 10–90 % width never falls below 1.17 dex, without the δ = 6 width.
Computed (`synchronism-site/maintainer/scripts/delta6_width_vs_compander_floor.py`, identity: simple μ = tanh at γ = ½ = 1.908 dex):
δ = 6 is 1.006 dex, δ = 4 1.060, δ = 1 1.725. The width gap is small. The larger difference is the Newtonian return: power law
(1 − C ≈ 2x^(−2γ)) against exponential; at the same a₀ and g_N = 2a₀, ν − 1 is 6×10⁻⁵ (δ = 6) against 3×10⁻³ (γ = 3). The site now
prints both. Whether the 11-07 "no sharp member" reading should say "no exponentially returning member" is a wording choice for the
B2/TEST-25 notes (gates on the project lead with 10-16 (i)).

## 5. WAKE: the reinvention finding needs an ablation (research direction)

The researcher persona named the one A2ACW datum it would cite: four published theories re-derived by same-corpus adversarial
sessions and never cited in 3,308 sessions, every match arriving through a lane with literature retrieval. It also named the
missing control. The 10-12 CBP note already showed the noticing lanes share a *task* ("what does this reduce to?"), which
A2ACW never asks. Together these say the cheapest discriminating experiment is a 2×2, not the post-cutoff arm:

| | no retrieval | retrieval |
|---|---|---|
| adversarial roles (A2ACW as is) | the record (0/4) | arm B |
| single agent, same code and budget, asked "what known result does this reduce to?" | arm C | arm D |

Items: the four documented identities (simple μ, Cardassian, Refracted Gravity, Collins LIV), each posed from the archive text
that preceded the identification, plus matched decoys with no published twin. Score: names the right paper (with a locator)
or not; decoys score false positives. Prediction to register first: D ≈ B > C > A. If C ≈ D, the question is the active ingredient
and retrieval is secondary; if B ≈ A, adversarial roles add nothing to novelty-checking. Cost: tens of sessions, in-family
(no cross-vendor scorer needed), so it does not wait on the gated items. **Gates on the project lead** only in the sense that
it competes with the post-cutoff arm's default date (2026-12-09); I recommend running this first because its outcome is
unknown and it separates two explanations the record currently confounds.

## 6. Frame, for the project lead

The four personas agree on one pattern: the derivation and ledger layers are current; the summary layers (landing equation
block, glossary entries, tool beginner boxes, breadcrumbs) lag them, and each visitor pass finds a new instance. Today's
sweep fixed seven and added 13 lint rules, plus a JSX check that caught the HTML5-only entities rendering as literal text on
/dark-energy. This is now the dominant defect class on the site, and lint only catches a wording after someone notices it.
The structural fix is the one `labelDefinitions.ts` already applied to badges: headline verdict sentences (γ = 2 pin, γ width,
density law, ρ input) single-sourced in `ledger.ts` and rendered, not retyped. I will continue moving them there; no decision
needed unless you prefer otherwise.

DR4 is 13 days out; the DR4 default and the dark-matter fork still wait (11-13, 09-22).
