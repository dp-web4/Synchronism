# Locality No-Go Prior-Art Audit — EXECUTED 2026-07-23

**Registered question**: `Research/proposals/locality_nogo_milgrom_surface_density_prior_art_audit_20260723.md`
(question and two-way verdict rule pre-fixed 2026-07-23, before this walk).
**Executor**: explorer track, synchronism-site.
**Method**: primary-source reads (arXiv PDFs via pdftotext; ADS scans; Scholarpedia; Living Reviews HTML)
plus targeted search for modern variable-comparison studies. NotebookLM unavailable (auth pending re-login).

## Verdict: Branch 2 — NO prior quantification found. With one mandatory credit enrichment (§ Template).

No published statement in the walked corpus quantifies the failure of **local volumetric density ρ(r)**
as the MOND organizing variable — not a cross-system offset magnitude, not a cluster-scale density
mismatch, not a ρ_t ∝ V⁻² scaling constraint. The phrase "volume/volumetric density" does not occur
in any walked source in an organizing-variable role. Every published organizing-variable
discrimination — from Milgrom 1983 to the 2023 machine-learning feature comparison — tests
acceleration against **{distance/radius, orbital frequency, size, surface density/brightness}**;
local volumetric density was never in any tested variable set.

Per the pre-fixed rule: the no-go gains a "prior-art audited 2026-07-23" line naming the corpus
searched, and citability **strengthens**. The preprint-strategy gate can be marked cleared on this axis.

## Per-paper findings

| # | Source | How read | Finding |
|---|--------|----------|---------|
| 1 | **Milgrom 1983a** (ApJ 270, 365) | ADS scan, OCR layer, full text | Candidate-variable reasoning present: r⁻² distance-law modification "can now be ruled out… discussed in Paper II §III." Modification keyed on acceleration. No density variable discussed. |
| 2 | **Milgrom 1983b §III** (ApJ 270, 371) | ADS scan is image-only (no OCR); content taken from three independent secondary witnesses: Sanders & McGaugh 2002 §2, Famaey & McGaugh 2012 §5, Banik & Zhao 2022 §10.2 | The famous exclusion: a modification attached to a **length scale** r₀ forces v² = GM/r₀, i.e., TF slope 2, vs observed M ∝ V⁴. Quantified, 1983. Variable is length, not density. Surface-density scalings (Σ ~ a₀/G) appear as MOND *predictions* (success side), not as a density-failure quantification. |
| 3 | **Sanders 1986** (via SM02) | secondary quote | Direction/sign statement for the **length** variable: "any modification attached to a length scale would imply that larger galaxies should exhibit a larger discrepancy — contrary to the observations." Structural twin of our V⁻² sign statement, for a different variable. |
| 4 | **Milgrom 2005, astro-ph/0510117** ("MOND as modified inertia") | arXiv PDF, read in full | Non-locality statement is *qualitative and structural*: "relativistic MI theories are non-local or non-metric," citing the Milgrom 1994 (Ann. Phys. 229, 384) Galilei-invariance theorem. The non-locality is **temporal** (kinetic action as trajectory functional), not the spatial enclosed-mass non-locality our no-go uses. Worked examples (dS embedding, Pioneer) contain no density content. Zero quantified density statements. |
| 5 | **McGaugh 2004** (ApJ 609, 652) | arXiv PDF, key sections | THE quantified organizing-variable study: discrepancy D vs **radius** (scatter), **orbital frequency** (weak), **acceleration** (tight), 1,145 points, 74 galaxies, robust to M*/L. Volumetric density not among tested variables. |
| 6 | **Sanders & McGaugh 2002** (ARA&A 40, 263) | arXiv PDF, §2 | Fig. 1: dynamical M/L vs **size** (no correlation) vs **acceleration** (M/L ∝ 1/a). "Any modification of gravity attached to a length scale cannot explain such observations." |
| 7 | **Famaey & McGaugh 2012** (Living Rev. Rel. 15, 10) | full HTML | §5: "there is no universal length scale apparent in the data (Figure 10)… This now observationally excludes all hypotheses that simply alter the force law at a linear length-scale." §4.3.3 full tested-variable list: acceleration ✓, baryonic surface density ✓, radius ✗, orbital frequency (weak). Volumetric density: absent (phrase does not occur). One footnote notes a theory can carry "a characteristic matter density ρ₀… as an additional order parameter" — a theory-construction remark, not a data exclusion. |
| 8 | **Milgrom Scholarpedia** (MOND paradigm) | full HTML | Laws list: quasi-isothermal systems must have mean **surface** densities Σ̄ ≲ Σ_M; CSDR (Milgrom 2016); central-surface-density law (Milgrom 2009b). All Σ-side successes. "A break… occurs at some critical acceleration… (and not, e.g., as one crosses a critical distance)" — qualitative, distance. Fig. 4 upper panel: η vs radius, "no correlation." |
| 9 | **Milgrom 2009b** (arXiv:0909.5184) | arXiv PDF | Entirely central *surface* density (CHSD ~ Σ_M); ρ₀ appears only inside the halo-fit product ρ₀r₀ (itself a surface density). No volumetric organizing statement. |
| 10 | **Scarpa 2006** (astro-ph/0601478) | arXiv PDF | Pressure-supported systems follow a₀; density mentions are central-density-Newtonian-regime and critical *surface* density for lensing. No ρ-vs-a discrimination. |
| 11 | **Banik & Zhao 2022** (Symmetry 14, 1331) | arXiv PDF (94k words, searched) | "It was still necessary to decide upon acceleration **rather than distance** as the crucial parameter." No volumetric-density organizing discussion. |
| 12 | **Stiskalek & Desmond 2023** (MNRAS 525, 6130) — modern extension beyond registered corpus | arXiv PDF, Table 1 | The strongest modern variable-comparison (NN + tree feature importance on SPARC): nothing improves on g_bar. Feature set: g_bar, r/R_eff, Σ_disk, Σ_bulge, Σ_tot, D, i, L₃.₆, R_eff, M_HI, T, e_N. **No volumetric density.** The nearest-miss study: had they included ρ(r), this audit would have gone the other way. |
| 13 | **Lelli et al. 2017** (ApJ 836, 152) | arXiv PDF, residual sections | RAR residuals tested vs radius, stellar *surface* density at R, baryonic mass, R_eff, effective surface brightness, gas fraction. No volumetric density. NOTE: the site's credit line "keyed on acceleration, not density (Lelli…)" is loose — Lelli tested columns (Σ), not ρ. Fix wording. |

## The Template point (mandatory credit in any writeup)

The site's one-line citable form — *a knee keyed on local volumetric density must fall as V⁻² to track
an a₀ threshold; the framework asserts V⁺²; inverted sign* — is an **instantiation of a 1983 argument
template**: pick a candidate organizing variable X, derive the mass–velocity scaling X forces, compare
to the observed BTFR/TF slope. Milgrom ran it for X = length (slope 2 vs 4, 1983b §III); Sanders 1986
stated its direction form ("larger galaxies ⇒ larger discrepancy — contrary to observations").
Nobody in the walked literature ran it for X = local volumetric density. The correct prior-art line is
therefore: *argument template Milgrom 1983b §III / Sanders 1986 (length variable); volumetric-density
instantiation first quantified here.* This is exactly the framing a referee will accept, and it
strengthens the claim by showing the method has a 40-year pedigree.

## Secondary corrections for the site (propagate)

1. **astro-ph/0510117 credit is imprecise**: its non-locality theorem is *temporal* (modified-inertia
   kinetic actions must be non-local in time, from Milgrom 1994), not the *spatial* enclosed-mass
   non-locality the no-go actually uses. Keep the citation, fix the gloss: the spatial statement
   ("the RAR is organized by g_bar, a non-local enclosed-mass quantity") is carried by the
   RAR/MDAR literature (McGaugh 2004; Lelli et al. 2016/17) + the elliptic Bekenstein–Milgrom
   field equation, not by 0510117's theorem.
2. **Lelli et al. "not density" gloss**: they tested surface columns, not volumetric ρ. Say "not
   radius, not surface brightness" or cite Stiskalek & Desmond 2023 for the systematic feature sweep.
3. **Add Stiskalek & Desmond 2023 to the audited-corpus line** — it is the modern strongest form of
   "acceleration organizes the data" and its feature set's ρ(r) gap is the sharpest evidence that
   the volumetric instantiation was genuinely unoccupied territory.

## Honest limitations of this walk

- Milgrom 1983b §III was read through three independent secondary witnesses, not the primary scan
  (ADS scan is image-only; no OCR tool available this session). 1983c not read directly. The three
  witnesses agree on §III's content; residual risk that 1983b contains an unquoted density remark
  is judged low but nonzero. If dp wants belt-and-braces before preprint submission: one manual
  read of the 1983b scan (~15 pages).
- Chameleon/symmetron screening literature (density-*screened* scalars vs MOND) is outside the
  registered corpus and was not walked; the site already cites it qualitatively as adjacent prior
  art for the wrong-variable diagnosis. If the no-go writeup claims priority for the *diagnosis*
  (not just the quantification), that vein needs its own walk. Queued as open thread.
- "Not found" is bounded by the corpus + modern extensions searched; obscure proceedings could
  still surprise. The audited-corpus line should name exactly what was searched (it does, above).

## Gate status

`stable_fixed_point_preprint_strategy.md` transferable-null #1 (locality no-go): prior-art gate
**cleared on the Milgrom/surface-density axis**, with the Template credit line required and the
screening-literature vein flagged if diagnosis-priority is claimed. Blocks lifted from
`explorer/topics/locality-nogo-standalone-writeup.md` (site repo).
