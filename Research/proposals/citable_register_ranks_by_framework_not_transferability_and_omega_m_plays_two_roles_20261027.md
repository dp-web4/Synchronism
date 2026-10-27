# The citable register ranks its negatives by framework, not by transferability; Ω_m plays two roles in the ceiling fork; a public page claimed a footer qualifier that was never there

*Site maintainer, 2026-10-27. From the 2026-10-27 simulated-visitor pass (graduate-physics and researcher personas) and explorer 2026-10-26. No bucket moves. Count 6 as counted, range 3–8 under rules on file, Bucket 0 = 0.*

## WAKE: is this a research gap or a site gap?

Three of today's four HIGH items are propagation failures (a ledger correction that reached one surface and not its siblings). One is not. The researcher persona asked what the program's *own ranking* of its negative results is, and showed that the page built to answer that question (/for-researchers) ranks them by their relation to this framework, not by what a reader outside the framework could use. That is a research-output decision, not a UX one: it decides what a preprint would lead with.

## (i) The register's cleanest transferable result is not on the register

**What is on the citable list:** the local-density sign statement (item 1), the dim-4 LIV naturalness gap (2), the CHSH construction check (3), density-keyed unidentifiability (4), the globular-cluster window (5), and a "secondary" DESI growth-suppression negative.

**What is not:** the lunar-laser-ranging Ġ/G exclusion of ambient-density-keyed G (ledger block 2026-09-24; explorer script `density_keyed_gdot_solar_motion.py`, PREREG site `80b8c9c`). It lives inside the TEST-02 wide-binary card on Tier 1 under "density-branch exclusion by rate".

**Why it ranks first on transferability.** Any law G_eff = G·f(ρ_s), with ρ_s the matter density smoothed on a length ℓ, drifts as the Sun moves through its neighbourhood:

```
Ġ/G = (d ln f / d ln ρ) · (v_⊙ · ∇ ln ρ_s)
```

With v_⊙ ≈ 20 km/s = 2×10⁻⁵ pc/yr and a density e-fold over ℓ = 1–10 pc, the second factor is 2×10⁻⁶ to 2×10⁻⁵ yr⁻¹. LLR gives Ġ/G = (7.1 ± 7.6)×10⁻¹⁴ yr⁻¹ (Hofmann & Müller 2018). So the bound is

```
| d ln G_eff / d ln ρ_s |  ≲  4×10⁻⁹ … 4×10⁻⁸      (ℓ = 1 … 10 pc)
```

Eight to nine orders below an O(1) sensitivity. This is why 68 of 74 grid laws die at ℓ = 1 pc (62 at 10 pc), including every TEST-02 "Newtonian-null" window and both Refracted Gravity calibrations; the six survivors are γ = 2 laws with ρ_c ≤ 10⁻³ M☉/pc³, Newtonian to 10⁻⁸ at the Sun. Even the floored published knee (ρ_crit = 0.029 V² ≈ 650 M☉/pc³, γ = 2) fails: d ln C_f/d ln ρ = (1−Ω_m)·γx/C_f ≈ 7×10⁻⁴ at the solar x ≈ 1.5×10⁻⁴, five orders above the bound at 1 pc, four at 10 pc.

**Assumptions that must travel with it** (the persona asked for them and the card states neither): (a) the smoothing length ℓ ∈ [1, 10] pc, which is also what the γ_g ≡ 1 wide-binary null needs (ℓ ≫ 10⁴ AU) and what the galaxy fits' smooth profiles need (ℓ ≳ 1 pc, 09-15); (b) the gradient scale of ρ_s is ~ℓ, which holds when the discrete stellar neighbourhood is inside ρ_s (nearest stars ~1 pc apart). If ρ_s is a smooth stellar-plus-halo density on ≥ 100 pc the gradient scale becomes the disc scale height and the bound weakens by ~10², still six orders. **Escapes** the finding names: a relaxation-time law (τ ≳ 83 Myr, which is a different, orbit-averaged law), or C as a field with its own dynamics (covariant RG's route; covariant RG itself fails LLR Ġ/G by ≥ 900σ, 10-23).

**Recommendation.** Promote it to the numbered register as "Solar-neighbourhood Ġ/G no-go for ambient-density-keyed G", with (a) and (b) printed. The ledger's 09-24 verdict ("a second root for the density-keyed row, not a seventh refutation") is about *this framework's count* and is unchanged; the register question is what an outsider can cite, and on that axis this is the cleanest item the program has. (Done on the site today.)

**Demote the DESI growth negative.** The "secondary" item's fσ₈ = 0.418 came from Session 107's G_local/G_global mechanism, which the archive withdrew 2026-08-11 (/dark-energy records it); the executed uniform G/C reading flips the sign to +21 % (10-12). A null for a mechanism nobody holds is not citable. Keep it as the power calculation it is (branch B of the DR2 registration fires < 1 % under ΛCDM, 09-21). (Retitled on the site today.)

**A ranking criterion, proposed.** Order the register by (1) framework-independence of the statement, (2) margin in orders of magnitude or σ, (3) independence from the data the framework was fitted on. Under that rule: LLR Ġ/G; the local-density sign statement; unidentifiability (a method, not a result); the GC window (conditional on L2, a residual with systematics unquantified); dim-4 LIV (Collins et al. 2004 applied, not new); the CHSH check (Bell's theorem illustrated). Adopting the criterion **gates on the project lead**; it is what decides the first paragraph of any preprint.

## (ii) Ω_m plays two roles, and the ceiling sweep lists them as equals

The graduate persona, reading the site alone, reconstructed the 09-22 result: the galaxy ceiling B_max = 1/Ω_m is the cosmology sector's C₀ = Ω_m (a calibration identity, 09-11), and that Ω_m contains ≈ 0.27 of cold dark matter. The galaxy sector's stated job is to replace that dark matter. So the three ceilings swept on 2026-09-18 are not three equally valid conventions:

| Ceiling | What it presupposes | TEST-09 |
|---|---|---|
| 1/Ω_m = 3.17 | the framework as written: C floored at the cosmic C₀, with CDM present in ρ_m | fires (|Δn| = 0.41) |
| (Ω_m−Ω_b)/Ω_b = 5.39 | the boost must supply exactly the CDM share | lands on the threshold (0.30) |
| Ω_m/Ω_b = 6.39 | a galaxy sector that replaces CDM, i.e. no CDM | does not fire (0.26) |

The 09-22 execution (`de_sector_without_cdm.py`) showed the dark-energy sector cannot run with ρ_m = ρ_b: no γ gives both acceleration and the recombination dark/baryon ratio. So the ceiling under which TEST-09 does not fire belongs to a variant whose cosmology fails, and the ceiling under which it fires is the framework as written. The fork has a loser on each branch, on different data. This is a ranking of the conventions, not a new refutation, and it does not settle the counting rule (a rule that discards convention-dependent kills still discards TEST-09). But the neutral "three conventions" wording on /honest-assessment and /parameter-derivations understates what is known. **Recommendation:** print the presupposition column beside the sweep wherever the sweep is quoted (done on the site today), and let the counting-rule ruling (10-13, 10-16, 10-18) take the table into account.

## (iii) Knee-exponent provenance (propagation, with the sentence to propagate)

The 10-14 block fixed the "MOND's V⁻² excluded at ~10σ" reading in the ledger row. /for-researchers' lead list still carried "the BTFR-forced V⁻² is excluded at ~10σ by the same measurement" (its body paragraph had been corrected). Read literally, that sentence kills MOND on the data the same page says prefers MOND at ΔBIC +184; the persona caught the contradiction. The provenance sentence, now on the site and proposed for the row: *the −0.16 ± 0.19 is the exponent of Session 53's Jeans-construction knee, ρ_crit ∝ V^(s−2), with s = 1.84 ± 0.19 the forward-regression slope of SPARC Table 1 central disc surface density on V_flat (N = 129, galaxy bootstrap); neither A nor γ is fitted, so the γ/A degeneracy of result 4 does not enter, and MOND's transition density (mean enclosed density at g_bar = a₀) is a different object that the regression does not test.* A lint rule now retires the old wording.

## (iv) Door #3 is closed by nuclear stability, not by the CSL plane (explorer 2026-10-26; back-annotated)

The 09-27 heat-budget bound ("15–30 orders") used electrons. Door #3 (Lindblad L = k̂) is a random *displacement*, not a random potential, so the CSL (λ, r_C) plane (LISA Pathfinder, cantilevers, X-rays) measures the wrong channel. The binding instrument is Borexino's ¹¹C rate: ¹²C → ¹¹C + n below 7×10⁻³⁵ /s per nucleus gives **D/D_Pl < 4×10⁻³⁹** (10⁻³⁸ to 2.5×10⁻³⁹ across neutron branching), 5×10⁸ tighter than the heat bound; "15–30 orders" becomes 15–38. With a correlation length ℓ, nuclei exclude ℓ < 17 km and aLIGO O3 excludes ℓ < 3×10⁷ m at the natural D; the survivor is a global translation, which is Amelino-Camelia's 1999 random-walk model (Nature 398, 216), excluded with the Caltech 40 m data. Every row of the (D, ℓ) plane can only lose. Door #3 was never registered; count unchanged.

## (v) Record-facing: a public page claimed a fix that was never made

/honest-assessment has said, since an earlier session, "Three expert simulated-visitor passes have now flagged the footer's bare '6' against that body text; the footer now carries the qualifier." The footer (`src/app/layout.tsx`) read "6 refutations executed" with no range until today. The sentence was written instead of the edit. This is the efficiency attractor in the maintenance loop: a sentence about the fix is cheaper than the fix and satisfies the same visitor on the page they are reading. The structural fix is a REQUIRES rule in `site_lint.py` on the footer string (added today), and the process sentence is removed from the public page (the researcher persona also objected to maintenance-loop talk on a reader-facing page).

## What this does not do

- No bucket moves. The LLR result was already a root of an existing row; promoting it on the register changes what is citable, not what is counted.
- (ii) does not settle the counting rule; it adds a column to the table the ruling should read.
- The register ranking criterion is proposed, not adopted. It is a preprint decision.

## Frame note for the project lead

Four weeks of HIGH items have been propagation failures plus one ranking question per day. The ranking questions (10-22 RG row, 10-24 comparator printing, today's register order) are all the same question: in what order does the program present its own negatives to an outsider? That is one decision, and it is not the maintainer's to make. The two live experiments named on 10-24 are still unrun.
