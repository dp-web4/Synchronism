# Proposal: Status of Synchronism's Cosmology Domain After Stress Tests

**Filed by**: Maintainer track, based on 2026-04-26 visitor feedback (researcher persona)  
**Date**: 2026-04-26  
**Type**: priority-shift

---

## The Question

After the Bullet Cluster structural failure and the galaxy rotation MOND-reparametrization diagnosis, what cosmology claims — if any — remain that are both (a) genuinely novel and (b) not yet falsified?

---

## Background

### What was claimed

The site's cosmology domain made three load-bearing claims:

1. **Dark matter mechanism**: Synchronism maps dark matter via CFD viscosity — low coherence regions (C = 1/μ_eff) are high-viscosity, which causes gravitational cohesion analogous to dark matter halos.

2. **Galaxy rotation curves**: Fits to 14,760 galaxies (SPARC + ALFALFA-SDSS) using the McGaugh 2016 RAR interpolating function plus an environmental scatter ansatz. The /galaxy-rotation page carries a "Strongly Supported" tile.

3. **a₀ = cH₀/(2π)**: The framework's preferred acceleration scale is claimed as a derivation from coherence principles.

### What audits found

**Bullet Cluster structural failure (March 2026)**  
The CFD viscosity mapping has a sign error. Synchronism predicts dark matter should be MORE viscous (high viscosity = low coherence = C = 1/μ_eff). The Bullet Cluster shows dark matter is LESS collisional than baryons: σ/m < 0.47 cm²/g (Harvey 2015), i.e., effectively inviscid. This is not a parameter mismatch — the mapping predicts the wrong direction. The failure is documented on /key-claims and /dark-matter-failure.

**Galaxy rotation: MOND reparametrization (ongoing)**  
The interpolating function used IS the McGaugh 2016 RAR function. The Synchronism contribution is an environmental scatter ansatz. R² for environment-dependent scatter = 0.14 — weak signal. ΔBIC vs. a MOND-only baseline has not been computed. Without ΔBIC, the galaxy rotation fits are evidence for RAR/MOND; they are not yet evidence for Synchronism over MOND.

**a₀: dimensional reparametrization**  
a₀ = cH₀/(2π) is shared with McCulloch 2007, Verlinde 2017, and Smolin 2017. The derivation routes differ but the result is the same expression. This is a dimensional coincidence or common underlying physics — it does not distinguish Synchronism from those frameworks.

### The propagation failure

The Bullet Cluster failure is documented internally (/key-claims, /dark-matter-failure) but has NOT propagated to /galaxy-rotation. A reviewer arriving at /galaxy-rotation — the natural entry point for the cosmology domain — sees "Strongly Supported" without any indication that the dark matter mechanism that underlies the viscosity-coherence mapping is structurally broken. This is not a site aesthetics issue; it is a scientific integrity issue.

---

## Why This Matters for the Research Direction

The cosmology domain is the site's highest-visibility domain for expert reviewers. It is also where the most failures have accumulated. The current state — "Strongly Supported" on a page that depends on a broken mechanism, ΔBIC not computed for the primary observational claim — puts the research program at risk of being dismissed not for the failures themselves (which are documented) but for the apparent unawareness of them.

More fundamentally: if both the dark matter mechanism (sign error) and the galaxy rotation MOND-baseline comparison (not computed) are unresolved, **the cosmology domain currently has no confirmed novel prediction**. The 500 Mpc cosmic interference scale (TEST-07) lacks a derivation (see `cosmic_interference_500mpc_derivation.md`). The BAO shift prediction lacks a derivation of the predicted magnitude. The /galaxy-rotation "Strongly Supported" badge is load-bearing credibility for a claim that has not cleared its MOND baseline.

This needs to be resolved, not labeled.

---

## Proposed Investigation

**Step 1: CFD viscosity mapping — does any dark matter mechanism survive?**

The sign error invalidates the specific mechanism (viscosity → dark matter stickiness). The question is whether any alternative coherence-based mechanism can account for dark matter behavior without requiring the same viscosity direction. Specifically:
- Is there a reformulation in which low-coherence regions are LESS collisional (consistent with Bullet Cluster σ/m < 0.47 cm²/g)?
- Does the coherence gradient — rather than the coherence level — provide such a mechanism?
- If no reformulation is consistent, document this as a refutation and remove dark matter from the active claims list.

Expected output: either a candidate reformulation (session + derivation) or a documented refutation that closes the dark matter arc.

**Step 2: Compute ΔBIC for galaxy rotation vs. MOND-only baseline**

The SPARC + ALFALFA-SDSS fits exist. The McGaugh 2016 RAR baseline is fully specified. ΔBIC is computable from existing data.
- Implement a MOND-only RAR fit on the same 14,760 galaxy sample
- Compare residuals, log-likelihoods, and BIC with the Synchronism environmental ansatz
- If ΔBIC < 2: the environmental ansatz adds no information — demote from "Strongly Supported" to "Reparametrization"
- If ΔBIC > 6: genuine improvement over MOND — the environmental scatter claim is real

Expected output: a ΔBIC value, properly interpreted, and a badge update on /galaxy-rotation.

**Step 3: Audit remaining cosmology predictions for novelty**

Survey the full cosmology claim list (cosmic interference, BAO shift, large-scale structure, a₀ derivation) and classify each:
- **Novel and unfalsified**: has a derivation, makes a symmetric prediction, not shared with MOND/ΛCDM/other frameworks
- **Reparametrization**: correct prediction but result follows from MOND, ΛCDM, or dimensional analysis without additional Synchronism content
- **Refuted**: contradicted by data (like Bullet Cluster)
- **Unanchored**: prediction stated without a derivation (like 500 Mpc scale)

Expected output: an updated cosmology scorecard — honest classifications, no softening of the failures.

**Step 4: Propagate known failures to /galaxy-rotation**

Regardless of Steps 1-3 outcomes, the /galaxy-rotation page needs to reflect the Bullet Cluster structural failure NOW. The mechanism that the page implicitly relies on is broken. Update the validation badge, add a prominent note linking to /dark-matter-failure, and specify what the page's "Strongly Supported" claim is actually contingent on (RAR fits only, not dark matter mechanism).

This step is site maintenance and can be executed by the maintainer track immediately; it does not require the research results from Steps 1-3.

---

## What Resolution Looks Like

**If the dark matter mechanism can be reformulated (Step 1 succeeds)**:  
A new mechanism is documented in the archive, the Bullet Cluster failure is relabeled as a refutation of the SPECIFIC CFD viscosity mapping (not dark matter in general), and the reformulated mechanism carries a "Speculative" badge pending new tests.

**If the dark matter arc is closed (Step 1 fails)**:  
Dark matter is removed from active cosmology claims. The Bullet Cluster entry is labeled "Refuted" and moved to the honest-assessment record. The cosmology domain's remaining scope is galaxy rotation (pending ΔBIC) and large-scale structure predictions (pending derivations).

**If ΔBIC is computed and is large (Step 2 confirms environmental scatter)**:  
/galaxy-rotation gets a narrower, more defensible badge: "Reparametrization with environmental extension" — honest about what is and isn't novel. The claim is real but not strong evidence for Synchronism over MOND.

**If ΔBIC is small or negative (Step 2 finds no improvement)**:  
/galaxy-rotation badge becomes "MOND Reparametrization." The fits are correct but they are MOND fits, and Synchronism's environmental scatter adds no detectable signal. This is a clean, honest result — more useful than a "Strongly Supported" badge that can't survive a journal reviewer's baseline comparison.

**Minimum outcome regardless of above**:  
Step 4 is non-negotiable. The /galaxy-rotation page must reflect the Bullet Cluster failure within the current site cycle. Leaving "Strongly Supported" on a page that depends on a broken mechanism is the one outcome that is straightforwardly wrong.

---

## Risk If This Proposal Is Wrong

If the audits are too harsh and some cosmology claim is genuinely stronger than this assessment suggests, the risk is: honest relabeling of a real result as "reparametrization," which can be corrected when ΔBIC or a derivation is produced. That's recoverable.

The risk of NOT doing this: a researcher lands on /galaxy-rotation, sees "Strongly Supported," digs in, finds the Bullet Cluster failure documented elsewhere on the same site, and concludes the research program lacks internal consistency. That's not recoverable in the same conversation.

The failures are already documented. The only remaining task is to propagate them honestly to where a reviewer would look.
