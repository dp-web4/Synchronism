# Proposal: SPINE's "instruments, not refutation" sentence; the galaxy floor *is* the dark-energy C today; five refutation scopes audited in both directions

**Date:** 2026-09-11
**From:** synchronism-site maintainer track (WAKE phase, written before any site edit)
**Status:** routed to dp. **No bucket moved. Executed count stays 6. Bucket 0 stays 0.**
**Evidence:** `simulations/floor_is_cosmic_C_and_w_sign.py` (+ `_output.txt`), mirrored from
`synchronism-site/maintainer/scripts/`. Site-side reads of `synchronism-site/explorer/findings/scripts/`
(`l2_sparc_core.py`, `l2_field_equation_on_sparc_output.txt`, `sparc_pinned_at_rg_knee_l2_output.txt`).

Input: the 2026-09-11 four-persona visitor log. The cross-persona pattern was one object, several surfaces:
one kill criterion with two verdicts, one environment run under three IDs, one field equation read two
ways, one no-go scoped wider than its computation. Most of that is site work and was fixed there. Seven
items reach the research core. Items 1–3 bear on direction; 4–7 are scope and ledger governance.

This session moved items in **both** directions, and I am saying so because a session that only moves one
way should be suspected. Items 1 and 2 cut against the framework. Items 4, 5 and 7 are over-refutation
corrections: they narrow what was claimed refuted. Item 3 hardens a no-go that was already expected.

---

## 1. SPINE's summary sentence contradicts Bucket 2 (anti-oscillation, undersell of the negatives)

SPINE §"An invitation" says: *"The physics is **untested for lack of means, not refuted on the merits.**"*
The site had compressed this to *"the honest reason is lack of instruments, not refutation"* on the landing
page and `/fundamentals`. A researcher persona read it against the scoreboard directly beneath it and was
right: TEST-09 refuted, on the merits, what the ledger calls the framework's only structural difference
from MOND. SPINE's own paragraph draws the correct line (borrowed data *can* refute a derivation, and
has; it cannot *confirm* a novel prediction), but the bolded sentence erases that line. It is the sentence
that gets quoted.

**Proposed SPINE wording:** *"The quantitative probes were refuted where borrowed data could reach them.
The one move, the ontology, is untested for lack of instruments, not refuted on the merits."*

The site now carries this two-clause reading. **The SPINE edit gates on dp.** SPINE is canonical and I did
not touch it.

## 2. The galaxy floor is the dark-energy sector's C at today's mean density, identically

**Identity (verified at five γ; exact to printed precision).** Session 100 closes the DE sector on flat
ΛCDM: ρ_DE/ρ_m = (1−C₀)/C₀ = Ω_Λ/Ω_m. That gives **C_DE(ρ̄_m,0) = Ω_m for every γ**, the same number the
galaxy sector uses as its floor. (PREDICTIONS already calls Ω_DE = 1−Ω_m a tautology of the calibration;
the identity itself is not new. What is new is noticing that the floor and the calibration are one number.)

**It reconciles a live contradiction in the ledger.** The PREDICTIONS TEST-09 row says *"Its floor Ω_m is
genuinely derived from cosmology."* Site `/parameter-derivations` and `/for-researchers` say *"1/Ω_m is
nowhere derived."* The identity says exactly what is true: the floor is **identified with** the DE
sector's calibration constant, and that constant is itself set by hand. So the floor is not a derived
number. Proposed ledger wording: *"floor identified with C_DE at the present mean density (a calibration
identity), not derived."*

**If the identification is taken literally, it makes a prediction.** Read the floor as C_DE at the
*ambient* cosmic density, not as a fixed 0.315, and the boost ceiling evolves:

| z | B_max (γ = 0.487) | f_DM,max | B_max range, γ ∈ [0.3, 2] |
|---|---|---|---|
| 0 | 3.175 | 0.685 | 3.175 (identity) |
| 1 | 1.279 | 0.218 | 1.06–1.45 |
| 2 | 1.085 | 0.078 | 1.00–1.20 |
| 4.2 | 1.017 | 0.017 | 1.00–1.07 |

At fixed z the same reading makes the ceiling environment-dependent. Voids at δ = −0.8 give B_max = 11.8,
which is what SPARC's dwarfs demand. Group-scale overdensities give B_max ≈ 1.1.

**Two gates, both named.** (a) *Two knees.* The DE knee is 4.2×10⁻⁸ M☉/pc³ and the galaxy knee is 0.161
(3.9×10⁶× higher). A single C with the DE knee is exactly Newtonian in every disc (C = 0.999999 at
0.1 M☉/pc³). The identity therefore needs two coherence functions, one as the other's floor. That is the
existing "the coherence function is not one function" problem arriving from cosmology. (b) *The MRH of
the floor.* Nothing in the framework says which smoothing scale defines a galaxy's "ambient" density.
Without it, the environment version is not a bet.

**Why this matters for direction.** SPINE says door #3 (secular / time-domain) is "an untested direction,
not a registered bet — no specific falsifiable prediction yet." This is a specific number on that door.
It is conditional on gate (a), and at z ≳ 1 it barely depends on γ. It also cuts against the framework.
The galaxy sector is already refuted on SPARC with the floor too high at z = 0, and an evolving floor
makes high-z discs *more* Newtonian, so this is a test of what the floor *means*, not a rescue.

**Cheapest check, on existing data:** published dark-matter fractions within R_e at z ≈ 1–2.5 (Genzel et
al. 2020; Price et al. 2021; Nestor Shachar et al. 2023). Any robust disc with f_DM(<R_e) above
1 − C_DE(z) (0.22 at z = 1, 0.08 at z = 2) refutes the literal reading. **Recommendation:** do not
register this. The explorer runs the literature check first (topic seeded). If it dies, record it as an
elimination of the only floor interpretation that puts a number on door #3. If it survives, dp decides
whether gate (a) is acceptable.

## 3. TEST-26: the no-go is exact for the substituted sector

A visitor researcher persona derived it; I verified it numerically against the continuity equation (not
the closed form). With u = 1 + ρ_m/ρ_crit and C = tanh(γ ln u):

    1 + w = F(u)/(u^{2γ} − 1),   F(u) = u^{2γ} − 1 − 2γ(u−1)u^{2γ−1},   F(1) = 0,
    F′(u) = −2γ(2γ−1)(u−1)u^{2γ−2}   ⇒   sign(1 + w) = sign(1 − 2γ) at every redshift.

The check found 0 sign violations over 11 γ × 51 z, with max |numeric − closed| = 1.4×10⁻⁶. At γ = 0.487:
w₀ = −0.992 and wₐ = +0.015. The family moves along (1+w₀) ∝ +wₐ, which is why the 2026-08-12 likelihood
put it at Δχ² ≈ 0. It also answers the persona's question 8: **the Ω_m floor never binds on the past light
cone** (C ≥ Ω_m for all z ≥ 0), so flooring opens no escape. **Recommendation:** put this into the TEST-26
statement as the model-class clause for the algebraic sector. That turns "kill-or-tie" from a class
argument into a theorem. The covariant completions (A, B) are unaffected. Gates on dp, like TEST-26 itself.

## 4. The Refracted-Gravity comparison: the charge was wrong, the scope was still too wide

The researcher persona charged that "3 to 17× worse at RG's f = 0.089" came from the algebraic g_bar/C
reduction, making it a straw version of RG. **That is false, in two ways.** (i) The numbers come from
`l2_sparc_core.py`, a finite-volume solve of ∇·[C∇Φ] = 4πGρ on the axisymmetric (R, z) half-plane,
refraction term included, validated against exact Hankel-transform discs. (ii) RG at its **own**
published parameters was also run with that solver (2026-08-28, not refit). With Υ profiled, χ²/N is
188 / 240 for the two DiskMass sets and 716 / 911 / 1252 for the Cesare+2020 E0-floor sets. MOND is 21.25
and Newton 465, and RG beats MOND in 10–17 % of galaxies.

**But the site stated neither**, and one sentence was over-scoped. `/for-researchers` artifact 1 limited
the local-density no-go to "algebraic coupling — the class this framework belongs to", and the same page
names the field-equation form (which has a gradient coupling) as the framework's completion. Both are
now scoped on the site. **What remains open, and should be on the ledger as open:** a SPARC *refit* of
RG's three parameters (the steepness exponent carries a 2.30× ln-vs-log ambiguity across RG papers, per
the 08-27 Publisher note), and the striction force of the variational completion (08-26), which neither
RG as published nor any of these runs include. Until both run, "density-keyed permittivity on SPARC" is
an executed grid, not a class theorem.

## 5. Bell B1: the class claim is three constructions; the missing primitive has a known price

The site said *"No construction reaches the Tsirelson bound without signaling"* and counted the substrate
as refuted "by execution." **Toner & Bacon (PRL 91, 187904, 2003):** local hidden variables plus *one
hidden bit* of communication per trial reproduce the singlet correlations exactly, with uniform marginals,
so no signaling is observable. No such construction was built here. The global-clock construction's
signaling was observable (`04_global_clock_chsh.py` measures Alice's marginal against Bob's setting), so
that result stands as recorded.

Two consequences. **(a) Ledger precision:** B1's criterion ("violation only appears with signaling")
should say *observable* signaling. Hidden-communication substrates are untested, not refuted. SPINE's
own question asks about a **purely local** substrate, which excludes Toner–Bacon, so SPINE's framing is
already right; only the site's class wording was wrong. **(b) Direction, and the reason this is in a
proposal rather than a site note:** SPINE names the gap as *"a non-relabelable, conditional
setting-dependence — a primitive the ontology does not contain and would have to derive."* Toner–Bacon is
the minimal known instance of that primitive, and it costs exactly one bit per trial. A single substrate
updating in parallel on a global tick does not obviously forbid a hidden conditional bit. Whether the
ontology *wants* one is a choice about Bohm-like nonlocality, not a loophole in Bell. That is now the only
unexplored branch of "the one test that matters," and it has a price tag. Worth dp's attention.

## 6. Test-ID collision between this ledger and the site

PREDICTIONS Bucket 2 files the 2026-07-14 SPARC × Cosmicflows-4 environment run as **"TEST-08"**. On the
site, TEST-08 is the **Freeman-law** card, and the environment run is TEST-03s (a substitute for the
never-run ALFALFA TEST-03, run-as-registered against Session 177). The site's own 2026-08-10 rule is one
flat namespace, no ID denoting two things. A grad-student persona's "run as registered vs NEVER RUN"
contradiction was this collision plus two different registrations. The site now states both. **Proposed:**
relabel the ledger row "S177 environment run (site TEST-03s)". Archive-lane label fix; routed.

## 7. Globular-cluster fork: two disclosures the P611.2 note should carry

Verified in `explorer/findings/globular-cluster-knee-test-executed-…md`:

1. The verdict cutoffs (ok ≤ |MOND+EFE|, marginal < 2×, excluded otherwise) were **set at execution** on
   2026-09-07. The Session 611 registration has no quantitative criterion.
2. **Mass profiles were not refit per gravity law.** They are catalogue profiles normalized to Baumgardt &
   Hilker's Newtonian N-body masses. The slope statistic cancels normalization, not profile shape.

The ordering claims survive: density-keyed at γ = 0.489 is about 2× MOND+EFE under the same bias. The word
"marginal" should not be read as a pass. At γ = 2 the density law does no better than Newton or MOND+EFE.
**Proposed:** append "(post-hoc benchmark; Newtonian-fitted mass profiles)" to the P611.2 note in
PREDICTIONS. Routed.

---

## What I did *not* do

- I did not edit SPINE.md or PREDICTIONS.md. Every item above is a site fix or a question to dp.
- I did not register the evolving floor (item 2). It has a named gate, and a literature check that could
  kill it cheaply comes first.
- I did not count anything. The environment row, Bell, the RG grid and the GC fork all stay where they are.
