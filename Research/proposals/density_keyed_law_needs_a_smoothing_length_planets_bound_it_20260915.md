# A Density-Keyed Law Needs a Smoothing Length, and the Planets Bound It From Below

**Filed**: 2026-09-15 (site maintainer track, WAKE phase — before any site fix)
**Status**: proposal. Scope condition on an existing statement plus a missing parameter. **No bucket moves. Count stays 6. Bucket 0 stays 0.**
**Raised by**: visitor 2026-09-15, Pass 3 (grad student, "has the density-keyed law ever met the Solar System?") and
Pass 4 (researcher, "interplanetary density at Saturn sits below RG's ρ_c").
**Bears on**: PREDICTIONS.md Bucket 1 note "📌 UPDATE 2026-09-09" ("The Sun and the clusters do not close this sector"),
the 09-08 P611.2 GC fork, the joint local window published on the site (/honest-assessment#gc-fork, /for-researchers),
and the Refracted Gravity identification (08-26/08-27).
**Pre-registration**: `synchronism-site/maintainer/scripts/density_keyed_law_vs_interplanetary_medium_PREREG.md`
(commit `0ad68f3`, made before running anything)
**Script**: `synchronism-site/maintainer/scripts/density_keyed_law_vs_interplanetary_medium.py` (+ `_output.txt`)

---

## The question nobody had asked

Every density window in the archive evaluates ρ as a **smooth** density:
- the Oort limit is a disc average;
- the GC window uses cluster-averaged densities;
- SPARC uses exponential-disc mass models.

None of them says over what length the density is averaged. Primary-layer grep for "interplanetary|solar wind" finds
only the explorer's 08-05 coarse-graining finding. That finding treats the framework's **virial** law ρ_crit = A·V²,
where ℓ cancels and the Cassini kill is ℓ-independent. It does not treat a **universal** knee ρ_c, which is what the
09-09 window, RG, and the GC fork use.

Read pointwise, the Sun's own neighbourhood is the solar wind:
- ρ ≈ 0.14 M☉/pc³ at 1 AU (n_p = 5 cm⁻³, He included);
- falling as r⁻², to 1.6×10⁻³ at Saturn and 1.6×10⁻⁴ at Neptune.

**That range brackets every knee in the published joint window.**

## Executed (rule fixed before running)

Spherical flux conservation for ∇·[C(ρ)∇Φ] = 4πGρ is exact outside the Sun: g(r) = GM☉/(C(ρ(r)) r²). Statistic
D = |C(1 AU)/C(Saturn) − 1|, the fractional disagreement between GM☉ inferred from Earth and from Saturn. Excluded if
D > 10⁻⁶. That threshold is roughly 10⁴ looser than ranging allows.

| Family | Grid | Excluded | Typical D |
|---|---|---|---|
| Framework floored `f + (1−f)tanh(γ ln(1+ρ/ρ_c))` | f ∈ {0.089, 0.315}, γ ∈ {0.489, 2}, ρ_c ∈ [3.2×10⁻⁴, 0.161] + window edges, n_p ∈ {3, 5, 10} | **192/192** | 0.4–7 at the window edges |
| RG Eq. 4.1, ρ_c = 0.0083 | ε₀ ∈ {0.20, 0.25}, Q ∈ {0.1, 0.5, 1, 2}, ln and log₁₀ | **48/48** | 0.12–4 |

At the joint window's own edge (γ = 0.489, ρ_c = 0.0079, f = 0.315): C = 0.93 at Earth and 0.38 at Saturn, so the
two orbits disagree on G by a factor of 2.5. For RG to pass at ρ_c = 0.0083 it needs **Q ≲ 3.7×10⁻⁷**, i.e. no step at all.

**Coarse-graining.** With a ball of radius ℓ centred on the planet, D stays O(1) until ℓ exceeds the planet's
heliocentric distance. Once the ball swallows the Sun, every planet sees the same Sun-dominated density, C becomes
uniform, and the effect is absorbed into GM☉. So **the Solar System requires ℓ ≳ the orbit of the outermost ranged body
(≳ 30 AU; spacecraft tracking pushes it further out).**

## Prior art: RG anticipated this and never did it

Matsakos & Diaferio 2016 (arXiv:1603.04943, §2.2.1), read from the full text:

> "by introducing the permittivity ε as a function of density ρ, we unavoidably introduce a smoothing length D over
> which we average the density field. We … postpone the investigation of this crucial topic to future work; here, we
> simply speculate that the scale of D is expected to be of astronomical interest, for example tens of astronomical
> units or larger, because the mass discrepancy problem appears to be relevant on scales larger than the solar system."

So this is **not a new no-go against RG**. It turns their speculation into a lower bound and confirms its order.
Whether later RG papers (Cesare+2020; Sanna+2023) fixed D is **not checked**, and is listed below.

## What it changes in this archive

1. **The 09-09 line "The Sun and the clusters do not close this sector" is conditional.** It holds for a smoothing
   length at which the solar neighbourhood reads as the Oort disc average, ℓ ≳ interstellar spacing (~1 pc). **Read
   pointwise, the planets close it at O(1).** Neither reading appears anywhere in the archive. Proposed back-annotation:
   add the condition to that note.
2. **The window has a two-sided constraint on ℓ that has never been written down:**
   - *Lower:* ≳ 30 AU from the planets. SPARC's smooth stellar profiles implicitly need ≳ 1 pc: a 30 AU ball in a
     disc almost never contains a star, so ρ would be ISM gas, not the modelled stellar density.
   - *Upper (plausible, not computed):* ≲ GC half-mass radii (a few pc), else cluster densities dilute into the field
     and the GC window moves.
   - *A detail that may matter:* at ℓ ≈ 1 pc centred on the Sun, the Sun itself contributes ≈ 0.24 M☉/pc³, which
     exceeds the Oort-limit disc density (~0.1). The "solar neighbourhood density" is then ℓ-dependent inside the
     window's own range.
3. **For Synchronism specifically, the MRH is the concept that ought to supply ℓ.** The framework says relevance is
   bounded by a Markov Relevancy Horizon, and a density-keyed law needs exactly one such length. This is the first
   place a SPINE-level concept is forced to carry a **number**. If the MRH cannot be given a length that sits in
   [~1 pc, ~few pc] for stars, and is consistent with the Solar System, the density branch has a free parameter it
   never declared. *Stated as a question, not a claim.*

## Not claimed

- **Not a seventh refutation.** No archive document commits to a pointwise ρ. Count 6.
- Not a claim about RG beyond the quoted 2016 text.
- The upper bound on ℓ from GCs is not computed.
- Whether a smoothed-ρ law remains Φ-independent (so EFE = 0 survives; see 08-24) is presumed but not re-derived.
  Smoothing ρ does not touch Φ, so it should.

## Asks

- **Explorer:** (i) compute the GC-side upper bound on ℓ, and whether [ℓ_min, ℓ_max] is non-empty at the 09-09 window;
  (ii) check whether Cesare+2020 / Sanna+2023 fix RG's D; (iii) recompute the Oort window at ℓ = 1, 3, 10 pc with the
  Sun included.
- **Back-annotation (maintainer, today):** the 09-09 note gets its condition.
- **dp (gated):** whether "declare the smoothing length" becomes a required field of any density-keyed registration
  (TEST-26 does not use ρ locally; P611.2 does).
