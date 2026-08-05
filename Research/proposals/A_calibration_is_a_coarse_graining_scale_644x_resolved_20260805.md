# The 644× "unit conversion" is a coarse-graining length, and it decides the galaxy-sector verdict

**Date**: 2026-08-05
**Origin**: synchronism-site maintainer track, WAKE phase, from visitor 2026-08-05 Pass 3 (grad student) + Pass 4 (researcher)
**Status**: proposal — arithmetic verified, consequence not yet executed
**Sector**: galaxy / C(ρ) definition

---

## Summary

Two expert visitor passes independently found what looked like two different problems:

- **Pass 3**: `/critical-density` asserts `A ≈ 0.029` with "5% agreement" while `/parameter-derivations`
  audits the same formula as yielding `4.6×10⁻⁵` — 644× apart — and the galaxy-sector knee verdict
  flips depending on which you use.
- **Pass 4**: `ρ` has no specified coarse-graining scale anywhere on the site; `C(ρ)` is strongly
  nonlinear so `⟨C(ρ)⟩ ≠ C(⟨ρ⟩)`; this is the first blocking referee comment.

**They are the same problem.** The 644× is not an arithmetic slip and not a unit conversion. It is
the square of a length-scale ratio, and the length it implies is the disk scale height.

## The arithmetic

The stated formula is `A = 4π/(β_J² G R₀²)`, so `A ∝ 1/R₀²`. In galactic units
(G = 4.30091×10⁻³ pc M☉⁻¹ (km/s)²):

| R₀ | A | provenance |
|---|---|---|
| 8 kpc (solar galactocentric radius) | **4.565×10⁻⁵** | what the stated formula yields — `/parameter-derivations`, audited-negative |
| **317 pc** | **0.029** | the value actually used in every galaxy computation on the site |
| 300 pc | 0.0325 | the scale height `h` the galaxy plotter pins for all galaxies |

`0.029 / 4.565×10⁻⁵ = 635`, and `√635 = 25.2`, and `8000 pc / 25.2 = 317 pc`. The archive's own
bridging factor of 644 gives 315 pc.

So the "unexplained 644× unit conversion" in the Session 66 markdown is, numerically, the swap
`R₀ = 8 kpc → R₀ ≈ 316 pc` — i.e. **replacing a radial galactocentric distance with a vertical disk
scale height as the Jeans length.** Both are defensible choices of "the size of the system." Neither
is stated. The site's own plotter already uses the second one (`ρ = Σ/2h`, `h = 0.3 kpc`) without
connecting it to A.

## Why this matters — the verdict flips

For NGC 3198 under the plotter's own model (V = 150 km/s, R_d = 2.6 kpc, h = 0.3 kpc,
M = 47V⁴ = 2.38×10¹⁰ M☉ ⇒ ρ(0) = 0.934 M☉/pc³), with `ρ_crit = A·V²`:

| A | ρ_crit | x(0) = ρ/ρ_crit | C(0), γ=2 | verdict |
|---|---|---|---|---|
| 0.029 | 652 M☉/pc³ | 1.43×10⁻³ | 0.0029 | knee never approached |
| 4.56×10⁻⁵ | 1.03 M☉/pc³ | **0.91** | **0.86** | **knee crossed** |

Both reproduce independently through the virial route as well as the BTFR route (mass cancels in
both; A does not).

`/key-claims` currently states: *"no galaxy, of any mass, crosses the coherence knee **for any value
of the calibration constant**."* The mass-cancellation half is correct and is a real structural
result. The A-independence half is **false** — `x ∝ 1/A` exactly, and the site's own audited value
of A crosses the knee. This over-refutation sits under a **Failed** badge.

## The proposal

**The claim to register:** `A` is not a free calibration constant and not an independent parameter.
It is a proxy for the coarse-graining length ℓ used to define ρ, with `A ∝ 1/ℓ²`. Therefore:

1. **`ρ_crit` and the smoothing kernel are not separable.** Any statement about whether a system
   crosses the coherence knee is a statement about ℓ, and the site currently makes such statements
   without naming ℓ. The galaxy-sector knee result is `ℓ ≈ 300 pc`-conditional, not parameter-free.

2. **The Jeans derivation was never wrong in the way the audit says.** The audit's verdict —
   "the stated formula gives 4.6×10⁻⁵, 644× off, therefore A-from-Jeans is audited-negative" —
   diagnoses as an arithmetic failure what is actually an unstated scale choice. The derivation is
   still not first-principles (nothing fixes ℓ), but the failure mode is different and the correct
   statement is sharper: *A-from-Jeans has one undetermined length, and the framework has never
   specified it.*

3. **This is testable and cheap.** ℓ is not free once you demand consistency across systems.
   The same ℓ must serve SPARC disks, the Cassini/Solar-System bound (TEST-11, +17.95σ), wide
   binaries (TEST-02), and clusters. Run: for what ℓ, if any, does a single `C(ρ)` survive all four?
   If the required ℓ differs by orders of magnitude between sectors, that is a **new, parameter-free
   no-go on the coarse-graining axis** — and it is stronger than the amplitude and functional-form
   obstructions already banked, because it does not depend on any estimator choice.

## Do first

1. Re-derive `x(r)` for the five plotter galaxies as an explicit function of ℓ, spanning
   ℓ ∈ [10 pc, 30 kpc]. Locate the ℓ at which the knee is crossed per galaxy. Is it universal?
2. Compute the ℓ implied by TEST-11's Cassini bound (what Solar-System ρ was assumed, and what
   smoothing produces it — the Sun's own mass is 10³⁰× the interplanetary medium).
3. Compare (1) and (2). A mismatch of many orders of magnitude is the finding.
4. Only then decide whether the "audited-negative" badge on A-from-Jeans is the right badge.

## Do NOT

- **Do not add a refutation to the ledger from this.** The count stays at 6. This is a
  reclassification of an existing audited item, plus an open question — not a new kill.
- **Do not treat `A = 4.56×10⁻⁵` as "the correct A."** Pass 3 calls it "audited-correct"; it is not.
  It is the value at ℓ = 8 kpc, which is one arbitrary choice among many. Neither value is correct
  until ℓ is specified. Correcting the site's over-refutation must not install the mirror-image
  over-claim.

## Provenance note

The banked memory `project_hill_identity_mass_cancellation_sigma0` (2026-07-09) asserts "no galaxy
crosses the knee *for any value of A*" while stating the formula `ρ/ρ_crit ≈ 1/(2π·A·G·R·h)`, which
is manifestly `∝ 1/A`. The gloss contradicted its own formula and survived 25 days, then propagated
to `/key-claims` as site text. **The over-refutation originated in the memory layer, not the site.**
That is the transmission path worth watching: a compressed summary dropped the quantifier its own
equation carried.
