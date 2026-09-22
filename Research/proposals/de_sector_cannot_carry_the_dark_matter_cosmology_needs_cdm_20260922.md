# Proposal: the DE sector cannot also be the dark matter. The cosmology sector needs CDM. Which dark-matter story does the framework keep?

*2026-09-22, site maintainer track. Prompted by a visitor researcher persona. No bucket moves. Count stays 6, Bucket 0 = 0.
Asks dp for a frame decision.*

## The inconsistency

- `Session100_Modified_Friedmann.md` Part 2 calibrates `H² = 8πGρ_m/(3C(ρ_m))` with **C₀ = Ω_m = 0.3**. That is the Planck
  Ω_m, and ≈ 0.27 of it is **cold dark matter**. The 2026-08-12 DESI DR2 fit (γ = 0.487) also uses Planck distance
  priors, where ω_c enters through the sound horizon r_d.
- `Session241_Cosmological_Constant.md` Part 6 item 4 ("No Dark Matter Particles … direct detection experiments remain null") and
  Part 7 ("Coherence Explains Both Dark Sectors"; "Missing mass | DM particles | Coherence boost"), plus `Session277`
  P277.1 ("No dark matter particles needed"), say there are **no dark-matter particles**.
- The galaxy sector's boost ceiling `B ≤ 1/Ω_m = 3.17` uses the **same CDM-inclusive Ω_m**. So a ceiling on the boost
  that is supposed to replace dark matter is set by the dark-matter fraction.

A primary-layer grep (Research/, explorations/, manuscripts/ for baryon-only / without dark matter / no dark matter /
CDM / Ω_b) found no document that runs the DE sector with ρ_m = ρ_b, or says which ρ_m it uses. The no-particles reading
of the DE sector was **untested, not refuted**.

## Executed (pre-registered at site `71ba3e0` before the script existed)

`synchronism-site/maintainer/scripts/de_sector_without_cdm.py` (+ `_PREREG.md`, `_output.txt`, commit `5d32d70`).

With ρ_m = ρ_b, the sector's own dark component `ρ_b(1−C)/C` has to supply dark matter **and** dark energy. The
calibration C₀ = Ω_b = 0.0493 fixes x₀, which leaves γ as the only free parameter. I scanned γ = 10⁻⁴ … 3 in two
variants: the explicit argument as Session 100 writes it (x = ρ_b/ρ_crit), and an implicit one (x = ρ_tot/ρ_crit).

| | Variant A (explicit) | Variant B (implicit) |
|---|---|---|
| Identity control (CDM in, γ = ½): max \|H²/H²_ΛCDM − 1\| | 2.2×10⁻¹⁶ | — |
| γ with q₀ < 0 (accelerating today) | γ ≥ 0.0176 | γ ≥ 0.0318 |
| γ with dark/baryon at z = 1090 within 10 % of 5.36 | 0.0047–0.0058, **q₀ = +0.33…+0.36** | 0.0049–0.0062, **q₀ = +0.33…+0.36** |
| Both at once (10 %) / both at once (factor 2) | **0 / 0** | **0 / 0** |
| Best dark/baryon at recombination among accelerating γ | **1.54** (needs 5.36) | **0.79** |
| Closest joint point: max \|H/H_ΛCDM − 1\|, 0 ≤ z ≤ 1090 | 39 % | 42 % |

5 of 5 registered predictions held. The hand estimate in the PREREG (q₀ < 0 needs γ ≳ 0.016, which caps the ratio at
≈ 1.7) was right.

**Why it fails, in one line:** a log-argument tanh changes too slowly. Recombination needs a dark/baryon ratio of 5.4,
and today needs acceleration. Getting acceleration requires `d ln C/d ln x ≥ 1/3` today, and that steepness has already
used up the dark share by z = 1090. No single γ does both. This is a background-only, algebraic statement. It does not
need perturbations, CMB peaks or the fluid reading, and each of those would only add constraints (the unified-dark-fluid
P(k) argument of Sandvik+2004 is the next one).

## What it means

1. **The DE sector cannot also be the dark matter.** It works only with CDM (or some other w ≈ 0 dark component) put
   into ρ_m by hand. That makes it a dark-energy sector, not a "both dark sectors" sector. The claim that coherence
   explains both dark sectors (S241 Part 7) fails at the background level on standard numbers.
2. **This is a data-free inconsistency, not a new refutation of a registered prediction.** The count stays at 6 and
   Bucket 0 stays at 0. The DE sector stays in Bucket 3 as before; it was already ΛCDM where it lives.
3. **SPINE is not contradicted.** SPINE's dark matter is "patterns that interact **indifferently** with our matter at
   our MRH (gravitational presence, no structural coupling)". That is a particle-like, CDM-compatible reading, and the
   cosmology sector needs exactly that. What is contradicted is the **coherence-boost** reading (S241 Part 7, S277
   P277.1, and the galaxy sector's "the boost replaces dark matter").
4. **So the archive carries two dark-matter stories, and they double count.** If indifferent patterns (≈ CDM) exist, they
   are in galaxies too, and a C_a boost on top of them adds to a halo that already exists. If they don't exist, this sector
   cannot supply the dark matter at recombination (this run).

## Asks (gate on dp)

- **Frame:** which dark-matter story does the framework keep? (a) SPINE's indifferent patterns, CDM-like. The cosmology
  sector is then consistent, the galaxy-sector boost is an addition to a real halo and needs re-deriving, and
  "no DM particles" (S241 Part 6 item 4, S277 P277.1) is withdrawn. Or (b) the coherence boost. Then the DE sector needs a
  different calibration, and on today's numbers none of the tanh-log family works.
- **Archive (maintainer lane, done today as notes, text preserved):** a pointer note at the top of Session 100 Part 2
  ("C₀ = Ω_m includes CDM"), Session 241 Part 7 and Part 6 item 4, and Session 277 P277.1.
- **Ledger:** a one-line block under the DE-sector notes in PREDICTIONS.md. No bucket moves.
- **Site:** /dark-energy should say plainly that ρ_m includes CDM. The landing page's "dark matter as incomplete
  decoherence" claim should carry this result.

## Open (not run)

- Floored C (C ≥ f) in the baryons-only sector. Session 100 is unfloored, so this was out of scope. A floor binds at
  low density (late times), not at recombination where C is already large, so I expect it cannot raise the dark share
  there. That is untested.
- Whether any non-tanh C(ρ_b) could do both. That needs a step in ln x, i.e. a new parameter; see the "interior maximum"
  escape condition in the 2026-08-12 DE block.
