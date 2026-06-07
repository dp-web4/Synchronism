# Proposal: A-from-Jeans Uses a Milky-Way-Specific Length Scale

**Submitted:** 2026-06-07
**Source:** Synchronism site visitor Pass 4 (Leading-Edge Researcher persona)
**Priority:** HIGH — this is the single surviving first-principles claim

## Finding

The derivation of A ≈ 0.0294 from the Jeans criterion uses the formula:

```
A = 4π / (βJ² · G · R₀²)
```

where **R₀ = 8 kpc** is the Sun's galactocentric radius in the Milky Way.

This value (R₀ = 8 kpc) appears inside a coefficient **claimed to be universal across all galaxies** (via ρ_crit = A · V_flat²). The problem: why would the Sun's location in the Milky Way set a universal galactic constant?

With βJ free and a single calibration point (the Milky Way), the 5% agreement (A_derived = 0.0294 vs A_empirical ≈ 0.028) is one-point calibration wearing a derivation's clothes. This is the same structural error as expressing a₀ ≈ cH₀/2π — a dimensional coincidence, not a derivation, because any O(1) dimensionless factor combined with the right dimensional scale will match.

## The Sharp Test

Re-derive A using only galaxy-intrinsic scales:
- Replace R₀ with V_flat / H₀ (Hubble scale), R_half (effective half-light radius), or R_disk (disk scale length)
- If the 5% agreement survives with galaxy-intrinsic scales → the derivation is real
- If it only works with R₀ = 8 kpc → it's a Milky-Way coincidence

## Why This Matters

The current status of the framework's physics claims:
- a₀ = cH₀/2π: Reparametrization (40-year-old Milgrom coincidence)
- Σ₀ = cH₀/(4π²G): Reparametrization (Freeman's law re-expressed)
- R₀ = V²/(3a₀): Reparametrization (follows from a₀)
- Γ = γ²(1−c): Reparametrization (textbook Palma-Suominen-Ekert 1996)
- RAR shape γ=2: Refuted (ΔBIC=+184 on SPARC ensemble)
- **A from Jeans: Only surviving first-principles claim**

If A from Jeans fails the re-derivation test, the framework has zero first-principles predictions with independent derivation. That is the honest endpoint of the physics audit.

## Recommended Site Update

1. Add explicit note to /parameter-derivations: "R₀ = 8 kpc is the Sun's galactocentric radius — justify why a Milky-Way-specific length sets a universal galactic constant, or replace with a galaxy-intrinsic scale."
2. Change ValidationBadge from "Jeans Criterion | 5% Agreement" to something acknowledging the open question (e.g., "active-mrh" status, not "validated")
3. Show the Jeans derivation steps so the dimensional bookkeeping (and βJ degree of freedom) is visible and attackable by a referee

## Session Archive Reference

Site visitor Pass 4 (2026-06-07), parameter-derivations page. Consistent with prior observation in memory: "A-from-Jeans hides R₀=8 kpc in a 'universal' constant" (memory captured 2026-06-07).

See also: `test04a_s8_receding_baseline.md` for the pattern of post-hoc calibration being framed as derivation.
