# ρ_crit Velocity Exponent Is Sign-Inverted for MOND-Matching (V⁻² required, V⁺² asserted)

**Source:** synchronism-site explorer session 2026-07-02
**Finding:** `synchronism-site/explorer/findings/rho-crit-velocity-exponent-mond-requires-minus2.md`
**Relates to:** `a_from_jeans_chain_of_custody_closure.md`, `density_compander_nogo_locality_classification.md`

## Result

The critical density `ρ_crit = A·V_flat²` used in the C(ρ) galaxy sector has **three
mutually-incompatible velocity-exponent provenances and no derived one**:

| Provenance | Exponent | Status |
|---|---|---|
| Stated "Jeans" derivation (Session 66) | ρ_crit ∝ V^**+0.5** | only reading that hits the empirical coefficient; uses galaxy-intrinsic length |
| Framework code / usage (`equations.ts`) | ρ_crit ∝ V^**+2** | asserted; coefficient carried from the V^0.5 law with mismatched units |
| **MOND-matching requirement (this note)** | ρ_crit ∝ V^**−2** | derived; profile-independent |

### Derivation of the MOND requirement

For C(ρ)'s density threshold (ρ = ρ_crit) to coincide with MOND's acceleration threshold
(g_bar = a₀), identify ρ_crit with the baryonic density at the transition radius r_t where
g_bar(r_t) = a₀. Using BTFR (`M_bar = V⁴/Ga₀`) and `r_t = V²/a₀`:

```
mean-density:  ρ ~ M/r_t³ = (V⁴/Ga₀)/(V²/a₀)³ = a₀²/(GV²)   ∝ V⁻²
isothermal:    ρ = V²/(4πG r_t²) = a₀²/(4πGV²)               ∝ V⁻²
```

Both give exactly −2. Structural: r_t ∝ V² and M ∝ V⁴ ⇒ ρ ∝ V⁻², independent of baryon profile.

### Magnitude

Required ρ_crit ~ 0.01–0.3 M⊙/pc³ (V = 50–300 km/s) — galactic-outskirt / mean-disk densities,
exactly where a modification should switch on. Framework's 0.029·V² gives 70–2600 M⊙/pc³
(inner-galaxy values), **240×–300,000× too high, gap growing with V**. With ρ_crit = 652 M⊙/pc³
(V=150) the entire luminous disk (ρ ≈ 0.01–100) sits at C ≲ 0.28: the coherence knee is never
crossed, so C(ρ) crossing ρ_crit is not what produces the rotation-curve transition.

## Why it matters for the research core

This is the **locality no-go (Milgrom-non-locality instance) resolved on the velocity axis**. The
existing statement is a magnitude one — a global ρ_crit(V) leaves a ~1.7 dex cross-system offset
between local ρ(r) and the non-local g_bar the RAR tracks. The velocity-exponent form is sharper
and harder to dismiss:

> **A modification keyed on local volumetric density must have a per-object threshold falling as
> V⁻² to track an a₀ acceleration threshold. This is forced by BTFR + the a₀ threshold and is
> independent of the modification's functional form. C(ρ) asserts V⁺² — the wrong sign.**

That one sentence constrains not just C(ρ) but *any* ρ(r)-keyed MOND mimic, and belongs in the
same transferable-nulls drawer as the locality classification criterion. It is a candidate line
for PREDICTIONS.md's galaxy-sector self-audit and for any preprint built around the local-density
no-go.

## Recommended PREDICTIONS.md / catalog update

Add to the galaxy-sector closure notes: *"ρ_crit ∝ V² is not derived (Jeans gives V^0.5) and is
anti-derived by MOND-matching (which requires V⁻²). The exponent has three incompatible values and
the sign required by the only test it faces is opposite to the one asserted."*
