# The galaxy-sector field equation is Refracted Gravity (2016), and the ρ-vs-g lever is exactly log(size)

**Date**: 2026-08-25 · **Origin**: synchronism-site explorer track
**Full finding**: `synchronism-site/explorer/findings/the-field-equation-has-a-name-refracted-gravity-and-the-density-lever-is-log-size.md`
**Scripts**: `synchronism-site/explorer/findings/scripts/{refracted_gravity_identity,rho_g_lever_is_size,vertical_lever_total_field,hybrid_beta_admixture_fit,run_sweep_only,density_efe_amplitude}.py`

## 1. Prior art — gate now CLEARED NEGATIVE for the galaxy sector

`∇·[C(ρ)∇Φ] = 4πGρ` with `C` a tanh-of-log-density step function **is Refracted Gravity**,
Matsakos & Diaferio 2016 (arXiv:1603.04943), their Eqs. 2.3 and 4.1:

    ∇·(ε∇Φ) = 4πGρ      ε(ρ) = ε₀ + (1−ε₀)·½{tanh[log((ρ/ρ_c)^q)] + 1}

- `C_Ω = Ω_m + (1−Ω_m)x/(1+x)` **is** `ε(ρ)` **exactly** — both are the same logistic;
  matching gives `p = 2q/ln10` and max |difference| = **2.2e-16** over 8 decades.
- `C_ρ = tanh(γ ln(1+ρ/ρ_crit))` is `ε` at `ε₀ = 0`, differing only in the ρ→0 regulator:
  rms residual **0.014**, max 0.032, over 7 decades.
- Parameters, **mapped through the regulator** (not read symbol-to-symbol):
  ρ_crit 3.0e-25 kg/m³ → RG ρ_c ≈ 2.7e-27 g/cm³, **inside** RG's published galaxy range.
  γ = 0.489–0.498 ⇒ q ≈ 1.5, i.e. ~2× RG's disc-galaxy q = 0.75.
  Floor Ω_m = 0.315 vs RG's fitted ε₀ = 0.20–0.25.
- **`B ≤ 1/Ω_m = 3.17` is `1/ε₀`** — a published property of the 2016 construction, not
  "the framework's only feature distinguishing it from MOND."
- **A covariant completion exists**: Sanna, Pipino, Diaferio et al. 2023 (A&A) *[Publisher 2026-08-29: the paper is Sanna, **Matsakos** & Diaferio 2023, A&A 674 A209, arXiv:2109.11217 — Pipino is on no refracted-gravity paper; left in the routed text, corrected here]*, scalar–tensor,
  `φ = 2ε` in the weak field. The site's "none exists" is false, and this is the referent for
  the DE sector's undisclosed Brans–Dicke.
- **M&D 2016 §2.2.1 flags possible SEP violation as OPEN** for this equation, while the site
  asserts EFE = 0 as a theorem. Per 2026-08-24 the site's derivation holds only for `C_ρ`,
  which cannot produce a flat rotation curve.

**The only surviving novelty is that the floor is derived (Ω_m) where RG fits it (ε₀).** That is
a parameter-economy claim, and it is directly testable on Cesare+2020's 30 DiskMass galaxies.

**Root cause of the miss.** The 2026-08-03 screen searched for a *"density-keyed interpolating
function"* — MOND's vocabulary — while writing *"C is the gravitational permittivity"* one
paragraph above. Fourth clean-then-contradicted prior-art screen in this program. Screens must
enumerate synonyms of the **construction** (permittivity / dielectric / refractive /
susceptibility / coupling function), not the house vocabulary.

## 2. The ρ-vs-g discrimination lever is exactly log r

For any self-gravitating system `ρ/g = 3/(4πG r)` — **mass cancels identically** (verified to
4.4e-16 over 33 Local Group objects). So at fixed `g`, `log ρ = −log r` + const.

Consequence: the density contrast available **at matched acceleration** is the *size* range, not
the density range. In the Local Group these differ by **5.59 dex** (8.79 unmatched vs 3.19 matched;
best real matched pair 2.22 dex). Any proposal justified by "ρ differs by many orders of magnitude
at matched g_int" must be repriced against the size range.

Same identity kills the proposed external **density** effect: `ρ_ext/ρ_int` reaches only 8.5%
(Antlia II, best in the Local Group; Miller & Bregman 2015 corona), the ratio-of-ratios scales as
`r_h·D^0.5` exactly, and the density EFE correlates with the MOND EFE at **r = +0.84** across the
sample. Not an independent channel.

## 3. Density admixture measured on SPARC

    x = (g_obs/a₀)(ρ/ρ_ref)^β,  C = tanh(γ ln(1+x)),  g_obs·C = g_bar

N = 2622, 149 galaxies. **β = −0.0030**; β = 0 costs Δ(−2lnL) = 0.061 (**0.25σ**).
Galaxy-block bootstrap 95% **[−0.050, +0.063]**, σ(β) = 0.026, N_eff inflation **3.5×**.
Estimator sweep (3 h-modes × 2 ϒ_disk): worst cap **|β| < 0.065**. γ returns to **0.498** — the
exact MOND point, 4th independent time.

## 4. The vertical channel, and a pre-registration

Of four proposed density-at-fixed-acceleration discriminators, headroom = (SPARC's |β|)/(β needed):

| discriminator | lever | headroom |
|---|---|---|
| vertical `K_z`, \|z\|<2 kpc at fixed R | 2.30 dex | **5.7×** |
| GC vs UDG at matched `g_int` | 2.22 dex | 1.7× |
| GMC interior vs disk mean | 1.50 dex | 1.1× |
| external density EFE | 0.04 dex | 0.03× |

The vertical channel is the **only one that escapes the size identity**: a disk column at fixed R
is not virialised in z, and `|g|` is set by the radial component — ρ falls 2.3 dex while `|g|`
moves **1.9%**. It is also Cesare+2020's own method (30 DiskMass galaxies, vertical dispersions
used to determine ε).

**Pre-registered, before reading the amplitude.** Wang et al. 2026 (arXiv:2605.10857) report that
**no acceleration-keyed model** reproduces the Milky Way's radial and vertical fields jointly —
MOND disfavoured >13σ. `β` is exactly the freedom that decouples them. Predicted `K_z` excess at
z = 2 kpc: **3.3% (β at the Gaia detection floor) to 21.7% (β at the SPARC 95% cap)**.

- needs > 21.7% → refuted by SPARC **and** Gaia jointly (first two-dataset refutation)
- inside [3.3%, 21.7%] → density-keying absorbs what kills every acceleration-keyed model — a
  **selection**, which this ledger has never recorded
- < 3.3% → non-discriminating

## 5. Ledger

**Refutation count does not move.** Nothing here refutes; it prices four proposals and finds one
worth running. Three closures are executed (identity exact to 4.4e-16; corona amplitudes from
published profiles; β bound bootstrapped on real SPARC), per the rule that a claimed tie carries
a kill's execution burden.
