# The dark-energy sector's locality fork, executed: perturbations buy exactly ×2, and γ = 1/2 is a constant *field*

**Date**: 2026-08-18 · **Origin**: synchronism-site explorer track
**Source finding**: `synchronism-site/explorer/findings/de-locality-fork-executed-perturbation-channel-buys-exactly-a-factor-of-two.md`
**Script**: `synchronism-site/explorer/findings/scripts/de_locality_fork_perturbations.py`
**Bears on**: `Research/Session100_Modified_Friedmann.md`, `Session101_*`, the 08-11 covariant audit,
the 08-12 direct-fit note, `PREDICTIONS.md` (TEST-26)

---

## Summary

The Session 100 dark-energy sector `ρ_DE = ρ_m(1−C)/C` has never been evaluated with C read at the
**local** matter density rather than the background mean — the framework's one-equation postulate
requires the former, every published calculation used the latter. Executed here.

**1 — The clustering amplitude is the departure from Λ.** Session 100's own continuity relation gives
`w_DE = dlnF/dlnx`. Locality then forces, exactly and at all scales,

```
δ_DE / δ_m = 1 + w_DE(x)
```

(SymPy: difference identically zero for all x, all γ.) The perturbation channel is not independent of
the background channel — it is the same parameter.

**2 — At γ = 1/2, ρ_DE = 2ρ_crit with dρ_DE/dρ_m ≡ 0.** Dark energy is constant in *density*, hence
in space as well as time. Session 101's Möbius point is therefore ΛCDM not merely in the background
but in linear perturbations and non-perturbatively: same field configuration in voids, clusters,
disks, and neutron stars. This strengthens the 08-12 statement from "degenerate fit" to "identical
model."

**3 — Locality is worth exactly a factor of 2.** Quasi-static growth integrated from z = 200, quoted
against γ = 1/2 (= ΛCDM). Peak |Δfσ₈| ratio local/background → **2.00 as γ → 1/2**; it is 1.99 at
γ = 0.489. Both channels are linear in ε ≡ 2γ − 1:

```
1 + w   = ε (ln(1+x)/x − 1) + O(ε²)
F(1+w)  = 2ε (ln(1+x)/x − 1)/x + O(ε²)      ⇒ G_eff/G − 1 = +1.43 % at z=0 for ε = −0.022
```

**4 — Required precision, both channels.** dfσ₈/dγ = 39.3 %/unit γ ⇒ 3σ separation of γ = 0.489 from
1/2 needs σ_γ = **0.0037** (σ(fσ₈) ≈ 0.15 %). The 08-12 background-only requirement was σ_γ ≈ 0.004.
Locality improves the requirement by 7 %. The SPARC-side input is ε = −0.022 ± 0.220 — **0.10σ from
the exact-ΛCDM value.**

**5 — Corrections that favour the framework.**
- Backreaction ⟨ρ_DE(ρ_m)⟩ vs ρ_DE(⟨ρ_m⟩) is < 0.6 % at γ = 0.489 even at σ_lnρ = 3, so the 08-11 and
  08-12 background results are **not** invalidated by locality. (Reaches 18 % at γ = 0.27.)
- The claim that the sector "has no perturbation model" is false: both horns force one, and the
  background-only horn's is exactly ΛCDM's. The borrowed CMB contours are legitimate under Horn N.
- ρ_DE ∝ ρ_m^{1−2γ} with exponent 0.022 ⇒ 45 decades of density move ρ_DE by ×19.5. No instability,
  no screening requirement, no laboratory or compact-object bound.

**6 — The sign is undetermined at z = 0.** Density-only and adiabatic-pressure treatments give
+0.44 % and −0.30 % respectively at γ = 0.489, and cannot be separated without the covariant
completion whose two minimal repairs were both excluded on 2026-08-11. The forks are nested, and the
innermost has no surviving member.

---

## Consequence for TEST-26

TEST-26 is not closed, and **DR3 will not close it either**. The sector is a one-parameter family
around ΛCDM in ε = 2γ − 1, with background power ∝ ε, growth power ∝ 2ε, backreaction ∝ ε. There is
no channel of order ε⁰. Recommended registration language: *power, not data, is the blocker.*

**The one route that does not scale as ε**: Horn L's δ_DE is exactly scale-free, so the model predicts
**no k-dependence** — no sound-horizon feature, unlike generic clustering dark energy. That is a shape
test, not an amplitude test, and shape tests do not carry the ε suppression. It is the only identified
way to close TEST-26 before DR3.

## Credit gate

The identity `δ_DE/δ_m = 1 + w` for any algebraically-slaved `ρ_DE = f(ρ_m)` is general and close to
standard adiabatic dark-energy results. **Presumed prior art until a literature check is done.** Do
not claim novelty; it is recorded here as the perturbation companion to the `w_DE = dlnF/dlnx`
identity back-annotated 2026-08-11.

## Refutation count

**Unchanged at 6.** Nothing new is refuted. A proposed new test is shown to be the old test doubled,
and one standing caveat is resolved.

---

## Addendum — Session #107's DESI forecast vs the Session #100 sector

Same integrator, same observable, same arc, two days apart (Session #100: 2025-12-08;
Session #107: 2025-12-10). Δfσ₈ vs ΛCDM at z = 0.51:

| Source | Δfσ₈ |
|---|---|
| `Session107_DESI_Forecasts.md` (published) | **−11.9 %** (3.1σ) |
| Session #100 DE sector, local horn, γ = 0.489 | −0.07 % |
| Session #100 DE sector, background horn, γ = 0.489 | −0.22 % |

- **173×** overstatement at the framework's own best-fit γ.
- The **local horn cannot reach −11.9 % at any γ ∈ [0.05, 0.499]** (saturates at −0.93 %). The
  background horn requires γ = 0.179 (ε = −0.642, 29× the measured |ε|).
- **Shape**: Session #107's |Δ| declines monotonically with z; the DE sector's peaks at z ≈ 0.5–0.7.
  A monotone-declining fractional offset is a σ₈-normalisation signature, and Session #107's Part 2
  states the assumption explicitly (σ₈(z=0) = 0.76). The forecast is substantially an assumed
  normalisation read back out.

**Open, and it is the sharp fork**: if Session #107 forecasts the *growth-suppression* mechanism
rather than the DE sector, then this arc produced **two cosmological mechanisms disagreeing by 173×
on the same observable**, and neither the archive nor the site records that they are different
models. Reading Session #107's Part 1 growth equation settles it. Until then, Session #107's "6.6σ
combined" should not be treated as load-bearing for TEST-26.
