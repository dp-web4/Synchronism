# The modified-G reading of the DE sector fails lunar laser ranging; the density fit's response is g_bar's

*Site maintainer, 2026-11-05. Sources: a simulated-visitor researcher persona (site visitor log 2026-11-05, Pass 4,
HIGH) and site explorer 2026-11-04 (`synchronism-site/explorer/findings/density-response-is-g-bar-leakage-and-the-08-02-null-never-reached-the-site.md`,
PREREG site `4891f7e`). No bucket moves; count 6; range 3–8; Bucket 0 = 0.*

## 1. The mean-density horn has two readings, and LLR removes one

The DE sector is `H² = 8πGρ_m / (3C(ρ̄_m))` with C₀ = Ω_m (Session 100). The 2026-09-15 locality trilemma treated
"mean-density only" as one horn. It is two:

- **Fluid (Cardassian) reading:** gravity is unmodified, the extra term is a dark-energy fluid. Ġ = 0. Survives; it is
  Cardassian-class prior art (2026-09-14) and ΛCDM-degenerate (2026-10-09).
- **Modified-G reading:** G_eff = G/C(ρ̄_m(t)) is *the* gravitational coupling, the reading the "uniform G/C" growth row (U)
  uses (2026-10-07). If that G binds the Earth–Moon system, it drifts with the cosmic mean density:

  Ġ/G = −Ċ/C = 3H₀ · d ln C / d ln ρ.

  At γ = ½, C = x/(x+2) and d ln C/d ln x = 1 − C = 1 − Ω_m exactly, so Ġ/G = 3(1 − Ω_m)H₀ = **1.42×10⁻¹⁰ yr⁻¹**
  (H₀ = 67.4). LLR: (7.1 ± 7.6)×10⁻¹⁴ yr⁻¹ (Hofmann & Müller 2018, as quoted on the site; not re-read here). Excluded by
  ~2,000× the central value and by **≥ 527× the 2σ upper edge for every γ in 0.3–3** (d ln C/d ln ρ runs 0.57–0.88).

Script: `synchronism-site/maintainer/scripts/de_modified_g_reading_vs_llr.py` (+ `_output.txt`). Predictions written in its
header before it ran (the γ = ½ number by hand); 3/3 held; identity control 5×10⁻¹¹.

**Consequence.** (U) is now excluded twice: by the ISW sign and lensing (10-07, needs a light-deflection choice) and by the
Moon (needs nothing but the background). Every G/C reading that touches local dynamics is LLR-excluded: the local-density
one (2026-09-24, 68/74 laws) and the mean-density one (here). What remains of the DE sector is ordinary dark energy in
Cardassian clothing.

**Scope.** This assumes the cosmic G/C is the G of bound systems. A reading in which G/C acts only on the Friedmann
equation and on nothing else is the fluid reading by another name. No archive document asserts the local reading, so this is
not a count change under the current rule; it joins the 10-12 question of whether executed exclusions of unasserted
readings count.

**Open (for the explorer):** prior art. Modified-gravity Cardassian variants may already note the Ġ problem (I recall
Gondolo & Freese discussing fluid vs modified-gravity interpretations; that is recall with no locator, not verified). Also: binary-pulsar and BBN Ġ/G bounds beside LLR
(the persona's request).

## 2. "Weak, not off" is a property of the fit; the data's density response is g_bar's

Correction to the 2026-11-04 block (i). The explorer re-ran the question on the exact sample behind the +945 / +2843:

- Partial slope of log B on log ρ at fixed g_bar: −0.001 (accel law refit around it; ΔBIC +7.7 *against* adding it) and
  +0.003 (cubic baseline), galaxy-bootstrap 95% CI −0.022…+0.027. Positive controls at −0.05 and −0.10 are recovered.
- A model that never sees ρ reproduces the data's marginal ρ slope (−0.211) to 0.2%. corr(log g_bar, log ρ) = +0.84.
- The γ = 0.046 density fit's ρ^−0.10 tracks about half of that leakage.

The 2026-08-02 form-free null (σ(log B | ρ)/σ(log B | g_bar) = 1.37×, 1.34–1.77× across nine treatments; r_partial
+0.0012) already said this. It reached the whitepaper (section 11) but no site page, for three months. It is now on the
site's /honest-assessment (`#density-form-free-null`) and /rar-scatter, and a lint rule requires every "ρ^−0.10" to name
g_bar on the same line.

## 3. Process: the orphan rate is unmeasured

Two orphaned explorer results have now been found by accident (the 09-02 lensing P0 at 23 days; the 08-02 null at three
months). The explorer showed a filename-based drain check is mostly wrong (5/6 spot-checked "orphans" had reached the site by
content). The right instrument extracts each finding's headline numbers and greps `src/`. Seeded to the explorer.

## For the project lead

1. Whether to print the executed exclusions outside the count (11-04 (ii): §D.3 EdS, completion B, formula-A fit; now also
   the modified-G reading here) as a descriptive line beside the 6. The site now prints the 6 *by kind* and names the lensing
   exclusion against rows #1–#2; it does not print the outside-the-count list.
2. Whether the DE sector's status should be stated archive-wide as "Reparametrization; fluid reading only" (the site's
   /dark-energy now says so, replacing three inconsistent statuses).
