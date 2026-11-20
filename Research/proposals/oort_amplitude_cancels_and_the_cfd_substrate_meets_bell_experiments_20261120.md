# Proposal: the Oort amplitude cancels against laboratory G, and the CFD substrate meets the Bell experiments, not only Bell's theorem

**Date:** 2026-11-20 · **Author:** site maintainer (Claude, synchronism-site maintainer track)
**Sources:** simulated-visitor researcher and graduate personas (site `visitor/logs/2026-11-20.md`, Passes 3–4); site
explorer 2026-11-19 (dSph); maintainer run `synchronism-site/maintainer/scripts/oort_vertical_gradient_density_keyed.py`
(predictions committed at site `21b842f` before running, 3/4 held). Abstracts read 2026-11-20 by a fetch agent, quoted
verbatim: Hensen et al. 2015 (arXiv:1508.05949), Bancal et al. 2012 (arXiv:1110.3795), Salart et al. 2008
(arXiv:0808.3316), Bovy & Rix 2013 (arXiv:1309.0809), McKee, Parravano & Hollenbach 2015 (arXiv:1509.05334).
**Bucket moves:** none. Count 6, range 4–8, Bucket 0 = 0.

---

## 1. The SPINE question this touches (frame first)

SPINE's "one test that matters" asks whether a **local** substrate, measured through an observer-pattern's phase-lock,
can reproduce quantum correlations without signalling, and says that a "no" makes Synchronism "an ontology /
interpretation, not a local-hidden-variable physics". Ledger row #6 (Bell/CHSH) is booked as a *construction check*
(the substrates built give S = 1.85–2.0) and the hidden-communication class (Toner–Bacon) as *untested*.

A researcher persona pointed out that both bookings are weaker than the record allows.

- **The local substrate is excluded by data, not only by theorem.** Loophole-free experiments measure S > 2:
  "We perform 245 trials testing the CHSH-Bell inequality S ≤ 2 and find S = 2.42 ± 0.20 … p = 0.039"
  (Hensen et al. 2015, Nature 526, 682; photonic tests the same year agree, not re-read here). A substrate capped at
  S ≤ 2 is excluded by that measurement, the same way TEST-25 is excluded by Cassini: inherited, but on data.
- **The hidden-communication class is constrained, not untested.**
  - Bancal et al. 2012 (Nature Physics 8, 867): models with influences at "any finite speed v with c < v < ∞" predict
    correlations "that can be exploited for faster-than-light communication"; they "exclude any possible explanation of
    quantum correlations in terms of influences propagating at any finite speed" (given no signalling).
  - Salart et al. 2008 (Nature 454, 861): in any privileged frame in which Earth moves slower than 10⁻³ c, the
    influence "would have to exceed that of light by at least 4 orders of magnitude".
- **What the framework's own update rule allows.** FUNDAMENTALS: "all cells simultaneously evaluate neighbor tensions
  and step forward together". A nearest-neighbour parallel update carries influence at most one cell per tick, i.e. at
  ≤ c on a Planck grid. So the CFD substrate *as defined* is local, and is excluded by the experiments above. A
  Toner–Bacon bit across spacelike separation needs a different update: global and instantaneous in the preferred
  frame (finite speeds are closed by Bancal). That is the Bohm class with a preferred foliation. It is prior art, and
  it is compatible with the absolute tick, but it is not what FUNDAMENTALS writes.

**So the honest form of SPINE's "no" branch is sharper than "it is an interpretation".** An interpretation built on a
nearest-neighbour substrate is empirically excluded. The interpretation survives only if the substrate gets a global,
instantaneous update rule, at which point its entanglement story is Bohm's. SPINE's "distance is scan-delay" sentence
reads like a global scan, while FUNDAMENTALS (and SPINE's own note) say the grid updates in parallel and the CRT models
only the observer's sampling. Those two readings now give opposite verdicts, and the framework should pick one.

**For the project lead (gates):**
1. Which update rule is the framework's: nearest-neighbour parallel (local, excluded by Hensen et al.), or global
   instantaneous in the absolute frame (Bohm-class, prior art)? A SPINE sentence either way.
2. Re-label ledger row #6 from "construction check" to "excluded on data (Bell experiments), inherited from Bell's
   theorem", parallel to row #5. This changes the count's *texture*, not the count. The site now prints that reading as
   a description (landing table, `ledger.ts` COUNT_BY_KIND_SENTENCE) and leaves the count at 6.
3. The B1 row's "hidden-communication substrates untested" should read "constrained by Bancal 2012 / Salart 2008 to an
   instantaneous update". Row edits gate on the project lead (with 10-16 (i)).

## 2. The Oort argument had the wrong mechanism (verdict partly survives)

The site said the density-keyed boost at the Sun (2.7–10.5 in g at the archive's 0.028·V^0.5 knee, unfloored) is
excluded "by the Oort limit". The persona's objection: under the 1–10 pc smoothing that the wide-binary null and the
LLR item need, laboratory G is measured at the same ρ_s, so G_lab = G/C(ρ_s). An Oort analysis that divides by
4πG_lab recovers the baryonic midplane density exactly. The amplitude cancels, as it does in the wide-binary ratio
and in ephemerides (an argument the site itself makes for those two). Only the boost's change with height is visible.

**Run (estimate, not a registered test).** 1D slab, K_z = 2πG Σ_b(<|z|)/C(ρ_b(z)) (algebraic and L2 coincide in 1D),
lab G at the midplane, no dark matter. Baryons: midplane 0.084 M☉/pc³ (McKee+2015: all matter 0.097 ± 0.013, DM
0.013); 38 shapes (gas + thin + thick exponential layers) with Σ_b(<1.1 kpc) in 40–56 M☉/pc². Datum: Bovy & Rix 2013,
Σ(|z| < 1.1 kpc) = 68 ± 4 M☉/pc².

| law | C(ρ₀)/C(ρ(1.1 kpc)) | Σ_dyn(1.1 kpc) | vs 68 ± 4 |
|---|---|---|---|
| archive knee 0.028·V^0.5, γ = 2, unfloored | 18.5–43.9 | 886–1783 | +9.5…+11σ |
| same, γ = 0.489 | 19.3–45.7 | 924–1858 | +9.5…+11σ |
| Jeans A = 4.6×10⁻⁵·V², γ = 2 | 20.7–49.3 | 992–2001 | excluded |
| units hybrid 0.029·V², γ = 2 | 21.2–50.2 | 1012–2042 | excluded |
| one global knee (+2843 fit: 4.4×10⁻⁶, γ = 0.046) | 1.40–1.59 | 61–67 | passes (−1.0…−0.1σ) |
| archive knee, γ = 2, affine Ω_m floor | 1.70–1.74 | 69–81 | passes (+0.2…+1.7σ) |
| archive knee, γ = 0.489, affine floor | 1.18–1.19 | 48–57 | short by 2–3σ |
| control C = const (Newton, no DM) | 1.00 | 40–48 | short by 3.6–4.9σ (Bovy & Rix's DM detection) |

Identity control: the inferred *midplane* density equals ρ₀ for every law (the amplitude cancels), verified
numerically. Predictions: P1 (archive knee unfloored excluded > 5σ) held; P2 (global-knee fit within 3σ) held; P4
(identity) held; P3 ("floored variants 1–5σ high") **failed**: γ = 2 floored passes, γ = 0.489 floored falls short.

**Reading.** The exclusion of the archive-knee unfloored law at the Sun survives, through the gradient (≈ 10σ, a
significance set by the baryon census, not by the 18–44× ratio). The floored variants and the free-global-knee fit are
*not* Oort-excluded; they die elsewhere (SPARC +2843 and LLR for the density-keyed family). No row moves. The site's
"Oort-excluded amplitude" wording is replaced on five pages.

Caveats: a slab, not a Jeans analysis of tracers at several heights; one datum; the gas column (10–14 M☉/pc²) is
recalled; Σ error scaled with the ratio. A proper K_z(z) run is routed to the explorer.

## 3. Record items (site, already fixed)

- **Quadrature does not bound the boost.** Under v² = v_b² + (V_flat·C)², C ≤ 1 bounds the *added velocity* by V_flat,
  but B = 1 + (V_flat·C/v_b)² is unbounded as v_b → 0. The 3.17 cap is a property of floored division only;
  quadrature's failure is circularity. /key-claims and /parameter-derivations said otherwise.
- **The Fisher +1.000000 (γ, A) non-identifiability is a units-slip result.** It was measured at A·V², A = 0.029,
  where every disc sits in the linear regime and only γ/A enters. With one free global knee γ and ρ_c are separately
  identified (corr +0.89, 11-04 block); at 0.028·V^0.5, 45/122 discs cross the knee. Untested there, not refuted.
- **The LLR register item does not test Refracted Gravity as published.** It smooths over 1–10 pc; RG's ε is pointwise
  and sits inside the divergence, and read pointwise at the Earth–Moon system it keys on the interplanetary medium. The
  pointwise check on file is the 09-15 Earth-vs-Saturn GM test (RG passes only for Q ≲ 3.7×10⁻⁷, one in-repo reading of
  Eq. 4.1). Outsiders were about to cite "LLR excludes RG".
- **dSph leg (explorer 2026-11-19).** Crater II is observed disrupting (Limberg et al. 2025, as returned by search: a
  ≈7σ velocity gradient along its debris track; central σ 2.51 km/s), so it leaves the equilibrium leg; the cap-family
  exclusion by matter dynamics rests on Draco alone (≈ 7 Walker – ≈ 11 Wolf). Dwarfs do not close 1/Ω_b = 20.3
  (2.5–2.8σ).
- **TEST-04a** carried "Failed", "Untested", and "Kill Triggered" on three pages; now one badge (untested,
  underpowered as registered) and one sentence for its three σ figures.
- **C_a is not a tanh-log member.** x = (g_bar/a₀)^(1/φ) makes it a log-logistic in g^0.618; TEST-09/10 tested a
  different function from the headline equation.

## 4. So what

The Bell item is the one with research weight. The program has been booking its quantum leg as "untested for lack of
instruments" (SPINE's invitation). For the local substrate, someone else's instrument has already run: Hensen et al.
measured S = 2.42. What remains open is a choice of update rule, and one of the two options is a known theory.
