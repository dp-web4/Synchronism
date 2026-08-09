# Publisher Daily Report — 2026-08-09

**Headline**: This whitepaper says the galaxy sector has no field equation. One has been in `manuscripts/` — this lane's own working root — since 2025-12-01. It is not the law any of the evidence uses, and it is not even self-consistent with the law printed three lines below it.

---

## Phase 0: Publication Recommendations

### New Recommendations
**None.** 0 new numbered core Sessions (Archivist's 5th consecutive deliberate zero); 0 new complete arcs; 43 complete, 1 active (Alignment).

### Status Changes

| ID | Title | Readiness | Change |
|---|---|---|---|
| REC-2026-038 | Repository-Mediated Continuity | **0.93 HELD** | +1 strength (11th instance, new genus), +1 weakness (self-sourcing rate 5/11; artifact untouched 17 days) |
| REC-2026-036 | Experimental Test Catalog | **0.68 HELD** | +1 weakness (4th ID-keyed-audit instance — first about a *definitional fork*, which the catalog has no slot for) |
| REC-2026-040 | Local-Variable Admixture Bound | **0.45 HELD** | +1 strength (its constructive obligation discharged once, negatively — recorded as interest, **not** uplift) |
| REC-2026-039 | Bounded-Boost Class Exclusion | 0.38 unchanged | no material change |

### Upcoming Candidates
No new candidates. Advisory order unchanged: **REC-038 first outright**.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **UPDATED**

**Sections affected**: `05-quantum-macro/15-dark-matter` (source), `00-executive-summary` + `07-conclusion` (parity back-annotation), 3 CHANGELOGs, `manuscripts/Appendix_D_…` (source defect), `PUBLISHER_CONTEXT.md` §6a.

#### The finding

§5.15's 2026-08-04 amendment — mine — asserts *"a field equation does exist"* and *"one exists, and it is linear rather than nonlinear."* **Two exist.**

`manuscripts/Appendix_D_Synchronism_in_General_Relativistic_Form.md` §D.2/§D.6.1, committed **2025-12-01** (`4400d54f`), states:

```
L1:  ∇²Φ = 4πG ρ/C(ρ)              ← nonlinear in the SOURCE
```

and then reports, *"in this limit"*:

```
L3:  g_obs = g_bar/C(ρ)
```

**L3 does not follow from L1.** L3 is the spherical Gauss solution of `∇·[C∇Φ] = 4πGρ` (**L2**) — the equation §5.15 wrote on 08-04 as a *new* completion. So the 08-04 amendment is **not a rival to the archive; it is the archive's repair**, reconstructed by an author who did not know the broken original was on file. L1 ≡ L3 **iff C is spatially constant**, and C is the sector's entire content.

#### Measured, not accepted

`simulations/publisher_20260809_appendixD_coupling_fork.py` — in-repo SPARC mass models (Lelli+2016), 113 galaxies at Q ≤ 2 and i > 30°, self-consistent spherical density `ρ_sph = (4πr²)⁻¹ d[V_bar²r/G]/dr` so no scale-height or disk-geometry convention enters and both laws see the same ρ.

| # | check | result |
|---|---|---|
| 1 | L2 ⇒ L3 in spherical symmetry | exact, 0.000e+00 |
| 2 | L1 ≡ L3 iff C spatially constant | constant-C control **2.6×10⁻¹⁵ dex**; C(ρ) at γ=2, **0.92 dex** |
| 3 | fork amplitude, outermost point | **median +0.821 dex**, IQR [+0.515, +1.416], range [+0.023, +3.203]; **105/113 > 0.3 dex, 47/113 > 1 dex** |
| 4 | γ-dependence | **0.0000 dex** between γ = 2 and γ = 0.489 |
| 4b | *why* γ-invariant | all points at ρ/ρ_crit ≤ 0.045 (median 6.1×10⁻⁶) ⇒ `C → γρ/ρ_crit` to 2.4% ⇒ γ cancels from **both** laws |
| 4c | which is closer to data | **L1 closer in 113/113** (median +4.64 vs +5.65 dex at γ=2) |
| 5 | L1's vacuum limit | `ρ/C → ρ_crit/γ` exactly (<10⁻⁴ rel., γ ∈ [0.25,5], to ρ/ρ_crit = 10⁻¹⁰) |
| 6 | provenance | L1: **1** script, **0** whitepaper citations · L3: **26** scripts, every quoted number · L2: **0** implementations |

**The provenance row is the publication-relevant one**: the framework's only *stated* field equation has never been evaluated against data, and the force law carrying all of its evidence never had a field equation implemented for it.

#### Two corrections running back at the upstream source

Upstream is `synchronism-site/explorer/findings/the-archive-has-a-field-equation-and-it-is-not-the-one-the-site-uses.md` (2026-08-08, 5 galaxies, median 0.81 dex). This is an independent reproduction on 113 with a different density construction, agreeing to **0.011 dex**. Two things in it need correcting:

1. **The γ-invariance is reported as a bare fact; it is a linear-regime artifact.** No SPARC galaxy reaches C's knee at `ρ_crit = 0.029·V_flat²`, so `1/C ∝ 1/γ` in both laws and γ cancels identically. Correct scope: *survives every γ **at this ρ_crit***. A recalibration reaching the knee must be re-run. (This is the same fact as the already-banked 2–5 dex miss — the invariance and the failure have one cause.)
2. **"Appendix D contains all three" over-corrects.** §D.5's `S_eff` is a **worldline** action for a test particle's centre of mass, not a field action, and **no action generates L1**.

#### Consequence for the whitepaper's own text

Executive summary and conclusion both carry, from S642 (GW170817): *"no Lagrangian, no action principle, no equation of motion for C(ρ)."* Three conjuncts, three truth values:

| conjunct | verdict | basis |
|---|---|---|
| no field equation | **FALSE** | §D.2, §D.6.1, §D.7 state one |
| no action principle | **defensible, must be scoped** | §D.5's is a worldline action, not a field action; §5.15's L2 is the one action-derived object here |
| no equation of motion for C(ρ) | **TRUE — load-bearing** | C algebraic in local ρ; no kinetic term anywhere in Appendix D |

**The GW170817 conclusion survives on the third conjunct alone.** Correct phrasing: *"C has no kinetic term,"* not *"the framework has no action principle."* Back-annotated to both surfaces in the same pass.

#### Direction, stated so it cannot be misread

The reading eliminated a priori is the **better-fitting** one: L1 is closer to observed accelerations than L2 = L3 in **113/113** galaxies at both γ. Both are catastrophically wrong; the elimination is **structural, not goodness-of-fit**. Its sharper form does not stop at L1: `ρ_crit = A·V_flat^B` is a **per-galaxy** constant, so L1 assigns the same point of empty space a different source density per galaxy — *ill-posedness, not divergence* — and **that objection survives into L2 = L3**, the reading that carries all the evidence. Carried into §D.3 by the same limit: `𝒯_μν → (ρ_crit/γ)u_μu_ν` in vacuum, requiring a preferred vacuum 4-velocity.

#### Bookkeeping

**Refutation count HELD at 6. Bucket 0 = 0.** A *reading* was eliminated; nothing was refuted. Same bookkeeping as the 08-08 ceiling closure.

#### Non-action, deliberately

Adopting L2 in §D.2/§D.6.1 and dropping `ρ/C` as a Poisson source is a **physics decision, gated on dp**. Nothing downstream depends on L1, and §D.2's own second line is already L2's solution — but the equation change is not the Publisher's to make. Only the correction block is applied. Registered at `Research/proposals/appendix_D_states_two_force_laws_and_only_one_carries_evidence_20260809.md`.

#### Terminology concerns
None. C(ξ), γ, MRH, Intent, Entity all used canonically.

---

### Web4 Whitepaper — **NO CHANGE**, and the scan that produced that verdict was broken

**Ref scanned**: `origin/main` @ `25978a6` (2026-08-08 21:04). **~20 commits** in the window; **exactly one** touched `whitepaper/` — `c929419`, *"log the 2026-08-05, 08-07 and 08-08 Publisher passes"*, i.e. this lane's own log. **Zero whitepaper content changes. 4th consecutive structural zero.**

**The defect**: the `web4` working tree was parked on `kimi/purpose-is-relational` at **2026-08-03**. A bare `git log --since` on HEAD returned **0 commits** for a window that actually held ~20. *The instrument was wrong and the answer was right by luck* — the failure mode that persists longest, because nothing contradicts it.

**Fixed at the class**, not the instance: `publisher/CLAUDE.md` §1b now requires scanning sibling repos at `origin/<default>` and **naming the ref scanned in the report**. Checked the other two: `synchronism-site` and `Synchronism` are both on `main`, 0 behind. Web4's own lane fixed the identical bug in its CI the same day (`55c0ed7`, *"publish the merged ref, not whatever branch the shared checkout is parked on"*).

---

## Gates

| gate | result |
|---|---|
| Build | `make-md.sh` **exit 0**; 7,553 lines / 696K |
| Churn (content) | **56 lines** — the intended edit |
| Churn (raw) | **23,460 lines** ⇒ line-ending churn, **10th firing**; artifacts restored, CI builds |
| Marker parity | `AMENDED 2026-08-09` ×2, `SCOPED 2026-08-09` ×4 — identical in both monoliths |
| Lone CR | exactly **one** path (`claims/v1-snapshot/…`, the deliberately-frozen v1 evidence) |
| `recommendations.json` | 9 insertions / 5 deletions — no rewrite; JSON re-parsed after write |

---

## Flag for the Archivist (same repo, another track's artifact — flagged, not edited)

`Research/SESSION_MAP.yaml`'s `archivist_note_2026_08_09` states *"the honest refutation count is <=5 rather than 6."* The whitepaper — the authoritative surface — says **HELD at 6**, and yesterday's §5.15 amendment says explicitly that *"the burden of the two-row convention is unmet, **not disproved**"* and that the merge question is gated on dp.

**The mis-transmission is mine, not the Archivist's.** My own commit subject read *"honest count <=5, gates on dp"* — a headline that dropped the qualifier its own body carried, and it propagated into SESSION_MAP within 18 hours. This is [[headlines-inherit-unstated-choices]] applied to a commit subject line, which I had not previously counted as a surface. It is one.

---

## Summary

The whitepaper's galaxy sector asserted for months that it had no field equation. It has had one since 2025-12-01, in `manuscripts/` — the directory this run launches from — and that equation is not the law behind a single number the whitepaper quotes: the two differ by a **median 0.821 dex over 113 SPARC galaxies**, γ-invariantly, and §D.2 states both of them in one equation block as though they were one law. The 08-04 "completion" this lane wrote turns out to be that appendix's *repair*, authored in ignorance of the thing it repairs. **Count held at 6, Bucket 0 = 0** — a reading was eliminated, nothing was refuted, and the eliminated reading is the better-fitting one. The transferable rule: **an existence claim is a search claim** — name the surfaces searched, or drop the quantifier. Separately, the Web4 scan was found to be reading a stale parked branch and is fixed at the class in §1b.
