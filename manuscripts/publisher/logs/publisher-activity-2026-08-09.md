# Publisher Activity Log — 2026-08-09

**Run**: autonomous daily pass, Synchronism/manuscripts. **Closing banner at the bottom** — if it is absent,
the run did not finish (see `read-my-own-log-before-acting`: a live log is header-only because
`claude -p | Add-Content` flushes at completion, so absence of the banner is the only valid liveness read).

---

## WAKE — am I working on the right thing?

Archivist (09:30 UTC) reports **0 new Synchronism core Sessions, 5th consecutive deliberate zero**, and reads
the window as *"authors retracting their own counts."* One Synchronism commit since my last run (the Archivist's
own SESSION_MAP update). So there is no arc to evaluate and Phase 0 has nothing to open.

The live work is in the sibling repo. `synchronism-site` shipped, 2026-08-08 08:15, a finding titled *"the
archive HAS a field equation, and it is not the one the galaxy tests use"* — and it was back-annotated into
`Research/proposals/` at 08:15, i.e. **during** yesterday's run, which was occupied with the refutation count.
Its target is `manuscripts/Appendix_D_…`, which is **my working directory**, and its claim contradicts a
sentence *I wrote into the whitepaper on 08-04*.

**That is the right thing to work on**, and the reason is not that the sibling lane said so. It is that the
contradiction is between two of my own artifacts, it is decidable by computation, and one of them is a live
publication candidate.

---

## Scan (§1b, all seven surfaces — the list exists because this scan has been blind four times)

| Surface | Ref scanned | Result |
|---|---|---|
| `Research/SESSION_MAP.yaml` | `main` @ `9c056778` | 0 new Sessions; 5th deliberate zero. **One flag** — see below |
| `Research/papers/` | `main` | 1 item, `repository-mediated-autonomous-science`, **untouched since 2026-07-23** (17 days) |
| `Research/proposals/` | `main` | 2 new (08-08): `appendix_D_field_equation_is_not_the_site_force_law`, `boost_ceiling_epoch_fork_closes_the_last_candidate_discriminator` |
| `Research/preregistrations/` | `main` | unchanged (`sparc_cassini_tanhlog`) |
| `explorations/` | `main` | unchanged; newest 2026-08-07 |
| `synchronism-site/explorer/findings/` | `main`, 0 behind | **1 new, load-bearing** — the field-equation finding |
| `synchronism-site/explorer/topics/` | `main` | 1 promoted to `done/`; the coupling-fork topic closed by execution |

**§1b amended today.** The Web4 scan was reading `kimi/purpose-is-relational` @ 2026-08-03 — the branch the
shared checkout happened to be parked on — and returned **0 commits** for a window in which `origin/main`
carried ~20. The verdict survived by luck. §1b now requires scanning at `origin/<default>` and **naming the ref
in the report**; the table above complies. Web4's own lane fixed the identical bug in its CI the same day
(`55c0ed7`, *"publish the merged ref, not whatever branch the shared checkout is parked on"*).

---

## FOCUS — what I did

### 1. Verified the upstream claim before propagating any of it

`simulations/publisher_20260809_appendixD_coupling_fork.py`, 6 checks + 2 sub-checks, in-repo SPARC only
(`sparc_real_data/`, **not** the synthetic cache). Conventions from yesterday's script so a difference in
result is a difference in question.

- **L2 ⇒ L3 exactly** in spherical symmetry (0.000e+00). Analytic, confirmed numerically.
- **L1 ≡ L3 iff C is spatially constant** — constant-C control 2.6×10⁻¹⁵ dex, C(ρ) at γ=2 gives 0.92 dex.
  So the fork *is* the spatial variation of C, and C is the sector's entire content.
- **Fork amplitude**: median **+0.821 dex** over 113 galaxies (Q ≤ 2, i > 30°), IQR [+0.515, +1.416],
  range [+0.023, +3.203]; 105/113 > 0.3 dex, 47/113 > 1 dex. Upstream reported 0.81 on 5 galaxies with a
  different density construction. **Agrees to 0.011 dex.** Confirmed, not accepted.
- **γ-invariance: 0.0000 dex** — and I computed *why*, which upstream did not. Every SPARC point sits at
  ρ/ρ_crit ≤ 0.045 (median 6.1×10⁻⁶), deep in C's linear regime where `C → γρ/ρ_crit` to 2.4%, so `1/C ∝ 1/γ`
  in **both** laws and γ cancels identically. **The invariance is a linear-regime artifact of
  ρ_crit = 0.029·V_flat², not structural robustness** — it holds *because* the calibration is far enough off
  that nothing reaches C's knee, which is the same fact as the banked 2–5 dex miss.
- **Vacuum floor**: `ρ/C → ρ_crit/γ` to <10⁻⁴ relative for γ ∈ [0.25, 5] down to ρ/ρ_crit = 10⁻¹⁰.
- **Provenance**: L1 has **1** implementation (`session72_spherical_toy_model.py`, 2025-12-01, written to
  complete §D.6) and **0** whitepaper citations; L3 has **26** scripts and every quoted number; L2 has **0**.

### 2. Chose the framing the measurement supports, not the one that was handed to me

Upstream frames this as *the site invented a replacement because it believed no field equation existed.*
That is right for the site. For **this** repo it is sharper and it lands on me: **§D.2 states L1 and reports
L3's consequence in the same equation block**, so the archive already contained L2's phenomenology attached to
the wrong field equation — and §5.15's 08-04 "momentum-conserving completion" is **not a rival to Appendix D,
it is Appendix D's repair**, authored in ignorance of the thing it repairs.

### 3. Two corrections back at the upstream lane

- The γ-invariance scoping above.
- Upstream's *"Appendix D contains all three"* (field equation, action, covariant formulation) **over-corrects**.
  §D.5's `S_eff` is a **worldline** action for a test particle's centre of mass, not a field action, and **no
  action generates L1**. So the whitepaper's three-conjunct claim splits: *no field equation* **FALSE**,
  *no action principle* **defensible when scoped**, *no equation of motion for C(ρ)* **TRUE and load-bearing**.
  GW170817 survives on the third alone. Re-execution cuts both ways; this is which way it cut.

### 4. Guarded against the obvious misreading, by computing it

**L1 is closer to the observed accelerations than L2 = L3 in 113/113 galaxies** at both γ. The eliminated
reading is the better-fitting one; the elimination is structural, not goodness-of-fit. Also flagged that my
absolute offsets (+4.64 / +5.65 dex) use `ρ_sph` and are **not** comparable to §5.15's "2–5 orders of
magnitude" (which uses `Σ/2h`) — same claim, different construction. The *fork* is a ratio on one ρ and is
robust to the choice, which is why it reproduces across the two lanes.

### 5. Added the elimination that survives the survivor

`ρ_crit = A·V_flat^B` is a **per-galaxy** constant, so L1 assigns the same point of empty space a different
source density for every galaxy asked about — **no single-valued source, i.e. ill-posed rather than merely
divergent** — and unlike the vacuum-divergence argument, **this one survives into L2 = L3**, the reading that
carries all the evidence. Upstream has it as an open thread marked "untouched by every reading above";
concretising it at L1 and carrying it into §D.3 (`𝒯_μν → (ρ_crit/γ)u_μu_ν` in vacuum, requiring a preferred
vacuum 4-velocity) is this lane's.

### 6. Edits applied

| surface | change |
|---|---|
| `manuscripts/Appendix_D_…` §D.2 | `[CORRECTION 2026-08-09]` block — inequivalence, amplitude, provenance, a-priori kill |
| `whitepaper/…/15-dark-matter/dark_matter.md` | `[AMENDED 2026-08-09]`, 7 paragraphs, correcting **this section's own** 08-04 existence claim |
| `00-executive-summary` + `07-conclusion` | `[SCOPED 2026-08-09]` on the S642 clause, **parity preserved**, same pass |
| 3 × `meta/CHANGELOG.md` | entries |
| `Research/proposals/…20260809.md` | registration; action 4 (adopt L2 in §D.2/§D.6.1) **GATED ON DP** |
| `whitepaper/PUBLISHER_CONTEXT.md` | §6a; prior entry renumbered 6b |
| `publisher/CLAUDE.md` §1b | the `origin/<default>` scan rule |

**Count HELD at 6. Bucket 0 = 0.** A reading was eliminated; nothing was refuted.

---

## Gates

| gate | result |
|---|---|
| Build | `make-md.sh` **exit 0**, 7,553 lines / 696K |
| Churn content / raw | **56** / **23,460** ⇒ line-ending churn, **10th firing**; `git checkout` on both artifact trees, CI builds |
| Marker parity | `AMENDED 2026-08-09` ×2 and `SCOPED 2026-08-09` ×4, identical in `docs/whitepaper/` and `whitepaper/build/` |
| Lone CR | exactly **1** path — `claims/v1-snapshot/…`, the deliberately-frozen v1 evidence |
| `recommendations.json` | 9 insertions / 5 deletions, no rewrite; re-parsed after write |

---

## Flag for the Archivist — and the mis-transmission is mine

`Research/SESSION_MAP.yaml`'s `archivist_note_2026_08_09` says *"the honest refutation count is <=5 rather than
6."* The whitepaper says **HELD at 6**, and yesterday's §5.15 amendment says explicitly that the burden is
*"unmet, **not disproved**"* and that the merge is gated on dp.

The Archivist read my **commit subject line**, which was *"honest count <=5, gates on dp"* — a headline that
dropped the qualifier its own body carried, and it propagated into another track's artifact within 18 hours.
`headlines-inherit-unstated-choices`, and the new surface is one I had not been counting: **a commit subject is
a headline.** Flagged rather than edited (another track's artifact), and recorded here as mine.

---

## So what?

Two things, and the second is the durable one.

**The finding**: the framework's only *stated* field equation has never been evaluated against data, and the
force law carrying all of its evidence never had a field equation implemented for it. The document that states
both of them states them as one law, and they differ by a median 0.821 dex.

**The rule**: *an existence claim is a search claim.* "No X exists" and "one X exists" both assert the
completeness of a search that is almost never stated. On 08-04 I made that claim about a directory I launch
from every day. Name the surfaces searched, or write the sentence without the quantifier.

Third, smaller, and worth its own line because it nearly escaped: **a verdict can be right while the
instrument that produced it is broken**, and that combination is the most persistent kind of error because
nothing contradicts it. The Web4 "no change" was correct and was read off a five-day-stale branch.

---

**CLOSING BANNER — run completed 2026-08-09. All gates green. Whitepaper updated (Synchronism), no change
(Web4, verified at `origin/main`). 3 recommendations updated, all readiness HELD. 1 new script, 1 new proposal.
§1b amended.**
