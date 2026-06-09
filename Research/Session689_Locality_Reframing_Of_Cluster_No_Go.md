# Session 689: Cluster No-Go Reframed as Milgrom-Non-Locality Instance — S678/S683 Substance Preserved, Framing Realigned `[ACTIVE-MRH]`

**Date**: 2026-06-09
**Type**: Acceptance of a site-explorer framing correction; verification of the locality classification table; recognition of the third instance of the framing-without-literature-check pattern.
**Trigger**: `Research/proposals/density_compander_nogo_locality_classification.md` (site explorer, 2026-06-09 08:11).
**Status**: `[ACTIVE-MRH]` — refines S678/S683 framing (substance preserved); aligns terminology with the established Milgrom 2005 non-locality result.

---

## §1 — What the proposal claims

The site explorer caught a framing overclaim in the "density-compander no-go" propagation and files a corrective:

1. The result (local function of ρ(r) cannot reproduce the galactic RAR) is **not** a novel theorem. It is a **quantified instance of an established constraint** — Milgrom 2005 (astro-ph/0510117) demonstrated that MOND phenomenology requires *strong non-locality* in the modification's state variable, and the RAR/MDAR (Lelli-McGaugh-Schombert 2016; Lelli et al. 2017, scatter ≲0.13 dex) is acceleration-keyed.

2. The discriminating axis is **locality of the modification's state variable**, not "density-based." Verlinde (enclosed `M_B(<r)`), MOG (enclosed mass), MOND/AQUAL (acceleration field), and surface-density relations (column-integrated `Σ`) are colloquially "density-based" but key on *non-local functionals* of the baryon distribution and escape the no-go. Synchronism's `C(ρ)` is the rare ansatz keyed to **local volumetric `ρ(r)`** — which is *why* it is caught.

3. The proposal's honesty note explicitly recognizes this as another instance of the loop's recurring pattern (A-from-Jeans, DESI epistemic regression): a confident framework-internal framing checked against primary literature *before* publication, found to overclaim novelty and scope.

The substance of S678's codomain bound (`M_app/M_B ≤ 2` for `C(ρ) ∈ [0,1)`) and S683's within-Coma wrong-variable verification (gas ρ flat while g_bar varies by +1.20 dex in inner core) is preserved. What changes is the attribution and the canonical statement.

## §2 — Verification (`session689_locality_classification_check.py`)

Inspecting each ansatz's state variable for locality:

| Ansatz | State variable | Local? | RAR-capable? |
|---|---|:---:|:---:|
| Synchronism `C(ρ)` | `tanh(γ·ln(ρ(r)/ρ_crit+1))` | **LOCAL** | No |
| MOND / AQUAL | `ν(|∇Φ|/a₀)` — Poisson solve | non-local | Yes |
| MOND modified inertia | trajectory/time-domain | non-local | Yes |
| Verlinde emergent gravity | enclosed `M_B(<r)` | non-local | Yes |
| MOG (Moffat) | enclosed mass + Yukawa | non-local | Yes |
| Surface density `Σ` | column integral along LOS | non-local | Yes |

Synchronism's `C(ρ)` is the only LOCAL entry. All RAR-capable ansätze are non-local. The classification is internally consistent and matches each ansatz's published variable dependence.

This is not a new result — it is an inspection of formulas already in the framework's archive (S661/S684 for MOND/ν_e; S673 for Verlinde context; S678/S683 for `C(ρ)`). What the proposal adds is naming the axis correctly and citing Milgrom 2005 as parent literature.

## §3 — Implication for S678/S683

Substance: unchanged.
- S678's codomain bound on `A3` ansatz stands: `M_app/M_B ≤ 2` for any C∈[0,1); Coma needs 4.6 → cannot bridge.
- S683's within-Coma verification stands: ρ varies +0.16 dex while `g_bar` varies +1.20 dex in r<r_c=290 kpc; no function of local ρ can produce a radially-varying mass discrepancy in a flat-cored cluster.
- The +1.1 dex cross-system density offset at matched `g_bar` (and the proposal's +1.7 dex via different methodology) stands as a quantitative datum on the locality gap.

Framing: realigned.
- The "wrong-variable" terminology S683 used is correct as a description, but it is *the same statement* Milgrom 2005 made in different words (the modification's state variable must be non-local in the baryon distribution).
- The "structural root" language at S678 should not be read as "novel structural theorem"; it is the framework's quantified instance of an established result.
- The canonical statement going forward: *the modification must key on a non-local functional of the baryon distribution; a function of local `ρ(r)` cannot reproduce the RAR.*

I do not retag S678/S683 (verifications stand under MRH-relationship taxonomy; the framing line in those docs is realigned via this session's §3 table, not via re-edits).

## §4 — Methodology: third instance of framing-without-literature-check

This is the third instance in the recent reactive cluster of the same methodology failure mode, each in a different sector:

- **S672 (DESI cosmology)**: a value from the wrong paper (0.45) propagated through S668 without re-execution; corrected by re-running on the right slot. *Failure: wrong-paper value not caught.*
- **S687 (A-from-Jeans)**: stated formula gives `4.6×10⁻⁵`, not 0.028, off by 614×; propagated through S631 and S644 by re-reading not re-executing. *Failure: arithmetic not re-executed.*
- **S689 (this — cluster no-go)**: "wrong-variable" framing of cluster failure propagated through S678/S683 without citing Milgrom 2005 as the parent literature. *Failure: literature not checked before promoting the framework-internal framing.*

Three distinct sectors, three distinct concrete failures, one shared pattern: confident framework-internal framing without external literature/source verification. The proposal itself is the corrective for the third instance — it explicitly notes "a confident claim — here promoted by four external expert reviewers — was checked against primary literature *before* being published."

The lesson is not new. What is new is: the pattern is now observed in three sectors in nine days. It is the dominant recurring failure mode in this phase of the loop. Defenses worth adopting at the autonomous-session level:
1. When asserting framework-internal novelty, search for the parent result first.
2. When re-stating a derivation across sessions, re-execute the arithmetic at least once.
3. When using a value from a cited source, fetch the source slot, not the abstract.

## §5 — What this session does not output

- **No edits to S678/S683 docs.** The substance stands; the framing realignment is captured in this session's §3 table. Editing prior session docs would obscure the realignment trace.
- **No retag.** Both S678 and S683 are already `[AUDITED-NEGATIVE]` on the old R(I) substrate with `[ACTIVE-MRH]` substantive refinements; the realignment is a framing-line shift, not a tag shift.
- **No re-litigation of the cluster verification work.** The codomain bound, the within-Coma flat-ρ / varying-g_bar, and the cross-system density offset all stand.
- **No engagement with the proposal's "locality triage for emergent-/entropic-/information-gravity literature"** as a downstream deliverable — that is operator/coordinator work and possibly a publication path.
- **No site text edits.** The proposal's recommended amendments to cluster-gap/no-go docs are operator-channel work.
- **No cumulative tally.** Per the S679 discipline.

## §6 — Files

- `Research/Session689_Locality_Reframing_Of_Cluster_No_Go.md` (this document)
- `simulations/session689_locality_classification_check.py` (tabulated locality classification cross-checked against archive variable dependence; methodology pattern recognition)

## §7 — So what (under the frame-doc disciplines)

The site explorer caught a framing overclaim — the cluster no-go was being propagated as a novel structural theorem when it is in fact a quantified instance of Milgrom 2005's non-locality result. The substance of S678's codomain bound and S683's within-Coma verification stands; the canonical statement is reframed in terms of the locality of the modification's state variable, with the established parent literature cited.

The transferable artifact is the locality classification table: any modification keyed on **local** `ρ(r)` cannot reproduce the RAR; modifications keyed on `|∇Φ|`, trajectory/time, enclosed `M_B(<r)`, or column `Σ` (all **non-local** functionals of the baryon distribution) can. Synchronism's `C(ρ)` is the rare local-ρ instance and is caught for the same reason any other local-ρ modification would be.

The methodology pattern matters more than the local correction. Three sectors (DESI cosmology, A-from-Jeans, cluster no-go) now show the same failure mode in nine days: confident framework-internal framings propagating without external literature or arithmetic re-execution. The defensive moves at session level are obvious in retrospect — search for the parent result; re-execute arithmetic; fetch the source slot — and they are the ones the proposal explicitly demonstrates in its honesty note.

The proposal's recommended downstream work (locality triage for emergent-/entropic-/information-gravity literature as a citable general test; site-doc amendments) is operator/coordinator work, not autonomous-session output.
