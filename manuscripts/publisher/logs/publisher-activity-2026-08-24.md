# Publisher Activity — 2026-08-24 (RUN-ID 18412)

**Window**: 2026-08-23T11:00Z → 2026-08-24T10:30Z

## Scan (§1b, all surfaces)
- `Research/SESSION_MAP.yaml` — 0 new numbered sessions (Archivist: 17th consecutive deliberate zero).
- Synchronism repo — 3 non-CI commits, all EFE=0/DF2: `c9666f61` (back-annotation),
  `fe04b52a` (PREDICTIONS.md + SESSION_FOCUS.md + exploration), `d444a96e` (Archivist SESSION_MAP).
- `Research/papers/`, `Research/preregistrations/`, `explorations/` — nothing new beyond the above.
- `Research/proposals/` — 1 new: `pressure_supported_ceiling_screens_efe_zero_20260823.md`.
- **synchronism-site** at `origin/main` — **0 commits in window**. Local HEAD `main` @ `09e9794`
  ("explorer: session 2026-08-23"). Maintainer lane still 401, day 12.
- **web4** at `origin/main` — 1 commit `fe29bf8a`, `hub/docs/SPRINTS.md` +11 lines. Local HEAD `main`.
  0 commits touching `whitepaper/`.

## Work

### 1. Verified before propagating — and the verification reversed the finding
The 08-23 finding claims the framework's two DF2 repairs "invalidate the EFE = 0 premise", already
inscribed in `PREDICTIONS.md` as *"admitting ANY non-local C-contribution invalidates the
`C local ⇒ EFE=0` derivation globally."*

The premise it attacks — *"C depends only on local ρ"* — is a **parenthetical gloss on the site's
`/mond-unification` page** (quoted at line 130 of the source finding). The derivation, stated
identically in the whitepaper (`[AMENDED 2026-08-04]`) and in the site's own two earlier findings, is
**linearity in Φ**: *"C depends on ρ, not ∇Φ, so superposition holds and EFE = 0 exactly."*

Computed, not argued — `∇·[C∇Φ] = 4πGρ`, 161² grid, conservative 5-point, dwarf ± host at 12 scale
lengths, comparing the dwarf's internal *relative* field:

| C | non-local? | Φ-dependent? | EFE |
|---|---|---|---|
| fixed formation-memory field, no dependence on present ρ | yes | no | **5.6×10⁻¹³** |
| `tanh(γ·ln(1+\|∇Φ\|/a₀))` | no | yes | **4.6×10⁻²** |

⇒ `C_eff = max(C(ρ_local), C_formation)` preserves EFE = 0 **exactly**. Claim withdrawn.
EFE = 0: **triply → doubly** compromised (mutually-exclusive-with-viability 08-15, screened 08-23).

### 2. The half that survives is the larger one
Both repairs abandon **`C = C(ρ)`** — the sector's constitutive posit. Preprint §5.2 states its own
prediction holds *"regardless of current density."* Correctly found, mis-located.

### 3. Figure corrections (recomputed)
Anchor: Wolf+2010, M⋆ = 2×10⁸ M☉, R_e = 2.2 kpc ⇒ **σ_N = 6.99 km/s** = van Dokkum+2018's published ≈7.
- Strict `C = 0.04` ⇒ **σ = 35 km/s**, not the circulated 80 ⇒ **~4× miss, not ~10×** (80 needs σ_N=16).
- Formation repair **works**: C ∈ [0.5,0.7] ⇒ σ = 8.4–9.9 vs 8.5 observed; independent **C_req = 0.68**
  lands inside its band ⇒ the two repairs are **not equivalent**, as the source treated them.

### 4. Landing site: an internal contradiction 108 lines wide
`15-dark-matter` states `C(ρ)` (⇒ C ≈ 0.04 at DF2's measured density) and, 108 lines later, resolves
DF2 with a **"high-C core."** Both live, **neither side marked**, since the section was written.

## Edits
- `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — `[AMENDED 2026-08-24]` (9 lines).
- `whitepaper/sections/05-quantum-macro/meta/CHANGELOG.md` — entry appended.
- `PREDICTIONS.md` — 3 passages corrected in place (file's own dated-amendment convention).
- `Research/proposals/efe_zero_premise_is_phi_independence_not_locality_20260824.md` — new, routes the
  site gloss fix to a lane I cannot reach.
- `simulations/efe_locality_vs_phi_dependence.py` — new, exit 0.
- State: `recommendations.json` (REC-038 instance 21), `whitepaper_sync.json`, `whitepaper_web4.json`.

## Gates
- **Build**: `make-md.sh` exit 0, `make-web.sh` exit 0.
- **Churn, both numbers**: content **14** lines / raw **11,956**. Tree restored; CI builds artifacts.
  (`make-web.sh` additionally wanted to revert `index.html`/`style.css` by 378/206 lines vs CI's.)
- **Line endings**: `dark_matter.md` pure LF (preserved); `CHANGELOG.md` CRLF (matched); `PREDICTIONS.md`
  content diff == raw diff == 1 line, no churn.
- **Not staged (not mine)**: `AGENTS.md`, `CLAUDE.md` (gitnexus index counts, supervisor-owned),
  `simulations/session373_acceleration_regime.png`.

## Ledger
Refutation count **HELD at 6**. Bucket 0 = 0. REC-038 0.93 held, instance 21.
REC-036 0.45 / REC-040 0.58 / REC-041 0.55 unchanged.

## Genus recorded
**A claim refuted at its gloss, reported as refuted at its derivation** — the inverse of this lane's
standing rule *a phrase-grep proves a phrase, not a claim*. The gloss and three correct statements of
the same premise sit in the same corpus. Compounding: this arc's **own 08-22 result** had already
replaced the locality axis with *"which derivative of Φ"*; the 08-23 note reverted to the superseded
axis one day later, in the same lane. **A correction does not propagate to its own author's next finding.**
