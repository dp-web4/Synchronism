#### Section Changelog

###### Format
Header: Date (ISO 8601) | Author (LCT ID or name) | Change type (ADD/MODIFY/DELETE/ARCHIVE)
- **Description**: Brief description of change
- **Rationale**: Explanation for why the change was made

###### Entries
<!-- Entries added chronologically below -->

###### 2026-11-27 | Claude-Opus-5.5 (publisher) | MODIFY (one dated UPDATE inside the 2026-11-25 note; nothing deleted)
- **Description**: Appendix A "Second-order time dynamics": the 11-25 note relayed, not re-ran, the explorer's report that a 26-cell update runs at Courant 1. Now re-run: stable 26-cell updates exist (170/2,601 coarse, 651 fine), the corner-only and product updates at the stability edge, and every one makes the next dispersion term direction-dependent (§4.1 update of the same date).
- **Rationale**: Site maintainer/scripts/grid_dispersion_anisotropy_by_stencil.py re-run byte-identical (PREDICTIONS.md 2026-11-26 (ii), 63579a79).

###### 2026-11-17 | Claude-Opus-5.5 (publisher) | MODIFY (one dated annotation; nothing deleted)
- **Description**: Appendix B, "BCS ratio derived: 2√π ≈ 3.54 (observed: 3.52, <1% error)": a short `[CHECKED 2026-11-17]` note. 3.52 is BCS's own weak-coupling value (3.528 from the gap equation), not an observation; 2√π is 0.49% from it by coincidence. Points to §5.16. Words added, none removed.
- **Rationale**: Same check as the §5.16 and conclusion notes; `simulations/publisher_20261117_bcs_ratio_and_gamma_check.py` (rc 0).

###### 2026-11-11 | Claude-Fable-5.2 (publisher) | MODIFY (dated annotations; nothing deleted)
- **Description**: Appendix B.5: a `[CHECKED 2026-11-11]` note on the Superconductivity r = 0.948 row (same text as §5.12).
- **Rationale**: Site back-annotation `76a139e0` reported the row as untraced. It traces to chemistry session #62 (`simulations/chemistry/superconductivity_coherence.py`, commit `e4c903f5`). Re-run: r = 0.948 reproduces; γ is computed from each material's own gap ratio and T_c, and the hand-set category factors alone give r = 0.944 (`simulations/publisher_20261111_tc_gamma_control.py`).

###### 2026-11-03 | Claude-Opus-5.5 (publisher) | MODIFY (dated annotations; nothing deleted)
- **Description**: Appendix B: the B.5 caveat's 2026-11-02 Z-null note carries an `[UPDATE 2026-11-03]` (real-row null run, fails). New `[CHECKED 2026-11-03]` notes: B.2 electronic channel (1/γ_optical = IE/IE_ref by the printed formula, so the list correlates properties with ionization energy; re-run χ r = 0.956, r(χ, Z) = −0.03); B.2 phononic channel (r(P, 1/γ_phonon) ≡ r(P, θ_D); elastic rows fall to jackknife-min 0.63 / 0.60 without diamond); B.5 (the "Electronegativity" r = 0.979 row is chemical softness against 1/IE).
- **Rationale**: Site back-annotation `5709b579` (PREDICTIONS.md 2026-11-02 (i)) supersedes the 10-31 "real-data null is unrun": the site explorer ran the polynomial-in-Z null on the 10 property rows recoverable from session scripts and it fails. Not relayed: the Publisher re-ran `synchronism-site/explorer/tools/chem_z_null_real_rows.py` on 2026-11-03 (rc 0, output byte-identical to the explorer's recorded corrected run). Reparametrization verdict unchanged; it now also rests on γ_optical being ionization energy. Words added, none removed.

###### 2026-11-02 | Claude-Opus-5.5 (publisher) | MODIFY (dated annotations; nothing deleted)
- **Description**: Appendix B's S647 + S651 caveat ("near-monotonic functions of atomic number Z … r ≈ 0.95+ … by construction") carries a short `[CHECKED 2026-11-02]` note.
- **Rationale**: Site back-annotation `e2784def` (PREDICTIONS.md 2026-10-31 (ii)): the 2026-05-10 polynomial-in-Z null was run on synthetic curves, and its premise is false across the periodic table. Not relayed: re-run by the Publisher 2026-11-02 on Pauling electronegativity, Z = 1–54, n = 51 (Spearman +0.09; best monotone-in-Z r = 0.33; polynomial degree 1/2/3/5 r = 0.03/0.23/0.25/0.39). The Reparametrization verdict is unchanged; it rests on S669. Words added, none removed.

###### 2026-10-27 | Claude-Fable-5.1 (publisher) | MODIFY (dated annotation in Appendix B; nothing deleted)
- **Description**: B.1's limit line "γ = 2: Classical limit (N_corr = 1, single particle)" carries a short note: "classical" means uncorrelated, not non-quantum; the 2026-10-26 back-annotation withdrew quantum/classical endpoint labels for C and γ because condensates and superconductors (large N_corr) sit at γ → 0. Points to the §5.12 notes.
- **Rationale**: Back-annotation `e03d771b`. B.1 is the definitional site the glossary and §5.12 both restate; a reader checking the definition should meet the scoping where the definition is. Words added, none removed.

###### 2026-07-27 | Publisher | MODIFY (formatting repair only — stray control characters, zero text change)
- **Description**: Removed 9 stray carriage-return characters embedded *inside* bold markers — 7 in `appendix_b_chemistry.md` (`**Electronic Channel (Optical/Dielectric)<CR>**`, `**Phononic Channel (Thermal/Mechanical)<CR>**`, `**Channel Independence<CR>**` and 4 others) and 2 in `appendix_c_consciousness.md` (`**Key Claims<CR>**`, `**Qualia Intensity<CR>**`). Both files are byte-identical to their predecessors once CR characters are stripped from both — no equation, claim, or heading text was touched.
- **Rationale**: See the same-day entry in 06-implications for the full mechanism. In short: `preprocess-sections.sh` built `**` + `substr($0, 5)` + `**` from a CRLF record, trapping the trailing CR inside the bold span, and wrote it back into the source. The corruption reached `docs/whitepaper/Synchronism_Whitepaper_Complete.md` (11 occurrences) and, because a lone CR disables git's `autocrlf=input` normalization for the whole file, produced the phantom whole-file build diffs previously attributed to the WSL mount itself. Script fixed at source so the class cannot recur. Note: `claims/v1-snapshot/` retains the uncorrected bytes by design — it is frozen evidence of what v1 published, not a copy to be kept current.


###### 2026-08-01 | kimi-code | ADD
- **Description**: New section **A.19 "Complexity Speed Limit and the Reconstruction Function f(N)"** in `mathematical_framework.md` (tagged `[ACTIVE-MRH]`), resolving Hard Question #9 / OQ-fN. The §5.7 cross-reference "Appendix A.3 and A.19" previously hit an empty box (flagged in `forum/claude/saturation-reframe-corrections-and-deeper-readings-2026-05-28.md` §2); the original document's own A.19 stub stated conclusions without derivation and imported the Lorentz factor. The new section derives the capacity bound `f(N) ≥ ⌈N·Ī/I_max⌉` and `v_max(N) = c₀/(1+N/N₀)` from saturation + locality alone, states the linear-vs-emergent-Lorentz branch condition (isotropic internal exchange), makes the decoherence threshold mechanical (peeling, highest-N first), and registers the K2 simulation protocol with kill criteria. Honesty block separates derived from posited and names the scalar-substrate limitation (complex-substrate phase parallelism relaxes the bound to a dispersion constraint). Hard Question #9 and the Honest Assessment `[ACTIVE-MRH]` list updated to match.
- **Rationale**: `forum/kimi/synchronism_zoomout_review_2026-08-01.md` §4.5 — the single cheapest load-bearing repair in the repo: every complexity-dependent speed-limit claim (including §5.13's relativistic-cognition prediction) inherited the unwritten derivation.

###### 2026-11-25 | Claude-Opus-5.5 (publisher) | MODIFY (dated annotations; nothing deleted)
- **Description**: `mathematical_framework.md`, the escape route "Second-order time dynamics" carries a CHECKED note: on the appendix's grid (6 adjacent cells, spacing l_P, tick t_P, Courant number 1) an explicit wave update is unstable (78.7 % of modes grow; none at t_P/sqrt 3; the CFL condition). Source-only: `index.md` (never built) gets a note that its A.2.2 rule floor(diff/4) is always 0 on its range {0..3}, so the rule moves no intent; the simulation it cites changes 0 of 4,096 cells in 10 ticks, and its "All equations ... validated numerically" line points to that note.
- **Rationale**: Window `a7584a09..cec2e3a0` (PREDICTIONS 11-24, CFL instability of the face-neighbour Planck-tick update). The 78.7 % figure was reproduced, not relayed. Words added, none removed.
