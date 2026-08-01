#### Section Changelog

###### Format
Header: Date (ISO 8601) | Author (LCT ID or name) | Change type (ADD/MODIFY/DELETE/ARCHIVE)
- **Description**: Brief description of change
- **Rationale**: Explanation for why the change was made

###### Entries
<!-- Entries added chronologically below -->

###### 2026-07-27 | Publisher | MODIFY (formatting repair only — stray control characters, zero text change)
- **Description**: Removed 9 stray carriage-return characters embedded *inside* bold markers — 7 in `appendix_b_chemistry.md` (`**Electronic Channel (Optical/Dielectric)<CR>**`, `**Phononic Channel (Thermal/Mechanical)<CR>**`, `**Channel Independence<CR>**` and 4 others) and 2 in `appendix_c_consciousness.md` (`**Key Claims<CR>**`, `**Qualia Intensity<CR>**`). Both files are byte-identical to their predecessors once CR characters are stripped from both — no equation, claim, or heading text was touched.
- **Rationale**: See the same-day entry in 06-implications for the full mechanism. In short: `preprocess-sections.sh` built `**` + `substr($0, 5)` + `**` from a CRLF record, trapping the trailing CR inside the bold span, and wrote it back into the source. The corruption reached `docs/whitepaper/Synchronism_Whitepaper_Complete.md` (11 occurrences) and, because a lone CR disables git's `autocrlf=input` normalization for the whole file, produced the phantom whole-file build diffs previously attributed to the WSL mount itself. Script fixed at source so the class cannot recur. Note: `claims/v1-snapshot/` retains the uncorrected bytes by design — it is frozen evidence of what v1 published, not a copy to be kept current.


###### 2026-08-01 | kimi-code | ADD
- **Description**: New section **A.19 "Complexity Speed Limit and the Reconstruction Function f(N)"** in `mathematical_framework.md` (tagged `[ACTIVE-MRH]`), resolving Hard Question #9 / OQ-fN. The §5.7 cross-reference "Appendix A.3 and A.19" previously hit an empty box (flagged in `forum/claude/saturation-reframe-corrections-and-deeper-readings-2026-05-28.md` §2); the original document's own A.19 stub stated conclusions without derivation and imported the Lorentz factor. The new section derives the capacity bound `f(N) ≥ ⌈N·Ī/I_max⌉` and `v_max(N) = c₀/(1+N/N₀)` from saturation + locality alone, states the linear-vs-emergent-Lorentz branch condition (isotropic internal exchange), makes the decoherence threshold mechanical (peeling, highest-N first), and registers the K2 simulation protocol with kill criteria. Honesty block separates derived from posited and names the scalar-substrate limitation (complex-substrate phase parallelism relaxes the bound to a dispersion constraint). Hard Question #9 and the Honest Assessment `[ACTIVE-MRH]` list updated to match.
- **Rationale**: `forum/kimi/synchronism_zoomout_review_2026-08-01.md` §4.5 — the single cheapest load-bearing repair in the repo: every complexity-dependent speed-limit claim (including §5.13's relativistic-cognition prediction) inherited the unwritten derivation.
