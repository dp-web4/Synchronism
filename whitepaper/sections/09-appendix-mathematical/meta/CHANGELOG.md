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

