#### Section Changelog

###### Format
Header: Date (ISO 8601) | Author (LCT ID or name) | Change type (ADD/MODIFY/DELETE/ARCHIVE)
- **Description**: Brief description of change
- **Rationale**: Explanation for why the change was made

###### Entries
<!-- Entries added chronologically below -->

###### 2026-07-27 | Publisher | MODIFY (formatting repair only — stray control characters, zero text change)
- **Description**: Removed 2 stray carriage-return characters embedded *inside* bold markers in `04-open-questions/open_questions.md` — the published text read `**Post-Kimi consolidated open questions (2026-05-28)<CR>**` and `**Legacy categorical question lists<CR>**`. No word, sentence, claim or heading was altered: the file is byte-identical to its predecessor once CR characters are stripped from both. The same repair was applied to two files in 09-appendix-mathematical (11 occurrences repo-wide) and to the build script that produced them.
- **Rationale**: `preprocess-sections.sh` converts `### Header` to `**Header**` via `substr($0, 5)`. On a CRLF working tree the record's trailing CR is inside that substring, so it landed between the header text and the closing `**`, and the script wrote the result back into the *source* file. The stray CR then reached every published surface — the monolithic markdown at `docs/whitepaper/` carried all 11. It also had a second-order effect: a lone CR makes git classify a file as non-text, which silently disables `core.autocrlf=input` normalization for that file, so its CRLF endings began diffing against the LF blobs CI produces. That is the true cause of the '3,177-line content-free CRLF diff' recorded in the 2026-07-26 entries of the 00-executive-summary and 07-conclusion changelogs; that diagnosis named the symptom, not the mechanism, and is superseded here. The script is fixed at source (`sub(/\r$/, "", header_text)`, verified idempotent against a synthetic CRLF fixture), so this class cannot recur.

