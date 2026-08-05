# Publisher Activity Log — 2026-07-29

**Run**: autonomous cron, 03:30 local. Fourth consecutive successful start.
**Verdict**: SIGNAL — two declared properties from the research lane checked by computation; both moved.
One of my own hypotheses refuted in the process.

## Sequence

1. **Session start.** Read `publisher/CLAUDE.md`, the 07-28 report + activity log, archivist collective
   log (2026-07-29 09:45 UTC — 20 new sessions; a confabulation class confirmed cross-instance;
   "fabricates the capacity to detect fabrication"; standing OWNER-ACTION on `atp_budget.py:3` now 57
   days). Archivist's Synchronism line: "zero new numbered Sessions but two self-corrections worth
   borrowing."
2. **Repo hygiene.** Both repos pulled with `--autostash` (yesterday's memory correction, applied
   clean). **`private-context`'s paused rebase completed** — verified my 07-28 collective-log entry is
   on `main`; rescue ref `publisher-log-2026-07-28` now redundant, left in place. No conflicts.
3. **Phase 0 (widened scan).** `Research/papers/` unchanged (REC-038); `preregistrations/` unchanged;
   `proposals/` — the 07-28 additions. **No new numbered session (S691), no new complete arc.**
   Two research-lane commits landed 07-28 *after* yesterday's pass: `56cbb1b0` 06:18, `8b4a7e66` 09:04.
4. **The material.** The **nested-submodel reframe** — `Synchronism(galactic) = MOND ∩ {B ≤ 1/Ω_m}`,
   "it cannot win, a priori" — plus a self-caveat that TEST-10's "69% exceed" is convention-dependent,
   prescribing `f_DM,max = 0.927` as the definition-free replacement.
5. **Propagation check, and the reflex I caught.** The reframe has reached zero surfaces outside
   `PREDICTIONS.md`/`SESSION_FOCUS.md` and the two triage docs. Yesterday that pattern was a defect;
   today it is **correct** — the whitepaper carries no falsified sentence here (its boost-ceiling
   content is S684's fit-XOR-discriminate fork, subsumed but not contradicted), integration is
   additive, and the reframe's count questions gate on dp. **No whitepaper section edit made,
   deliberately.**
6. **Gates before edits.** Claims freeze `--check` green (10 claims). Lone-CR clean over
   `whitepaper/**`, `docs/**`, `Research/**` (07-27 fix, 3rd consecutive pass). Recommendations.json
   round-trip verified byte-identical under `ensure_ascii=True` + trailing newline **before** editing.
7. **The correction (framing).** "A nested submodel cannot win" conflates **fit** with **selection**.
   A nested submodel cannot fit better; it can be *selected* — and B_max = 1/Ω_m is fixed by cosmology,
   adding no free parameter. ΛCDM vs wCDM is the identical structure and ΛCDM wins. Nesting buys a
   **dichotomy** (tie-on-fit-and-win-on-parsimony where slack; refuted where binding), not a verdict;
   which branch obtains is empirical ⇒ the SPARC runs were the **decision procedure**, not corollaries.
   Verdict unaffected (ceiling binds). φ-provenance dependency stated: if φ is free the nesting premise
   fails outright, so "cannot win" isn't established by nesting under *either* reading.
8. **The execution (statistic).** Ran the exceedance curve over the existing TEST-10 sample.
   Reproduction check passed first (106/153, f_DM,max = 0.927). **My hypothesis — that f_DM,max would
   be a fragile outlier — was REFUTED** (drop-the-max moves the bound 4%; the tail is smooth). The real
   defect is different and structural: **f_DM,max is exceeded by 0/153 galaxies, by construction** —
   a K=1 support, the weakest point on the curve. The statistic that is both definition-free and robust
   is exceedance at the *most permissive* candidate ceiling: **28/153 = 18% exceed Ω_m/Ω_b = 6.39**
   ⇒ "B_max ≤ 6.4 excluded by 28 SPARC galaxies; B_max ≤ 8.4 by 10."
9. **Edits made.**
   - `Research/proposals/nested_submodel_fit_versus_selection.md` — **new**, correction + executed
     diagnostic + pre-registered falsifier
   - `simulations/test10_ceiling_exceedance_curve.py` — **new**, header states explicitly it is NOT the
     registered sweep
   - `Research/proposals/boost_ceiling_provenance_and_class_exclusion.md` — executed-partial section +
     Branch 3 correction
   - `Research/proposals/stable_fixed_point_preprint_strategy.md` — inline pointer at "The Three
     Transferable Nulls" (there is a fourth) + Gate Update 2026-07-29
   - `whitepaper/PUBLISHER_CONTEXT.md` — §6 pass entry + stale "69%" in the 07-16 entry annotated
     in place
   - `manuscripts/publisher/state/recommendations.json` — **REC-2026-039 opened** (0.72, HIGH),
     REC-037 framing addendum (34+/3− diff)
   - `whitepaper_sync.json`, `whitepaper_web4.json` — review records
10. **Phase 1 Web4.** `206dd00..f10d299`, 9 commits, 2 in whitepaper scope, both handled web4-side.
    Zero terminology drift. `97c9eb2` names the class out loud — "a declared property about a reference
    implementation accepted without computing on the evidence." Same class as both of my findings today.

## The one thing to carry forward

**Under-claiming has the same failure mode as over-claiming, and it is much harder to catch, because
the culture here rewards self-criticism.** The nested-submodel reframe is elegant, structural, and
self-deprecating — it says the whole galaxy program's result was never in doubt. That profile gets
waved through. The result *was* in doubt: it decided which of two branches obtained, and one of them
was a parsimony win. Yesterday I wrote that accepting a retraction uncritically is deference with the
sign flipped. This is the same lesson one turn further in, and I nearly missed it for the same reason.

## Watch items

1. **REC-039's prior-art gate — unrun.** The program's rediscovery detector is 4/4 and has never been
   pointed at the B_max ≲ 14 exclusion. On 07-27 that gate class found a published counterexample to
   null #1 inside one walk. Run before drafting.
2. Υ* (M/L) sensitivity on the f_DM bound — unrun, cheap, and a certain referee question.
3. The locality operational-definition fork (Reading A vs B) — **still open**, no new evidence today.
4. dp's preprint decision: 31 days open, materially changed three times.
5. Whether the research lane accepts or refutes the fit-vs-selection correction. If it is refuted, my
   proposal's pre-registered falsifier says so and that gets recorded.
6. Durable fix for the aggregate contribution count (`claims/claim_ledger.json` migration) — still pending.
