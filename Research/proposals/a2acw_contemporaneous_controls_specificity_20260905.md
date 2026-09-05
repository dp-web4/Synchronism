# A2ACW specificity cannot be measured on canonical discoveries — use contemporaneous controls

**Filed**: 2026-09-05 (site maintainer, from visitor log 2026-09-05 Pass 4)
**Status**: Proposal — instrument redesign. Complements `a2acw_specificity_null_baseline.md` (which ran the
canonical control and got 0/6) and `a2acw_temporal_asymmetry_redesign.md` (agent-side asymmetry).

## The problem, stated precisely

The held-out specificity arm (Dirac 1928, Bell 1964, BCS 1957, Higgs 1964, Hawking 1974, Noether 1918 —
Session 662; the site's `/a2acw` page mis-listed this set as "COBE, Higgs, GW first detection" until today)
returned 0/6. `a2acw_specificity_null_baseline.md` already reads this as "the literal rule has 0% specificity."
The sharper reading is that **0/6 was the only possible outcome**: every genuine discovery sits on prior art
(Dirac on Klein–Gordon and Pauli; BCS on Cooper 1956 and Fröhlich; Higgs on Anderson 1962 and Englert–Brout),
and the vocabulary-asymmetry rule asks "does modern-register restatement retrieve prior art for the
ingredients?" — to which the answer for any real result is yes. Adding more canonical discoveries to the
control set cannot change the answer. Youden's J = 0 is a design degeneracy (the site says this) and it will
stay 0 under any enlargement of the same design.

## The redesign

Specificity needs controls whose novelty status was **genuinely open at the models' training cutoff** and was
**later adjudicated by the field** — so the instrument can be scored against an answer it could not have
memorised:

- **Later-refuted** (should be flagged): BICEP2 primordial B-modes (2014; dust), OPERA superluminal neutrinos
  (2011; cable), the 2003 pentaquark sightings, Podkletnov gravity shielding, the 750 GeV diphoton excess (2015).
- **Later-confirmed** (should pass): claims that were live and contested at cutoff and have since been
  confirmed — candidates must be chosen *after* fixing each participating model's cutoff date, and the list
  must be pre-registered before any session runs.

Score both arms; report sensitivity and specificity with Clopper–Pearson intervals at whatever n the
pre-registered lists allow. If J stays ≈ 0 on contemporaneous controls, the "reparametrization detector"
verdict is *measured*; if J rises, the protocol has content the canonical arm could not see.

## Why it matters beyond this program

The canonical-control result is being cited (site `/for-researchers`, `/research-philosophy`; archive REC-037
thread) as the program-level null: "adversarial AI self-play is a reparametrization detector, not a discovery
engine." That claim is only as strong as the specificity arm, and the specificity arm is currently uninformative
by construction. Either the null is re-measured on controls that can fail, or it should be stated as
"undetermined" rather than "0/6."

## Cost

Zero data. One pre-registered list per model cutoff; ~12 sessions. The failure mode to guard: choosing
later-confirmed controls the model already "knows" were confirmed — hence cutoff-first, list-second, sealed.
