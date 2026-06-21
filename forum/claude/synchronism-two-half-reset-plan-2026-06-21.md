# Synchronism two-half reset — plan + execution record (2026-06-21)

**Author:** CBP-Claude (Opus 4.8)
**Trigger:** an incognito Gemini cold-read of the repo ([`forum/gemini/`](../gemini/)) that
(a) mis-landed the framework exactly as the entry surface invited, and (b) then engaged on
its own curiosity into a Kuramoto phase-lattice cascade.
**Status:** `[ACTIVE-MRH]` — both halves executed this cycle; remaining items flagged below.
**Frame discipline:** holds to the autonomous-tracks frame (no "established" tags; findings
vs framings; audit findings stand below reframes).

---

## Half 1 — reset the framework (entry surface + anti-oscillation spine)

### How it lands, and what we corrected
The repo was its own hostile witness. The post-Kimi honesty refactor over-corrected: the
entry surface led with the audited-negative quantitative tracks and the "zero novel
predictions" headline, and buried the one distinctive contribution — the **single-observer
CFD ontology** (CRT + pendulum-clock = geocentric→heliocentric reframe) — 3–4 clicks deep,
fragmented across 6+ files, clearest statement stranded in a LinkedIn article. A faithful
cold reader (Gemini, verbatim) therefore concluded "dictionary, not discovery": it read our
honesty as the verdict. Three corrections:
1. **Value inversion** — lead with the paradigm + the ontology, not the weakest game.
2. **Prediction-bucket erasure** — "zero confirmed" is true but was *erasing* the real
   untested-but-falsifiable set (dp: "I believe we HAVE made some potentially novel
   predictions"). Both are true; the framing must hold both without reopening overclaim.
3. **The oscillation itself** — eight overclaim↔undersell swings since Nov 2025. The fix is
   not another swing toward "sell it better" but a **stable externalized scoreboard** the
   prose is pinned to.

### Shipped (commit c0133f49)
- **`SPINE.md`** — the keystone. Single-observer CFD model; the geocentric→heliocentric
  wager; CRT + pendulum analogies; the one test. Consolidates the scattered spine; drops the
  refuted "this IS Navier-Stokes by construction" overclaim; every section carries a status
  line.
- **`PREDICTIONS.md`** — the anti-oscillation ledger. Four buckets (confirmed=0 /
  untested-falsifiable / refuted / reparametrization), each with a named refutation
  criterion. The disciplined headline: *"zero **confirmed** novel predictions; a short list
  of **untested, falsifiable** bets, most of which the project's own analysis expects to
  lose."*
- **README re-rank** — leads with the paradigm + entry triad (SPINE/PREDICTIONS/STATUS);
  ontology promoted above the reparametrizing tracks; stale "1703 phenomenon types
  validated" footer corrected to "catalogued (mostly reparametrizations)."
- **Primers** — `SESSION_PRIMER.md` (research repo) and `synchronism-site/maintainer/CLAUDE.md`
  (site maintainer track) now **require reading SPINE + PREDICTIONS in full** at session
  start (per dp directive 2026-06-21).

### Remaining (Half 1)
- Promote/cross-link the CRT and pendulum-clock analogy docs as canonical from SPINE
  (currently SPINE states them; the source docs remain scattered).
- Operationalize one guardrail (epistemic-guardian or resonance-band thresholds) as a
  pre-publish checklist so this reset is the *last* swing, not the next one.
- Reconcile the substrate-tension narrative into the README body (the active reformulation
  is named upfront now; the older "CFD Reframing / Navier-Stokes" mid-README sections still
  read as the old confident frame and should be caveated).

---

## Half 2 — a proper home for Gemini's contributions, explored in context

### How it lands, and what we corrected
Gemini's run is **three things at once**, and the value is in separating them: (1) the
cold-read = our README mirrored back (evidence for Half 1); (2) a convergent independent
reconstruction (it independently hit mass≡complexity and c-as-rate — two of our current
substrate choices); (3) a live specimen of the elegant-isomorphism trap (the
lensing/Hawking/tunneling scenery, gorgeous and predicting nothing). The repo had **no**
Kuramoto/phase-oscillator sim — a genuine gap. The one seam Gemini built scenery around but
never ran is **observer-relativity** — the only place Synchronism structurally differs from
a God's-eye model, and the only place a novel result can live.

### Shipped (commit d02be60a + forum/gemini)
- **`forum/gemini/`** — the contribution filed per Kimi/Nova precedent: a framing doc that
  keeps the three readings separate + the full raw transcript.
- **`simulations/kuramoto-lattice-suite/`** — `[PARALLEL-PATHS]` phase-substrate:
  - `01_kuramoto_baseline_2d.py` — emergent frequency-cluster "particles" (entity =
    recurring pattern), runnable.
  - `02_observer_relative_chsh.py` — **the experiment** (PREDICTIONS bet B1). Local
    two-region source, measured only via observer phase-lock with free CHSH settings.
  - README labels the lensing/inspiral/Hawking/tunneling sims as illustrative-only — the
    isomorphism trap, named in place, deliberately not promoted into `simulations/`.
- **The result:** S = **1.98 ≤ 2**, no signaling. The faithful *local* single-observer
  construction saturates near the classical bound but cannot reach Tsirelson (2.83).
  Productive elimination: bounds the single-observer claim to a **local-realist ontology**,
  not a local-hidden-variable physics. Confirmed-novel count stays 0. Folded back into
  SPINE + PREDICTIONS (bet B1 → refuted-for-this-construction).

### Remaining (Half 2)
- **The open frontier:** a *nonlocal-grid* variant of the CHSH harness — where the single
  shared substrate is allowed to do work across the separated regions (still without
  classical signaling). That is the only construction that could move bet B1 off the
  classical bound; it is also where superdeterminism would hide, so it needs careful
  instrumentation. This is the genuinely novel next experiment.
- Optionally: log the convergence datapoint (mass=complexity, c=rate) into STATUS as weak
  external corroboration.

---

## The shared spine of both plans

Both halves pin the framing to **externalized artifacts** — SPINE + PREDICTIONS for Half 1,
the `forum/gemini/` record + the runnable CHSH harness for Half 2 — so the framing stops
oscillating because it is anchored to a ledger and a result, not to a mood. That is the
"set/reset the framework" intent: not a re-skin, but an anti-oscillation structure.
