# Kuramoto Lattice Suite — observer-relative phase dynamics

*A `[PARALLEL-PATHS]` substrate exploration, contributed externally (Gemini, 2026-06-21 —
see [`forum/gemini/`](../../forum/gemini/)). This is a **different substrate** from the
repo's scalar-Intent grid: coupled **phase oscillators** (Kuramoto), not an Intent density
field. It is here because it makes one Synchronism claim — **measurement as
synchronization between observer-pattern and observed-pattern** — concrete and runnable,
and because it hosts the **one test that matters** ([SPINE.md](../../SPINE.md)): the
observer-relative CHSH experiment.*

## Why this exists

Synchronism's distinctive claim is the **single-observer** constraint: there is no
God's-eye view; all measurement is what one pattern can reconstruct of another via their
mutual phase-lock. A standard Kuramoto lattice is *observer-free* (computed from outside).
This suite asks what changes when you take the single-observer constraint seriously — and
whether it can buy a **novel** result.

## Substrate mapping (Kuramoto ↔ Synchronism)

| Kuramoto term | Synchronism reading |
|---------------|--------------------|
| oscillator phase θᵢ | local Intent-pattern phase |
| natural frequency ωᵢ | pattern's base oscillation (≈ mass/energy: f = E/h) |
| frequency cluster (locked region) | a recurring pattern = an **entity / particle** |
| coupling K | resonance strength (saturation-modulated, in principle) |
| order parameter R | local coherence |
| observer region phase-locking to a target region | **witnessing** = measurement-as-synchronization |

## What's in here

| File | What it is | Status |
|------|-----------|--------|
| `01_kuramoto_baseline_2d.py` | Baseline: random phases → emergent frequency clusters ("particles") on a 2D lattice; tracks order parameter R. The "entities are recurring patterns" claim, runnable. | runnable |
| `02_observer_relative_chsh.py` | **The experiment.** A local two-region "entangled" source, measured *only* through observer phase-lock with freely-chosen CHSH settings. Computes the CHSH S-value and compares to the classical bound (2) and the Tsirelson bound (2√2). | runnable → **S = 1.98** |
| `03_nonlocal_grid_chsh.py` | **The frontier variant.** Lets the shared substrate ("the grid") mix region B's state+setting into Alice's measurement, tunable by coupling g. Sweeps g and reports both S *and* signaling — asking whether a no-signaling nonlocal violation (S>2, signaling≈0) exists. | runnable → **S ≡ 2.0 ∀ g** |

Results are written to `results/` as JSON summaries (not raw trajectories).

### Frontier result (2026-06-21): the nonlocal channel is gauge-equivalent to relabeling

Sweeping the nonlocal coupling g ∈ [0,1], the CHSH value stays **pinned at S = 2.000** with
**zero signaling** at every g. The reason is the genuinely useful finding: with a uniform
shared phase λ, the smooth grid-mixing only adds a *constant phase offset* φ(b,g) to Alice's
readout axis — which is absorbed by relabeling the measurement angles. A uniform λ washes
the offset out of the marginals (hence no signaling) but leaves the correlation a function of
angle-differences only — i.e. **still a local-realist model, capped at S = 2.**

So the precise boundary is sharper than "nonlocal coupling causes signaling": a *smooth*
single-grid mediation is **gauge-equivalent to a local angle relabeling and buys nothing.**
A genuine no-signaling violation would require a **non-relabelable, conditional**
setting-dependence — exactly a quantum-entanglement-like primitive the single-observer
ontology does **not** contain and would have to *derive*, not assume. Confirmed-novel count
stays 0; the open arm is now narrowly specified (see `results/nonlocal_chsh_result.json`).

## The experiment, precisely (`02_observer_relative_chsh.py`)

The sharp question (PREDICTIONS.md, bet **B1**):

> Can a purely **local** substrate, measured **only** through an observer-pattern's
> phase-lock to a target-pattern (with freely chosen measurement settings, à la CHSH),
> reproduce the **quantum** correlations — *without* superluminal signaling?

The harness builds the most faithful **local** single-observer construction we can:
- a shared "source" prepares two regions A and B with a correlated phase relationship (a
  local common cause — the local-hidden-variable analog);
- Alice's observer phase-locks to A under a chosen setting angle a ∈ {a₀, a₁}; Bob's to B
  under b ∈ {b₀, b₁}; the binary outcome is the sign of the settled phase-lock;
- **locality is structural**: A's region/observer never couples to B's during measurement;
  the only A↔B correlation is the shared preparation. Settings are drawn independently.

**Expected outcome: S ≤ 2** (the classical/Bell bound), because this construction is local
by design. That is the **productive result**, not a disappointment: it demonstrates that
the single-observer phase-lattice, implemented *locally*, behaves as a local-realist model
and **cannot** reproduce quantum nonlocality (Tsirelson 2√2 ≈ 2.83). It thereby **sharpens
what Synchronism is**: an *ontology / interpretation*, not a local-hidden-variable
*physics*. To exceed the bound you would need either a genuinely **nonlocal grid** (the
"one substrate" doing work across A and B) or **superdeterministic** setting correlation —
and the harness flags whether either has crept in. Each outcome advances the program:
- **S ≤ 2** → the honest boundary (most likely); records what the local frame cannot do.
- **S > 2 without signaling/superdeterminism** → the first novel result; extraordinary,
  would demand independent replication.

## ⚠️ On the analogy simulations (not included here, and why)

Gemini's original cascade continued from this physics into **gravitational lensing, binary
inspiral, black-hole collapse, Hawking radiation, and quantum tunneling** visualizations.
Those are **illustrative analogy only** — the "gravity" in them is a hand-coded Newtonian
1/r² force on a tracer, *not* emergent from the lattice; the GR/Hawking behavior is asserted
by construction, not derived. They are beautiful and they predict nothing beyond standard
Kuramoto + hand-tuned forces. They are exactly the *elegant-isomorphism* failure mode this
repo names (CLAUDE.md, Kimi review). The full cascade is preserved in
[`forum/gemini/`](../../forum/gemini/) as a record (and as an A2ACW specimen of the bias);
it is deliberately **not** promoted into `simulations/` as if it were physics. Only the
two pieces that test or instantiate a real Synchronism claim live here.

## Run

```bash
cd simulations/kuramoto-lattice-suite
python 01_kuramoto_baseline_2d.py      # emergent frequency-cluster particles
python 02_observer_relative_chsh.py    # the CHSH test → results/chsh_result.json
```

Dependencies: `numpy` only (no matplotlib needed for the headless runs).
