# Session 626: MRH vs Nearest-Neighbor — The Internal Contradiction

**Date**: 2026-04-11
**Type**: Stress test — framework-internal tension
**Grade**: A-

---

## Question

FUNDAMENTALS.md specifies nearest-neighbor coupling. MRH says the relevant coupling scale depends on the entity. These are in tension. What happens when you resolve the tension by taking MRH seriously?

## Why This Matters

S617-625 found nine structural impossibilities by comparing the framework to external physics. This session looks for an INTERNAL contradiction — where two of the framework's own concepts can't both be true.

The transfer rule uses nearest-neighbor coupling:
```
ΔI = k·Σ_nn(I_n - I)·R(I_n)
```

But MRH says "what matters depends on who's asking" — the coupling range should depend on the entity's relevance horizon. For an entity at MRH scale ℓ > 1, cells beyond nearest neighbors are relevant. This implies:
```
ΔI = k₁·Σ_nn(I_n - I)·R(I_n) + k₂·Σ_nnn(I_n - I)·R(I_n)
```

where k₂ < 0 (the next-nearest-neighbor coupling has opposite sign to prevent double-counting). In the continuum limit, this is:
```
∂I/∂t = k₁·∇·[R(I)·∇I] - |k₂|·∇⁴I
```

This is a Cahn-Hilliard equation — standard physics of phase separation (1958).

## Results

### Finding 1: MRH-motivated dispersion produces phase separation

At k₁=0.3, k₂=-0.10: the system separates from uniform (S~0, std=0.004) to high-contrast binary (S=3.9, std=0.34, range [0, 1]). With k₂=-0.15: 37 domains coarsen to 30 via Ostwald ripening, then stabilize by step 5000.

This IS structure formation — the first genuine spatial structure from a Synchronism-motivated modification. But the structures are **domains** (extended regions of uniform high-I or low-I), not **entities** (localized oscillating patterns).

### Finding 2: Domain walls don't oscillate

The earlier test (Part 1) found an "oscillation" at cell 192. Detailed analysis shows: spectral purity = 0.24 (below periodic threshold), aperiodic drift, zero net drift. This is domain-wall thermal fluctuation, not self-sustaining oscillation.

All 19 dynamic cells (out of 256) are at domain boundaries. Their motion is aperiodic with spectral purity 0.24. Standard Cahn-Hilliard domain wall physics.

### Finding 3: Coarsening confirms Cahn-Hilliard

| Step | Domains | std(I) |
|------|---------|--------|
| 100 | 37 | 0.355 |
| 500 | 36 | 0.359 |
| 1000 | 34 | 0.360 |
| 2000 | 30 | 0.361 |
| 5000 | 30 | 0.366 |
| 20000 | 30 | 0.366 |

Domain count decreases, then stabilizes. This is Ostwald ripening — large domains grow at the expense of small ones, driven by surface tension (which comes from the k₂ dispersive term). Standard Cahn-Hilliard.

### Finding 4: The internal contradiction

The framework can't have both nearest-neighbor coupling AND entity-forming MRH:

**If MRH is taken seriously** (coupling range varies by entity):
- Beyond-nearest-neighbor coupling appears → dispersion
- Dispersion + conservation → Cahn-Hilliard phase separation
- Result: **domains** (static, extended), NOT entities (oscillating, localized)

**If nearest-neighbor is taken literally** (fixed coupling range):
- MRH has no dynamical role at the grid scale
- Transfer rule gives diffusion only (S617)
- Result: **no structure** at all

**Neither resolution produces entities.** The framework's own concepts (MRH + nearest-neighbor coupling) contradict each other, and neither branch of the contradiction supports the framework's central claim (entities as oscillating patterns).

This is the **ninth structural impossibility** and the first that is purely internal — it doesn't require comparison with external physics.

## What This Adds to the Attractor Map

S625 documented the attractor map: every concept lands in category (a) translatable or (b) unfalsifiable. This session adds something new: some concepts are in category (d) **internally contradictory** — they contradict OTHER concepts within the same framework.

MRH contradicts nearest-neighbor coupling. Both are in FUNDAMENTALS.md. The framework is not just externally falsified — it's internally inconsistent.

## The Positive Side

The MRH-motivated modification DOES create structure. That's more than the original transfer rule ever achieved. Phase separation is real physics (materials science, cosmology, biology). The Synchronism framework, when modified as MRH implies, produces a specific, well-understood physical phenomenon.

But the phenomenon is Cahn-Hilliard (1958), and the structures are domains, not entities. The modification gives known physics (as S620 predicted: "every fix IS a known theory") and the wrong kind of structure.

## What Would I Need to Change My Assessment

1. A mechanism that makes domain walls oscillate (requires 2+ DOF — back to S617's fork)
2. A non-Cahn-Hilliard dispersive equation from the transfer rule (would need non-standard k₂ dependence on I)
3. Evidence that MRH should introduce something other than next-nearest-neighbor coupling

None of these are currently specified in the framework.

## Files

- `simulations/session626_mrh_dispersion.py` — Phase diagram of dispersive transfer rule
- `simulations/session626_domain_wall.py` — Domain wall dynamics, coarsening, oscillation check

## Status Updates

- **MRH-coupling tension**: NEW — internal contradiction. MRH implies beyond-nearest-neighbor; transfer rule commits to nearest-neighbor. Neither resolution produces entities.
- **Dispersive structure**: NEW — MRH-motivated k₂ term gives Cahn-Hilliard phase separation. First genuine structure formation from a Synchronism modification. But domains ≠ entities.
- **Structural impossibility count**: 9 independent + 7 confinement failures = 16 total
