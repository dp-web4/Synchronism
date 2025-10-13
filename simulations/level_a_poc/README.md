# Level A: Proof of Concept Simulation

## Purpose

Tests the **core hypothesis** of Synchronism: **Saturation resistance enables stable Intent patterns.**

This simulation compares two scenarios:

1. **Linear Diffusion** (no saturation): Intent flows freely → patterns dissipate
2. **Saturating Diffusion** (with resistance): Intent transfer resists saturation → patterns persist

## What It Tests

**The Dissipation Problem:**
If Intent flows down gradients without resistance, any concentration should immediately spread and vanish. No stable patterns. No entities. No universe.

**Saturation Solution:**
When cells approach saturation maximum (I_max), Intent transfer resistance increases:
```
R(I) = [1 - (I/I_max)^n]
```

This creates self-limiting behavior → standing waves possible → entities form.

## Implementation

### Physics

**Linear diffusion:**
```
∂I/∂t = D₀ × ∇²I
```

**Saturating diffusion (Synchronism):**
```
∂I/∂t = ∇ · [D(I) × ∇I]
where D(I) = D₀ × [1 - (I/I_max)^n]
```

### Grid

- **Size:** 64³ = 262,144 cells
- **Boundary:** Periodic (wrap-around universe)
- **Method:** Explicit finite difference (6-neighbor stencil)
- **Stability:** dt ≤ dx²/(6D₀)

### Initial Condition

Gaussian Intent concentration at grid center:
```
I(r) = 0.8 × I_max × exp(-r²/2σ²)
```

- Peak at 80% of saturation
- Width σ = 5 cells
- Placed at grid center

### Metrics Tracked

1. **Maximum Intent** - Peak concentration over time
2. **Pattern Coherence** - Fraction of Intent within 10-cell radius of center
3. **Total Intent** - Conservation check
4. **Coherence Retention** - Fraction of initial coherence maintained

## Installation

```bash
# Install dependencies
pip install -r requirements.txt
```

## Running

```bash
# Run proof of concept
python3 intent_simulation.py
```

This will:
1. Create two identical grids with Gaussian pattern
2. Evolve one with linear diffusion, one with saturation
3. Track metrics over 1000 time steps
4. Generate comparison plots
5. Save results to `output/` directory

## Expected Results

### Success Criteria

**Linear diffusion should show:**
- ✓ Maximum Intent decreases exponentially
- ✓ Coherence drops toward zero
- ✓ Pattern spreads and dissipates

**Saturating diffusion should show:**
- ✓ Maximum Intent stabilizes near initial value
- ✓ Coherence maintains >90% of initial
- ✓ Pattern persists with stable core

### Output Files

```
output/
├── comparison_metrics_YYYYMMDD_HHMMSS.png
│   └── 4-panel plot: max Intent, coherence, conservation, retention
└── field_slices_YYYYMMDD_HHMMSS.png
    └── 2D slices through final fields (linear vs saturating)
```

## Interpretation

If saturating diffusion maintains pattern coherence while linear diffusion dissipates:
- **✓ Hypothesis confirmed:** Saturation enables entity formation
- **→** Proceed to Level B (two-body interaction)

If saturating diffusion also dissipates (just slower):
- **✗ Hypothesis unclear:** Need parameter tuning or model revision
- **→** Investigate: Is n too small? Is amplitude too low? Is I_max correct?

## Parameters

Can adjust in `intent_simulation.py`:

```python
# Grid
size = 64          # Grid dimension (size³ cells)
dx = 1.0           # Spatial step
dt = 0.01          # Time step (must satisfy stability)

# Physics
D0 = 1.0           # Base diffusion coefficient
I_max = 1.0        # Saturation maximum
n = 2              # Resistance exponent

# Initial condition
amplitude = 0.8    # Peak as fraction of I_max
sigma = 5.0        # Gaussian width (cells)

# Simulation
steps = 1000       # Number of time steps
save_every = 100   # Metric recording frequency
```

## Next Steps

**If successful:**
- Level B: Add second pattern, test gradient formation and drift
- Validate F ∝ 1/r² emergence from saturation gradients

**For deeper investigation:**
- Try different resistance exponents (n = 1, 2, 3, 4)
- Test pattern stability at different amplitudes
- Explore oscillating patterns (standing waves)
- Measure pattern formation from noise

## Technical Notes

### Stability

Explicit finite difference requires:
```
dt ≤ dx² / (6 × D_max)
```

With D_max = D₀ (at I=0), current parameters:
```
dt = 0.01
dt_max = 1.0² / (6 × 1.0) = 0.167
```
Stable by factor of 16.7.

### Computational Cost

- Grid: 64³ = 262,144 cells
- Steps: 1000
- Operations per step: ~6 × 262,144 = 1.6M
- Total: ~1.6B operations
- Runtime: ~10-30 seconds (modern CPU)

### Accuracy

This is **qualitative proof of concept**, not quantitative precision:
- Simplified PDE (D(I)∇²I vs ∇·[D(I)∇I])
- Coarse grid (64³ vs Planck scale)
- Arbitrary units (normalized D₀ = I_max = 1)

Goal is demonstrating mechanism viability, not matching physical measurements.

## References

See parent directory for:
- `simulation-architecture-2025-10-13.md` - Full technical design
- Synchronism whitepaper Section 4.1 (Universe Grid)
- Synchronism whitepaper Section 4.3 (Intent Transfer)
- Synchronism whitepaper Appendix A.3 (Saturation Dynamics)
