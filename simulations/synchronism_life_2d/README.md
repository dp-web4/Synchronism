# Synchronism Life 2D

**Inspired by:** Conway's Game of Life and Claus Emmeche's artificial life research

## Purpose

Simplified 2D Intent dynamics for rapid exploration of stable patterns.

**Why 2D instead of 3D:**
- 10× faster computation (128² vs 64³ cells)
- Easier visualization
- Faster iteration on pattern discovery
- Direct comparison to Game of Life patterns

**What makes it "Life":**
- Like GoL: Discrete cells, discrete time, local rules
- Unlike GoL: Continuous Intent values (0 to I_max), saturation dynamics
- Goal: Discover stable configurations (the "gliders" and "oscillators" of Intent dynamics)

## Pattern Types

### Test Patterns Implemented

**Gaussian Blob:**
- Spherically symmetric concentration
- Tests: static stability
- Analogy: GoL "still life" patterns

**Ring:**
- Circular Intent distribution
- Tests: oscillations, breathing modes
- Analogy: Could become "blinker" equivalent

**Glider Test:**
- Asymmetric shape inspired by GoL glider
- Tests: propagation, motion
- Question: Can Intent patterns move?

**Random:**
- Noise field (low amplitude, scattered)
- Tests: self-organization
- Question: Do patterns emerge from chaos?

## Installation

```bash
pip install numpy scipy matplotlib
# Or just numpy if scipy/matplotlib unavailable
```

## Usage

### Test Single Pattern

```bash
# Test Gaussian stability
python3 intent_life_2d.py gaussian

# Test ring oscillation
python3 intent_life_2d.py ring

# Test glider propagation
python3 intent_life_2d.py glider_test

# Test self-organization from noise
python3 intent_life_2d.py random
```

### Custom Experiments

```python
from intent_life_2d import IntentLife2D, run_pattern_test

# Create custom initial condition
grid = IntentLife2D(size=128)
grid.add_pattern('gaussian', center=(40, 40), amplitude=0.95, sigma=5)
grid.add_pattern('gaussian', center=(90, 90), amplitude=0.95, sigma=5)

# Evolve and watch interaction
for step in range(1000):
    grid.step()
    if step % 100 == 0:
        patterns = grid.detect_patterns()
        print(f"Step {step}: {len(patterns)} patterns detected")
```

## Expected Results

### Gaussian Blob (Strong Saturation)

**Parameters:** amplitude=0.95, n=3, sigma=5

**Prediction:**
- Should stabilize (>90% retention)
- Might oscillate slightly (breathing mode)
- Core stays near saturation

**If stable:** ✓ "Still life" equivalent found

### Ring Pattern

**Parameters:** amplitude=0.95, radius=15, width=3

**Prediction:**
- Might oscillate (pulsate)
- Could break symmetry
- Test for standing wave formation

**If oscillating:** ✓ "Blinker" equivalent found

### Glider Test

**Parameters:** Asymmetric shape, amplitude=0.95

**Prediction:**
- If stable: becomes symmetric (loses asymmetry)
- If moving: drift along gradient
- Unlikely to propagate like GoL glider (different physics)

**If moving:** ✓ "Glider" equivalent found (huge discovery!)

### Random Initialization

**Parameters:** density=0.15, amplitude=0.7

**Prediction:**
- Most noise dissipates
- Some regions reach saturation
- Emergent patterns form spontaneously
- Final state: scattered stable blobs

**If patterns emerge:** ✓ Self-organization validated

## Comparison to Game of Life

| Property | Game of Life | Synchronism Life 2D |
|----------|-------------|---------------------|
| **Grid** | 2D | 2D |
| **Cell state** | Binary (alive/dead) | Continuous (0 to I_max) |
| **Update rule** | Count neighbors | Saturation diffusion |
| **Conservation** | No (cells born/die) | Yes (Intent conserved) |
| **Stable patterns** | Still lifes, oscillators, spaceships | Still blobs, oscillators(?), drifters(?) |
| **Complexity** | Turing-complete | Unknown (to be discovered) |

## Output Files

```
output/
├── {pattern}_metrics_{timestamp}.csv
│   └── Time series: total Intent, max Intent, pattern count
├── {pattern}_final_{timestamp}.npy
│   └── Final Intent field (can reload for analysis)
└── {pattern}_results_{timestamp}.png
    └── Visualization: field + metrics plots
```

## Physics Parameters

```python
D0 = 1.0      # Base diffusion coefficient
I_max = 1.0   # Saturation maximum
n = 3         # Resistance exponent (higher = sharper cutoff)
dt = 0.01     # Time step (satisfies stability)
dx = 1.0      # Spatial step
```

**Why n=3 instead of n=2:**
- Level A used n=2 → weak saturation → dissipation
- n=3 gives sharper resistance cutoff
- R(0.95)^3 = 1 - 0.95³ ≈ 0.14 → 86% resistance
- Should enable stability

## Computational Cost

- Grid: 128 × 128 = 16,384 cells
- Steps: 2000
- Runtime: ~5-10 seconds (single core)
- Memory: ~5 MB

Much faster than 3D Level A (64³ = 262K cells, ~30 sec)

## Questions to Answer

1. **Do stable patterns exist?**
   - Is there any configuration that maintains >90% coherence indefinitely?

2. **What do they look like?**
   - Size? Shape? Symmetry?
   - Are they similar to GoL patterns?

3. **Do they oscillate?**
   - Static or periodic?
   - Natural frequencies?

4. **Can they move?**
   - Propagating solutions?
   - "Gliders" in Intent field?

5. **Do patterns emerge from noise?**
   - Self-organization?
   - Characteristic scales?

6. **Are properties quantized?**
   - Do only specific sizes/frequencies appear?
   - Natural units?

## Next Steps

**After pattern discovery:**

**Level 2D-B: Interaction Studies**
- Two stable patterns
- Measure gradients between them
- Test collision outcomes

**Level 2D-C: Pattern Catalog**
- Systematic parameter sweep
- Build library of all stable configurations
- Classify by properties

**Level 2D-D: Evolution**
- Long-term runs (10K+ steps)
- Track pattern births, deaths, mergers
- Emergent complexity?

## Connection to Synchronism Theory

This is **exploratory artificial life research** applied to Synchronism's computational model.

**Not trying to prove theory correct** - trying to discover what patterns this rule system can create.

Game of Life didn't predict gliders - they were discovered through exploration. What will we discover in Intent dynamics?

## References

- Conway's Game of Life: https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life
- Emmeche, C. (1994). *The Garden in the Machine*
- Wolfram, S. (2002). *A New Kind of Science*
- Parent exploration: `explorations/game-of-life-parallels-2025-10-13.md`
