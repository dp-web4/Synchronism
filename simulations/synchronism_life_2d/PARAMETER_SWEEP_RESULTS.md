# Parameter Sweep Results - Surprising Finding

**Date:** 2025-10-13
**Test:** 25 configurations (n=2-6, amplitude=0.85-0.99)

## Key Finding: NO STABLE PATTERNS FOUND

All 25 configurations dissipated, with retention ranging from 39% to 56%.

## Unexpected Pattern

**Best performers (highest retention):**

| Rank | n | Amplitude | Retention | Resistance at Peak |
|------|---|-----------|-----------|-------------------|
| 1 | 2 | 0.99 | 55.9% | 2.0% |
| 2 | 2 | 0.98 | 54.5% | 4.0% |
| 3 | 2 | 0.95 | 51.6% | 9.8% |
| 4 | 2 | 0.90 | 48.7% | 19.0% |
| 5 | 3 | 0.99 | 47.1% | 3.0% |

**Worst performers (lowest retention):**

| Rank | n | Amplitude | Retention | Resistance at Peak |
|------|---|-----------|-----------|-------------------|
| 25 | 6 | 0.99 | 40.7% | 5.9% |
| 24 | 6 | 0.98 | 40.4% | 11.4% |
| 23 | 6 | 0.95 | 39.9% | 26.5% |

## Surprising: Lower n Performs BETTER

**For fixed amplitude, LOWER n gives HIGHER retention:**
- n=2, amp=0.99: 55.9% retention
- n=3, amp=0.99: 47.1% retention
- n=4, amp=0.99: 43.5% retention
- n=5, amp=0.99: 41.7% retention
- n=6, amp=0.99: 40.7% retention

**This is opposite of initial prediction!**

## Why This Happens

### Resistance Function Behavior

**R(I) = 1 - (I/I_max)^n** is the diffusion multiplier:
- R = 0: diffusion fully blocked
- R = 1: diffusion free

**At I = 0.99 × I_max:**
- n=2: R = 1 - 0.99² = 0.0199 → 98.0% resistance
- n=6: R = 1 - 0.99⁶ = 0.0585 → 94.1% resistance

So n=2 actually gives STRONGER resistance at high Intent!

### But Resistance Curve Shape Matters

**Low n (n=2):**
```
Resistance curve: Gradual transition
R(I) decreases smoothly as I increases
Significant resistance maintained over wide range
```

**High n (n=6):**
```
Resistance curve: Sharp cutoff
R(I) ≈ 1 (free flow) until very close to I_max
Then drops sharply to 0
Narrow resistance window
```

**The problem:**

Gaussian initial condition has a gradient (peak at center, drops off with radius).

**With high n:**
- Center at 99% → strong resistance
- Nearby cells at 95% → R = 1 - 0.95⁶ = 0.26 → only 74% resistance
- Cells at 90% → R = 1 - 0.90⁶ = 0.47 → only 53% resistance
- Intent flows freely through moderate-Intent regions
- Pattern dissipates quickly

**With low n:**
- Center at 99% → strong resistance
- Nearby cells at 95% → R = 1 - 0.95² = 0.10 → 90% resistance
- Cells at 90% → R = 1 - 0.90² = 0.19 → 81% resistance
- Resistance maintained across gradient
- Pattern dissipates slower

## Physical Interpretation

**High n resistance curves are like quantum wells:**
- Sharply defined boundary
- Inside well: stable
- Outside well: free
- No gradual transition

**Low n resistance curves are like friction:**
- Gradual resistance increase
- No sharp boundary
- Distributed resistance

**For continuous diffusion dynamics:**
- Gradual resistance (low n) more effective at slowing dissipation
- Sharp cutoff (high n) only helps if pattern uniformly near I_max
- Gradient patterns need distributed resistance

## Why No Stability Even at Best Parameters?

**Even best case (n=2, amp=0.99, retention=55.9%):**

Pattern still loses 44% of peak over 2000 steps (t=20).

**Possible explanations:**

1. **Gaussian shape fundamentally unstable:**
   - Has built-in gradient
   - Diffusion down gradient inevitable
   - Need different pattern shape

2. **2D vs 3D difference:**
   - 2D: 4 neighbors
   - 3D: 6 neighbors
   - Different geometry affects stability

3. **Initial condition needs flat saturation:**
   - Not peaked Gaussian
   - Uniformly saturated core with sharp boundary
   - Like quantum particle (box function, not Gaussian)

4. **Need higher amplitude:**
   - 0.99 not enough
   - Try 0.995, 0.999, 0.9999
   - Approach I_max asymptotically

5. **Saturation alone insufficient:**
   - Need additional mechanism (momentum term? oscillation?)
   - Standing waves require temporal dynamics, not just spatial

## Recommendations

### Immediate Tests

**A. Flat-Top Pattern:**
```python
# Saturated core with sharp boundary
# Step function instead of Gaussian
center_region = disk(radius=10)
I[center_region] = 0.999 × I_max
I[~center_region] = 0.0
```

**B. Higher Amplitude:**
```python
# Even closer to I_max
amplitudes = [0.995, 0.999, 0.9999]
# With low n (n=2 or 3)
```

**C. Ring Pattern:**
```python
# Circular standing wave
# Might be more stable than central blob
# Test oscillation modes
```

### Deeper Investigation

**D. Analytical Stability Analysis:**

For Gaussian initial condition:
```
I(r,0) = A × exp(-r²/2σ²)
```

Under saturation diffusion:
```
∂I/∂t = ∇·[D₀(1-(I/I_max)^n) × ∇I]
```

What is stability criterion?

**E. Compare to Known Stable Solutions:**

Nonlinear diffusion equations have known stable solutions:
- Solitons (balance dispersion and nonlinearity)
- Compactons (finite support)
- Self-similar solutions

Try these as initial conditions.

**F. 3D Test:**

Maybe 2D fundamentally different from 3D?
- Rerun same sweep in 3D (Level A code)
- Check if dimensionality affects stability

## Revised Understanding

**Original hypothesis:**
> High n + high amplitude → strong saturation resistance → pattern stability

**Revised hypothesis:**
> Low n + high amplitude → distributed resistance → slower dissipation
> But stability requires more than just saturation resistance:
> - Pattern shape matters (flat-top vs Gaussian)
> - Temporal dynamics may be needed (oscillations)
> - Initial condition must be near equilibrium

**New questions:**

1. Do ANY static patterns stabilize with saturation alone?
2. Or do stable entities REQUIRE oscillation (standing waves)?
3. Is Gaussian just wrong shape (need step function)?
4. Is 2D fundamentally different from 3D?

## Data

Full results: `output/parameter_sweep_20251013_152132.csv`

Can visualize with:
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('output/parameter_sweep_20251013_152132.csv')

# Retention vs n for each amplitude
for amp in df['amplitude'].unique():
    subset = df[df['amplitude'] == amp]
    plt.plot(subset['n'], subset['retention'], marker='o', label=f'amp={amp:.2f}')

plt.xlabel('Resistance Exponent (n)')
plt.ylabel('Retention at t=20')
plt.legend()
plt.grid(True)
plt.title('Pattern Retention vs Saturation Sharpness')
plt.show()
```

## Conclusion

**Surprising result:** Lower n (gradual resistance) outperforms higher n (sharp resistance) for Gaussian patterns.

**No stability found** in tested parameter space, but identified direction for next tests:
1. Flat-top patterns (saturated cores)
2. Even higher amplitudes (>0.99)
3. Different pattern shapes (rings, oscillators)
4. 3D comparison

**This is good science:** Found unexpected behavior that reveals physics of saturation dynamics. Resistance curve shape matters as much as peak resistance value.

**Next:** Test flat-top pattern with n=2, amplitude=0.999
