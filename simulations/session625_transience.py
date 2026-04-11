"""
Session 625 Part 2: Is S624's Class 4 transient?

S624 reported Class 4 (edge-of-chaos) at (A=1.0, k=0.40).
But the coherence test shows ξ=1.0 and Δf≈0 after 2500 steps.
Is the "Class 4" just a long transient?

This matters because:
- If transient: the CA has NO sustained non-trivial regime. It always
  decays to trivial (Class 1-2). S624's finding is an artifact of
  short measurement windows.
- If sustained: the coherence test missed something. Check with
  the exact S624 setup.

Also: the deeper finding — spatial structure (ξ) and temporal dynamics
(Δf) are ANTI-CORRELATED. One field can carry one or the other,
not both. This is a new structural impossibility for entity formation.
"""

import numpy as np

I_MAX = 1.0


def R_nonmonotonic(I, n=2, A=0.5):
    x = np.clip(I, 0, I_MAX) / I_MAX
    base = 1.0 - x**n
    bump = 1.0 + A * np.sin(np.pi * x)
    return np.maximum(0, base * bump)


def step_1d(I, k, A=1.0):
    I_left = np.roll(I, 1)
    I_right = np.roll(I, -1)
    R_left = R_nonmonotonic(I_left, A=A)
    R_right = R_nonmonotonic(I_right, A=A)
    dI = k * ((I_left - I) * R_left + (I_right - I) * R_right)
    return np.clip(I + dI, 0.0, I_MAX)


def spatial_entropy(I, n_bins=20):
    hist, _ = np.histogram(I, bins=n_bins, range=(0, I_MAX))
    p = hist / hist.sum()
    p = p[p > 0]
    return -np.sum(p * np.log2(p))


def correlation_length(I, max_lag=50):
    I_centered = I - np.mean(I)
    var = np.var(I)
    if var < 1e-15:
        return 0
    acf = np.zeros(max_lag)
    for lag in range(max_lag):
        acf[lag] = np.mean(I_centered * np.roll(I_centered, lag)) / var
    for i in range(1, len(acf)):
        if acf[i] < 0:
            return np.sum(acf[:i])
    return np.sum(np.maximum(acf, 0))


# ============================================================
print("=" * 70)
print("TRANSIENCE TEST: Is S624's Class 4 sustained or transient?")
print("=" * 70)

N = 256
np.random.seed(42)
I0 = np.random.uniform(0.1, 0.9, N)

# Track entropy, ξ, and std over LONG runs
for A, k in [(1.0, 0.40), (1.0, 0.55), (1.5, 0.40), (0.8, 0.70)]:
    I = I0.copy()

    checkpoints = [100, 500, 1000, 2000, 5000, 10000, 20000]
    print(f"\n  A={A:.1f}, k={k:.2f}:")
    print(f"    {'Step':>8s} {'Entropy':>8s} {'ξ':>8s} {'std(I)':>8s} {'range':>14s}")
    print(f"    {'-'*50}")

    step = 0
    next_cp = 0
    while step <= 20000:
        if next_cp < len(checkpoints) and step == checkpoints[next_cp]:
            S = spatial_entropy(I)
            xi = correlation_length(I)
            std_I = np.std(I)
            print(f"    {step:8d} {S:8.3f} {xi:8.3f} {std_I:8.4f} [{np.min(I):.3f}, {np.max(I):.3f}]")
            next_cp += 1
        I = step_1d(I, k, A)
        step += 1


# ============================================================
print("\n" + "=" * 70)
print("SPATIAL-TEMPORAL ANTI-CORRELATION: The Entity Impossibility")
print("=" * 70)

print("""
In the 1-DOF CA, each cell has one real number (I). This number
encodes EVERYTHING about that cell's state. It must simultaneously:
  1. Carry spatial structure (vary between cells) — needs ξ > 1
  2. Carry temporal dynamics (change over time) — needs Δf > 0

With ONE field, these compete. Here's why:

- Spatial structure = gradients (cells differ from neighbors)
- Transfer rule = smoothing (I moves from high to low)
- Smoothing DESTROYS gradients → destroys spatial structure
- But smoothing IS the temporal dynamics

So the temporal dynamics (smoothing) destroys the spatial structure
(gradients). This is not a bug — it's the DEFINITION of diffusion.

The only escape: make the dynamics non-smooth (non-monotonic R, S624).
But then the dynamics destroys spatial structure even FASTER (chaos).

Low k: slow smoothing → gradients persist (high ξ) but no dynamics (low Δf)
High k: fast smoothing → gradients destroyed (low ξ) but active dynamics (high Δf)

Entities need both. This CA provides neither simultaneously.
""")

# Demonstrate: ξ vs Δf across k for A=1.0
k_vals = np.arange(0.20, 0.75, 0.03)
I0_test = np.random.uniform(0.1, 0.9, N)

xi_vals = []
df_vals = []

for k in k_vals:
    I = I0_test.copy()

    # Run 2000 steps
    for _ in range(1500):
        I = step_1d(I, k, A=1.0)

    # Measure ξ
    xi = correlation_length(I)
    xi_vals.append(xi)

    # Measure Δf from remaining 500 steps
    ts = np.zeros(500)
    for t in range(500):
        ts[t] = I[N//2]
        I = step_1d(I, k, A=1.0)

    fft = np.abs(np.fft.rfft(ts - np.mean(ts)))
    fft[0] = 0
    freqs = np.fft.rfftfreq(len(ts))
    power = fft**2
    tp = np.sum(power)
    if tp > 0:
        mf = np.sum(freqs * power) / tp
        vf = np.sum((freqs - mf)**2 * power) / tp
        df = np.sqrt(vf)
    else:
        df = 0
    df_vals.append(df)

xi_arr = np.array(xi_vals)
df_arr = np.array(df_vals)

print(f"  {'k':>5s} {'ξ':>8s} {'Δf':>8s}")
print(f"  {'-'*25}")
for i, k in enumerate(k_vals):
    marker = " ← both low" if xi_arr[i] < 3 and df_arr[i] < 0.01 else ""
    marker = " ← ξ high" if xi_arr[i] > 5 and df_arr[i] < 0.01 else marker
    marker = " ← Δf high" if df_arr[i] > 0.03 and xi_arr[i] < 3 else marker
    print(f"  {k:5.2f} {xi_arr[i]:8.3f} {df_arr[i]:8.4f}{marker}")

corr = np.corrcoef(xi_arr, df_arr)[0, 1] if np.std(xi_arr) > 0 and np.std(df_arr) > 0 else 0
print(f"\n  Correlation(ξ, Δf) = {corr:+.4f}")


# ============================================================
print("\n" + "=" * 70)
print("THE DEEPER PATTERN: Why Everything Maps to Known Physics")
print("=" * 70)

print("""
OBSERVATION: Every attempt to extract a prediction from Synchronism
maps to known physics. S617-625, nine sessions, nine mappings.

PATTERN: This isn't coincidence. It's structural. Here's why:

1. The framework CONCEPTS were derived from OBSERVATIONS of physics:
   - "Entity" → from observing particles
   - "Resonance/dissonance" → from observing interference
   - "MRH" → from observing scale-dependent phenomena
   - "Self-witnessing" → from observing stable structures
   - "C(ρ)" → from observing coherence in quantum systems

2. Concepts derived from physics observations → when formalized →
   must map back to the physics they were derived from.

3. This is CIRCULARITY, not confirmation:
   Physics → (abstract into vocabulary) → Synchronism concepts
   Synchronism concepts → (formalize) → Physics again

4. The only way to break this cycle:
   A concept that was NOT derived from physics observation.
   Something genuinely novel — not a reformulation.

5. Candidate novel concepts:
   - "Intent" — but defined as "unknowable greater force" → unfalsifiable
   - "Self-witnessing" — but maps to attractor dynamics → standard
   - "Indifference" — but maps to non-interaction → standard
   - "Paradigm shift over epicycles" — methodology, not prediction

6. RESULT: No genuinely novel concept exists in the framework.
   Every concept is either:
   (a) A reformulation of known physics (translatable)
   (b) A claim about unknowable things (unfalsifiable)

   Category (a) maps to physics. Category (b) makes no predictions.
   There is no category (c).

This is the ATTRACTOR MAP. Every probe lands in (a) or (b).
The space of novel predictions is not just empty — it's
structurally impossible given how the concepts were constructed.
""")

# ============================================================
# Check: did S624's Class 4 survive long runs?
# ============================================================
print("=" * 70)
print("S624 CLASS 4 REVISION")
print("=" * 70)

I = I0.copy()
# Run S624's exact Class 4 point for 20k steps
# Track entropy at regular intervals
entropies_long = []
for t in range(20000):
    if t % 100 == 0:
        entropies_long.append(spatial_entropy(I))
    I = step_1d(I, k=0.40, A=1.0)

S_arr = np.array(entropies_long)
print(f"\n  A=1.0, k=0.40 — entropy over 20,000 steps:")
print(f"    Steps 0-500:     S = {np.mean(S_arr[:5]):.3f} ± {np.std(S_arr[:5]):.4f}")
print(f"    Steps 500-2000:  S = {np.mean(S_arr[5:20]):.3f} ± {np.std(S_arr[5:20]):.4f}")
print(f"    Steps 2000-5000: S = {np.mean(S_arr[20:50]):.3f} ± {np.std(S_arr[20:50]):.4f}")
print(f"    Steps 5000-10k:  S = {np.mean(S_arr[50:100]):.3f} ± {np.std(S_arr[50:100]):.4f}")
print(f"    Steps 10k-20k:   S = {np.mean(S_arr[100:]):.3f} ± {np.std(S_arr[100:]):.4f}")
print(f"    Final step:      S = {S_arr[-1]:.3f}")

if S_arr[-1] < 1.0:
    print(f"\n  VERDICT: S624's Class 4 is TRANSIENT. System decays to Class 1.")
    print(f"  The edge-of-chaos regime exists but doesn't persist.")
    print(f"  Even the one positive finding from S624 is qualified.")
elif np.std(S_arr[100:]) < 0.01:
    print(f"\n  VERDICT: System reached steady state. Entropy sustained at {np.mean(S_arr[100:]):.3f}.")
    if np.mean(S_arr[100:]) > 2.0:
        print(f"  Class 4 SUSTAINED — complex but stable.")
    else:
        print(f"  Class 2 — periodic/simple.")
else:
    print(f"\n  VERDICT: System still fluctuating after 20k steps.")
    print(f"  Entropy {np.mean(S_arr[100:]):.3f} ± {np.std(S_arr[100:]):.4f}")
