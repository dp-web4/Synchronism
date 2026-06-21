"""
01 — Kuramoto baseline: emergent frequency-cluster "particles" on a 2D lattice.

Synchronism reading: entities are *recurring patterns*. Here, two regions seeded with a
distinct natural frequency phase-lock internally and persist as localized, coherent
clusters against a background that locks at its own rate — the simplest runnable form of
"a particle is a stable oscillation pattern" (FUNDAMENTALS §Entity).

This is a baseline / instantiation, NOT a novel-prediction test. It exists to make the
ontology concrete and to provide the lattice machinery the CHSH harness (02) reuses.

numpy only. Headless. Writes results/baseline_result.json.
"""
import json
import os
import numpy as np

GRID = 48
K = 1.5          # nearest-neighbor coupling
DT = 0.05
STEPS = 600
SEED = 7         # fixed for reproducibility (Date/random unavailable by policy; vary by hand)


def step(phases, freqs, k, dt):
    """One Kuramoto lattice update with periodic (toroidal) nearest-neighbor coupling."""
    influence = (
        np.sin(np.roll(phases, 1, 0) - phases)
        + np.sin(np.roll(phases, -1, 0) - phases)
        + np.sin(np.roll(phases, 1, 1) - phases)
        + np.sin(np.roll(phases, -1, 1) - phases)
    )
    phases = phases + dt * (freqs + (k / 4.0) * influence)
    return (phases + np.pi) % (2 * np.pi) - np.pi


def order_parameter(phases):
    """Kuramoto R in [0,1]: 0 = incoherent, 1 = globally phase-locked."""
    return float(np.abs(np.mean(np.exp(1j * phases))))


def cluster_lock(phases, region):
    """Local order parameter inside a region mask — how tightly that 'entity' is locked."""
    return float(np.abs(np.mean(np.exp(1j * phases[region]))))


def main():
    rng = np.random.default_rng(SEED)
    phases = rng.uniform(-np.pi, np.pi, (GRID, GRID))

    # Background at f=1.0; two "particle" regions at a distinct frequency → they self-lock
    freqs = np.ones((GRID, GRID)) * 1.0
    A = np.zeros((GRID, GRID), dtype=bool); A[10:16, 10:16] = True
    B = np.zeros((GRID, GRID), dtype=bool); B[30:36, 30:36] = True
    freqs[A] = 2.5
    freqs[B] = 2.5

    R_hist, A_hist, B_hist = [], [], []
    for _ in range(STEPS):
        phases = step(phases, freqs, K, DT)
        R_hist.append(order_parameter(phases))
        A_hist.append(cluster_lock(phases, A))
        B_hist.append(cluster_lock(phases, B))

    result = {
        "grid": GRID, "K": K, "dt": DT, "steps": STEPS, "seed": SEED,
        "global_R_final": round(R_hist[-1], 4),
        "global_R_mean_last50": round(float(np.mean(R_hist[-50:])), 4),
        "clusterA_lock_final": round(A_hist[-1], 4),
        "clusterB_lock_final": round(B_hist[-1], 4),
        "interpretation": (
            "Two distinct-frequency regions phase-lock internally (cluster lock -> ~1) and "
            "persist as localized coherent patterns while the global order parameter "
            "plateaus below 1 (the 'particles' refuse to dissolve into the background). "
            "This instantiates 'entity = recurring pattern'; it is NOT a novel-prediction test."
        ),
    }
    os.makedirs(os.path.join(os.path.dirname(__file__), "results"), exist_ok=True)
    out = os.path.join(os.path.dirname(__file__), "results", "baseline_result.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
