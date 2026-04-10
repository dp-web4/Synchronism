"""
Session 623: Verify whether 2D "gliders" are real or diffusion artifacts.

The 2D test showed COM drift of ~22 cells from (16,16) toward (32,32).
Hypothesis: this is just diffusion toward the grid center, not a moving structure.

Test: place pulse at grid center. If COM doesn't move, the previous drift was
just convergence to the uniform distribution's center of mass.
Also test: track the PEAK position, not COM. A real glider moves its peak.
"""

import numpy as np

I_MAX = 1.0
N = 64


def R(I, n=2):
    return 1.0 - (np.clip(I, 0, I_MAX) / I_MAX) ** n


def step_2d(I, k, n=2):
    up = np.roll(I, 1, axis=0)
    down = np.roll(I, -1, axis=0)
    left = np.roll(I, 1, axis=1)
    right = np.roll(I, -1, axis=1)
    dI = k * ((up - I) * R(up, n) + (down - I) * R(down, n) +
              (left - I) * R(left, n) + (right - I) * R(right, n))
    return np.clip(I + dI, 0.0, I_MAX)


steps = 1000

for label, cx, cy in [("corner (16,16)", 16, 16), ("center (32,32)", 32, 32), ("off-center (20,40)", 20, 40)]:
    print(f"\n=== Pulse at {label} ===")
    for k in [0.3, 0.8]:
        I = np.ones((N, N)) * 0.3
        y, x = np.mgrid[:N, :N]
        pulse = 0.6 * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2 * 3 ** 2))
        I = np.clip(I + pulse, 0, I_MAX)

        peaks_x = []
        peaks_y = []
        peak_vals = []

        for t in range(steps):
            I = step_2d(I, k)
            if t % 100 == 0:
                diff = I - I.mean()  # Remove background
                idx = np.unravel_index(np.argmax(np.abs(diff)), diff.shape)
                peaks_y.append(idx[0])
                peaks_x.append(idx[1])
                peak_vals.append(np.max(np.abs(diff)))

        print(f"  k={k}: peak positions x: {peaks_x[:5]}...{peaks_x[-3:]}")
        print(f"  k={k}: peak positions y: {peaks_y[:5]}...{peaks_y[-3:]}")
        print(f"  k={k}: peak values: {peak_vals[0]:.4f} → {peak_vals[-1]:.4f}")
        peak_moved = max(abs(peaks_x[-1] - peaks_x[0]), abs(peaks_y[-1] - peaks_y[0]))
        print(f"  k={k}: peak displacement: {peak_moved} cells → {'MOVED' if peak_moved > 2 else 'STATIONARY'}")
