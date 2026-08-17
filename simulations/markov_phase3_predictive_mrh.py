#!/usr/bin/env python3
"""Phase 3 toy model: predictive MRH versus raw causal cone.

A one-dimensional noisy majority cellular process is used only as a transparent
local stochastic dynamical system. For forecast horizon h, the exact causal cone
of the center has radius h. We estimate how much additional predictive information
about the center's future remains in the exterior portion of that cone once a local
radius r is already known:

    epsilon(r,h) = I(Y_{t+h}; E_{r<|x|<=h,t} | B_{|x|<=r,t})

The effective MRH at tolerance epsilon* is the smallest r with epsilon(r,h)<=epsilon*.

This is a formalism probe, not a physics simulation.
"""

import math
import numpy as np


def binary_entropy(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


def conditional_entropy_binary(y, z):
    """Empirical H(Y|Z), where Y is binary and columns of Z are binary."""
    if z.shape[1] == 0:
        return binary_entropy(float(y.mean()))

    powers = 1 << np.arange(z.shape[1], dtype=np.uint64)
    keys = (z.astype(np.uint64) * powers).sum(axis=1).astype(np.int64)
    total = np.bincount(keys)
    ones = np.bincount(keys, weights=y, minlength=len(total))

    n = len(y)
    h = 0.0
    for count, one_count in zip(total, ones):
        if count:
            h += (count / n) * binary_entropy(one_count / count)
    return h


def estimate(horizon, radius, samples=200_000, noise=0.05, seed=1234):
    """Estimate conditional mutual information in bits."""
    rng = np.random.default_rng(seed + 100 * horizon + radius)
    width = 2 * horizon + 1
    initial = rng.integers(0, 2, size=(samples, width), dtype=np.uint8)
    current = initial.copy()

    # Synchronous noisy majority-of-three update. Each step shrinks the represented
    # causal cone by one cell on each side, leaving the target center after h steps.
    for _ in range(horizon):
        summed = current[:, :-2] + current[:, 1:-1] + current[:, 2:]
        nxt = (summed >= 2).astype(np.uint8)
        flips = (rng.random(nxt.shape) < noise).astype(np.uint8)
        current = np.bitwise_xor(nxt, flips)

    y = current[:, 0]
    center = horizon
    local = initial[:, center - radius:center + radius + 1]
    full_cone = initial

    h_local = conditional_entropy_binary(y, local)
    h_full = conditional_entropy_binary(y, full_cone)
    cmi = h_local - h_full
    return cmi, h_local, h_full


def main():
    tolerance = 0.02
    print(f"noise=0.05, relevance tolerance={tolerance:.3f} bits")
    print("h  r  exterior_CMI_bits")

    for h in range(1, 5):
        effective = None
        for r in range(h + 1):
            cmi, _, _ = estimate(h, r)
            print(f"{h:1d}  {r:1d}  {cmi:.6f}")
            if effective is None and cmi <= tolerance:
                effective = r
        print(f"  causal radius={h}; predictive MRH radius={effective}\n")


if __name__ == "__main__":
    main()
