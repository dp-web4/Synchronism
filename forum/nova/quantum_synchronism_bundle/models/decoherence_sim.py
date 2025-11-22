"""decoherence_sim.py

Minimal toy simulation for ΔC-driven decoherence.

This is not meant to be physically exact. It is a conceptual scaffold to explore:

- how increasing D_c (complexity diffusion) destroys phase-coherent modes
- how noise in k, ω, or U(Ψ) affects pattern stability.
"""

import numpy as np
from typing import Tuple


def simulate_decoherence(
    C: np.ndarray,
    t_steps: int,
    omega: float,
    k: float,
    D_c: float,
    noise_std: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate a simple decoherence process over complexity C and time.

    Parameters
    ----------
    C : np.ndarray
        Complexity axis values.
    t_steps : int
        Number of time steps.
    omega : float
        Base temporal frequency.
    k : float
        Base complexity wave number.
    D_c : float
        Complexity diffusion coefficient (larger → faster decoherence).
    noise_std : float
        Standard deviation of Gaussian noise applied to k at each step.

    Returns
    -------
    phases : np.ndarray
        Phase values over time and complexity: shape (t_steps, len(C)).
    amplitudes : np.ndarray
        Toy amplitudes (currently unity, can be extended).
    """
    phases = np.zeros((t_steps, len(C)))
    amplitudes = np.ones_like(phases)

    current_k = k
    for t in range(t_steps):
        # Simple random walk in k as a stand-in for environmental noise
        if noise_std > 0.0:
            current_k += np.random.normal(0.0, noise_std)

        phase = omega * t - current_k * C

        # D_c can be interpreted as how fast phase spreads / randomizes
        decoherence_factor = np.exp(-D_c * t)
        phases[t] = phase
        amplitudes[t] = decoherence_factor

    return phases, amplitudes
