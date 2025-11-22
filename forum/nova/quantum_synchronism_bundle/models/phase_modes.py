"""phase_modes.py

Utilities for representing and manipulating phase-mode solutions of the form:

    Ψ_q(x, t, C) = A(x) * exp(i (ω t - k C))

This module is intentionally simple and symbolic. It can be extended with:

- sampling functions
- visualization helpers
- stability checks for given (ω, k).
"""

import numpy as np
from typing import Callable


def psi_mode(
    A_x: Callable[[np.ndarray], np.ndarray],
    omega: float,
    k: float,
    t: float,
    C: np.ndarray,
) -> np.ndarray:
    """Compute Ψ_q(x, t, C) for a given spatial envelope A(x).

    Parameters
    ----------
    A_x : callable
        Function that maps x-array → amplitude array A(x).
        (x itself may be implicitly encoded in the shape of A_x's output.)
    omega : float
        Temporal frequency ω.
    k : float
        Complexity wave number k.
    t : float
        Time.
    C : np.ndarray
        Complexity coordinates.

    Returns
    -------
    np.ndarray
        Complex array representing Ψ_q at the given t and C.
    """
    # For now we assume A_x already encodes x; we just broadcast over C.
    A = A_x(C)  # mild abuse; replace with more structured treatment later
    phase = omega * t - k * C
    return A * np.exp(1j * phase)
