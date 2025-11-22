"""measurement_collapse_sim.py

Simple conceptual simulation of "measurement" as a sharp ∂U/∂C event:

- Before measurement:
    pattern spans multiple C values.
- During measurement:
    large ∂U/∂C enforces collapse to a specific C_0.
- After measurement:
    pattern is localized in C (until further evolution).

This module is meant to be extended with:
- concrete U(Ψ, C) forms
- visualization hooks
- links to IIT-like measures (e.g., integrated information Φ).
"""

import numpy as np
from typing import Tuple


def collapse_to_C0(
    C: np.ndarray,
    C0: float,
    width: float = 0.1,
) -> np.ndarray:
    """Return a Gaussian-like mask peaked at C0.

    Parameters
    ----------
    C : np.ndarray
        Complexity axis values.
    C0 : float
        Collapse target.
    width : float
        Width of the "measurement" window.

    Returns
    -------
    np.ndarray
        Normalized mask over C.
    """
    mask = np.exp(-0.5 * ((C - C0) / width) ** 2)
    mask /= mask.sum() + 1e-12
    return mask


def apply_measurement(
    psi: np.ndarray,
    C: np.ndarray,
    C0: float,
    width: float = 0.1,
) -> np.ndarray:
    """Apply a toy "measurement collapse" to a field over C.

    Parameters
    ----------
    psi : np.ndarray
        Field values over C (1D or higher-dimensional with a C-axis).
    C : np.ndarray
        Complexity axis values.
    C0 : float
        Collapse target.
    width : float
        Collapse width.

    Returns
    -------
    np.ndarray
        Collapsed field.
    """
    mask = collapse_to_C0(C, C0, width=width)
    # If psi has extra dimensions, we broadcast the mask
    while mask.ndim < psi.ndim:
        mask = mask[np.newaxis, ...]
    return psi * mask
