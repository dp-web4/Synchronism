from typing import Callable, Sequence
import numpy as np


def superpose_modes(
    modes: Sequence[Callable[[float, np.ndarray], np.ndarray]],
    t: float,
    C: np.ndarray,
) -> np.ndarray:
    total = None
    for m in modes:
        psi = m(t, C)
        if total is None:
            total = psi
        else:
            total = total + psi
    if total is None:
        raise ValueError("No modes provided")
    return total


def coherence_measure(psi: np.ndarray) -> float:
    mag = np.abs(psi)
    mag_norm = mag / (mag.max() + 1e-12)
    var = np.var(mag_norm)
    return float(1.0 / (1.0 + var))
