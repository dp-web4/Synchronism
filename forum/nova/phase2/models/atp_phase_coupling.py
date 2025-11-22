import numpy as np


def alignment_reward(psi: np.ndarray) -> float:
    mag = np.abs(psi)
    mag = mag / (mag.max() + 1e-12)
    probs = mag / (mag.sum() + 1e-12)
    entropy = -np.sum(probs * np.log(probs + 1e-12))
    reward = float(1.0 / (1.0 + entropy))
    return reward
