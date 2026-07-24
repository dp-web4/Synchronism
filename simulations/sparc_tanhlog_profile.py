#!/usr/bin/env python3
"""Reconstruct the frozen SPARC likelihood profile for tanh-log gamma.

This is a reproducible extraction of the tanh-log portion of
``compander_family_aic_bic_real_sparc.py`` from synchronism-site commit
174eaf924b90cc6cc332e5c279560f325db02ae3.  It deliberately preserves that
analysis's row cut, fixed mass-to-light ratios, log-space residual, and
one-dimensional a0 profiling.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "simulations" / "sparc_real_data" / "MassModels_Lelli2016c.mrt"
KPC_M = 3.0856775814913673e19
KMS_M_S = 1.0e3
UPSILON_DISK = 0.5
UPSILON_BULGE = 0.7
ERROR_CUT = 0.10
SOURCE_COMMIT = "174eaf924b90cc6cc332e5c279560f325db02ae3"
SOURCE_SCRIPT_SHA256 = (
    "e90ea9a04290432039d0cf72f257dcbddb9b1f681125a7c24e29bc7476fa2bae"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def registered_gamma_grid() -> list[float]:
    """Return a grid satisfying the spacing and inclusion rules in registration."""
    values = set()
    values.update(np.arange(0.20, 0.35, 0.05).round(12))
    values.update(np.arange(0.35, 1.2500001, 0.025).round(12))
    values.update(np.arange(1.30, 5.0000001, 0.10).round(12))
    values.update({0.489, 0.5, 1.0, 2.0})
    return sorted(float(value) for value in values)


def load_accelerations(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load the exact 2026-07-22 row selection and baryonic prescription."""
    g_obs: list[float] = []
    g_bar: list[float] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) != 10:
                continue
            try:
                _, radius, v_obs, e_v_obs, v_gas, v_disk, v_bulge, _, _ = map(
                    float, parts[1:]
                )
            except ValueError:
                continue
            if radius <= 0 or v_obs <= 0 or e_v_obs / v_obs > ERROR_CUT:
                continue
            radius_m = radius * KPC_M
            v_bar_squared = (
                v_gas * abs(v_gas)
                + UPSILON_DISK * v_disk * abs(v_disk)
                + UPSILON_BULGE * v_bulge * abs(v_bulge)
            ) * KMS_M_S**2
            if v_bar_squared <= 0:
                continue
            g_bar.append(v_bar_squared / radius_m)
            g_obs.append((v_obs * KMS_M_S) ** 2 / radius_m)
    return np.asarray(g_obs), np.asarray(g_bar)


def tanhlog_prediction(g_bar: np.ndarray, a0: float, gamma: float) -> np.ndarray:
    """Invert g_bar = g_obs*tanh(gamma*ln(1+g_obs/a0)) vectorially."""
    target = g_bar / a0
    lower = target.copy()
    upper = np.maximum(2.0 * target + 1.0, 2.0)

    def residual(x: np.ndarray) -> np.ndarray:
        return x * np.tanh(gamma * np.log1p(x)) - target

    for _ in range(32):
        mask = residual(upper) < 0
        if not np.any(mask):
            break
        upper[mask] *= 2.0
    else:
        raise RuntimeError("failed to bracket tanh-log inversion")

    for _ in range(80):
        middle = 0.5 * (lower + upper)
        mask = residual(middle) < 0
        lower[mask] = middle[mask]
        upper[~mask] = middle[~mask]
    return a0 * 0.5 * (lower + upper)


def mcgaugh_prediction(g_bar: np.ndarray, a0: float) -> np.ndarray:
    y = g_bar / a0
    return g_bar / (-np.expm1(-np.sqrt(y)))


def ssr_log(g_obs: np.ndarray, prediction: np.ndarray) -> float:
    residual = np.log10(g_obs) - np.log10(prediction)
    return float(np.dot(residual, residual))


def fit_a0(
    g_obs: np.ndarray,
    g_bar: np.ndarray,
    prediction,
) -> tuple[float, float]:
    objective = lambda log_a0: ssr_log(
        g_obs,
        prediction(g_bar, 10.0**log_a0),
    )
    result = minimize_scalar(
        objective,
        bounds=(-11.0, -9.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    if not result.success:
        raise RuntimeError(f"a0 profile failed: {result.message}")
    return 10.0**float(result.x), float(result.fun)


def build_profile(data_path: Path) -> dict[str, object]:
    g_obs, g_bar = load_accelerations(data_path)
    sample_size = len(g_obs)
    if sample_size != 2807:
        raise RuntimeError(f"expected 2807 selected SPARC rows, found {sample_size}")

    reference_a0, reference_ssr = fit_a0(g_obs, g_bar, mcgaugh_prediction)
    reference_bic = sample_size * math.log(reference_ssr / sample_size) + math.log(
        sample_size
    )

    rows: list[dict[str, float]] = []
    for gamma in registered_gamma_grid():
        a0, ssr = fit_a0(
            g_obs,
            g_bar,
            lambda values, scale, shape=gamma: tanhlog_prediction(
                values, scale, shape
            ),
        )
        bic_fixed_gamma = sample_size * math.log(ssr / sample_size) + math.log(
            sample_size
        )
        rows.append(
            {
                "gamma": gamma,
                "a0_m_s2": a0,
                "ssr_log10": ssr,
                "rms_dex": math.sqrt(ssr / sample_size),
                "delta_bic_fixed_gamma_vs_mcgaugh": bic_fixed_gamma
                - reference_bic,
                "delta_bic_family_k2_vs_mcgaugh": bic_fixed_gamma
                + math.log(sample_size)
                - reference_bic,
            }
        )

    minimum_ssr = min(row["ssr_log10"] for row in rows)
    for row in rows:
        row["delta_bic_profile"] = sample_size * math.log(
            row["ssr_log10"] / minimum_ssr
        )

    best = min(rows, key=lambda row: row["ssr_log10"])
    return {
        "artifact": "frozen SPARC tanh-log gamma profile",
        "source_analysis": {
            "repository": "dp-web4/synchronism-site",
            "commit": SOURCE_COMMIT,
            "script": "explorer/scripts/compander_family_aic_bic_real_sparc.py",
            "script_sha256": SOURCE_SCRIPT_SHA256,
        },
        "data": {
            "path": str(data_path.relative_to(ROOT)),
            "sha256": sha256_file(data_path),
            "selected_rows": sample_size,
            "selection": "eVobs/Vobs <= 0.10; positive R, Vobs, and Vbar^2",
            "upsilon_disk": UPSILON_DISK,
            "upsilon_bulge": UPSILON_BULGE,
        },
        "likelihood": {
            "residual": "log10(g_obs) - log10(g_obs_model)",
            "objective": "unweighted sum of squared log10 residuals",
            "a0_bounds_log10": [-11.0, -9.0],
            "bic": "N*ln(SSR/N) + k*ln(N)",
        },
        "reference_mcgaugh": {
            "a0_m_s2": reference_a0,
            "ssr_log10": reference_ssr,
            "rms_dex": math.sqrt(reference_ssr / sample_size),
            "parameter_count": 1,
        },
        "grid_rule": (
            "registration section 3; includes 0.489, 0.5, 1, 2; "
            "spacing <=0.025 on [0.35,1.25] and <=0.10 elsewhere"
        ),
        "bic_conventions": {
            "delta_bic_profile": (
                "profile likelihood relative to best grid row; parameter count cancels"
            ),
            "delta_bic_fixed_gamma_vs_mcgaugh": (
                "each externally fixed gamma has k=1, compared with k=1 McGaugh"
            ),
            "delta_bic_family_k2_vs_mcgaugh": (
                "gamma treated as fitted family parameter (k=2), compared with k=1 McGaugh"
            ),
        },
        "best_grid_row": best,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = build_profile(args.data.resolve())
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
