#!/usr/bin/env python3
"""QUMOND Solar-System quadrupole for registered interpolation functions.

This instrument implements Eq. (12) of Desmond, Hees & Famaey (2024),
equivalent to Eqs. (9) of Park et al. (2026). Its benchmark mode is allowed
before target-family execution by the accompanying preregistration.
"""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq


G_SI = 6.67430e-11
M_SUN_KG = 1.98847e30


ArrayFunction = Callable[[np.ndarray], np.ndarray]


def nu_rar(y: np.ndarray) -> np.ndarray:
    """McGaugh-Lelli-Schombert/RAR nu function."""
    y = np.asarray(y, dtype=float)
    return 1.0 / (-np.expm1(-np.sqrt(y)))


def nu_n_family(y: np.ndarray, n: float) -> np.ndarray:
    """MOND n-family; n=1 Simple and n=2 Standard."""
    y = np.asarray(y, dtype=float)
    return (0.5 * (1.0 + np.sqrt(1.0 + 4.0 * y ** (-n)))) ** (1.0 / n)


def mu_tanh_log(x: np.ndarray, gamma: float) -> np.ndarray:
    """Registered tanh-log family in the mu convention."""
    x = np.asarray(x, dtype=float)
    return np.tanh(gamma * np.log1p(x))


class TanhLogNu:
    """Monotone numerical conversion from mu_gamma(x) to nu_gamma(y)."""

    def __init__(
        self,
        gamma: float,
        *,
        log_x_min: float = -32.0,
        log_x_max: float = 32.0,
        grid_size: int = 20001,
    ) -> None:
        if gamma <= 0:
            raise ValueError("gamma must be positive")
        self.gamma = float(gamma)
        self.log_x_min = float(log_x_min)
        self.log_x_max = float(log_x_max)

        log_x = np.linspace(log_x_min, log_x_max, grid_size)
        x = np.exp(log_x)
        mu = mu_tanh_log(x, self.gamma)
        y = x * mu
        nu = 1.0 / mu
        self._log_y_min = float(np.log(y[0]))
        self._log_y_max = float(np.log(y[-1]))
        self._log_nu_of_log_y = PchipInterpolator(
            np.log(y),
            np.log(nu),
            extrapolate=False,
        )

    def __call__(self, y: np.ndarray) -> np.ndarray:
        y = np.asarray(y, dtype=float)
        if np.any(y <= 0):
            raise ValueError("nu is defined here only for y > 0")
        log_y = np.log(y)
        result = np.empty_like(y)

        low = log_y < self._log_y_min
        high = log_y > self._log_y_max
        middle = ~(low | high)

        # mu ~ gamma*x gives y ~ gamma*x^2 and nu ~ 1/sqrt(gamma*y).
        result[low] = 1.0 / np.sqrt(self.gamma * y[low])
        # mu ~ 1 - 2*x^(-2 gamma); y and x coincide asymptotically.
        result[high] = 1.0 + 2.0 * y[high] ** (-2.0 * self.gamma)
        result[middle] = np.exp(self._log_nu_of_log_y(log_y[middle]))
        return result

    def mapping_relative_error(self, y: np.ndarray) -> np.ndarray:
        """Return error in y*nu*mu(y*nu) = y."""
        y = np.asarray(y, dtype=float)
        nu = self(y)
        x = y * nu
        reconstructed = y * nu * mu_tanh_log(x, self.gamma)
        return np.abs(reconstructed - y) / y


def solve_e_n(nu: ArrayFunction, e_tilde: float) -> float:
    """Solve e_N * nu(e_N) = e_tilde."""
    if e_tilde <= 0:
        raise ValueError("e_tilde must be positive")

    def residual(log_e_n: float) -> float:
        e_n = math.exp(log_e_n)
        return e_n * float(nu(np.array([e_n]))[0]) - e_tilde

    return math.exp(brentq(residual, -40.0, math.log(e_tilde) + 2.0))


@dataclass(frozen=True)
class QResult:
    e_tilde: float
    e_n: float
    q: float
    q_abs: float
    angular_order: int
    radial_order: int
    log_v_bound: float
    estimated_error: float


def qumond_q(
    nu: ArrayFunction,
    e_tilde: float,
    *,
    angular_order: int = 128,
    radial_order: int = 512,
    log_v_bound: float = 14.0,
) -> QResult:
    """Evaluate the dimensionless QUMOND quadrupole q.

    v is integrated in log space. An angle-independent nu reference is
    subtracted from the integrand; its contribution is identically zero
    because both angular polynomials integrate to zero. This is the
    convergence stabilization described by Hees et al. (2016).
    """
    if angular_order < 16:
        raise ValueError("angular_order must be at least 16")
    if radial_order < 64 or radial_order % 2:
        raise ValueError("radial_order must be an even integer of at least 64")
    if log_v_bound <= 0:
        raise ValueError("log_v_bound must be positive")

    e_n = solve_e_n(nu, e_tilde)
    xi, weights = leggauss(angular_order)

    def log_v_integrand(log_v: float) -> float:
        v = math.exp(log_v)
        v2 = v * v
        argument = np.sqrt(e_n * e_n + v2 * v2 + 2.0 * e_n * v2 * xi)
        reference_argument = math.sqrt(e_n * e_n + v2 * v2)
        delta_nu = nu(argument) - float(
            nu(np.array([reference_argument], dtype=float))[0]
        )
        angular_polynomial = (
            e_n * (3.0 * xi - 5.0 * xi**3)
            + v2 * (1.0 - 3.0 * xi**2)
        )
        # dv = v d(log v)
        return v * float(np.dot(weights, delta_nu * angular_polynomial))

    # The total Newtonian field can become small near v^2=e_N, xi=-1.
    # Split at that analytically known location rather than allowing the
    # adaptive integrator to spend its entire subdivision budget finding it.
    cancellation_log_v = 0.5 * math.log(e_n)
    split_points = sorted(
        {
            -log_v_bound,
            min(max(cancellation_log_v, -log_v_bound), log_v_bound),
            log_v_bound,
        }
    )
    def fixed_integral(order: int) -> float:
        nodes, radial_weights = leggauss(order)
        total = 0.0
        for lower, upper in zip(split_points, split_points[1:]):
            half_width = 0.5 * (upper - lower)
            midpoint = 0.5 * (upper + lower)
            log_v_nodes = midpoint + half_width * nodes
            values = np.fromiter(
                (log_v_integrand(value) for value in log_v_nodes),
                dtype=float,
                count=order,
            )
            total += half_width * float(np.dot(radial_weights, values))
        return total

    integral = fixed_integral(radial_order)
    lower_order_integral = fixed_integral(radial_order // 2)
    error = abs(integral - lower_order_integral)
    q = 1.5 * integral
    return QResult(
        e_tilde=e_tilde,
        e_n=e_n,
        q=q,
        q_abs=abs(q),
        angular_order=angular_order,
        radial_order=radial_order,
        log_v_bound=log_v_bound,
        estimated_error=1.5 * error,
    )


def q2_si(q: float, a0: float) -> float:
    """Convert dimensionless q to signed Q2 in s^-2."""
    if a0 <= 0:
        raise ValueError("a0 must be positive")
    return -3.0 * a0**1.5 * q / (2.0 * math.sqrt(G_SI * M_SUN_KG))


def benchmark() -> dict[str, object]:
    """Run only literature benchmarks, never the target tanh-log family."""
    expected = {1.0: 0.094, 1.5: 0.159, 2.0: 0.221}
    rows = []
    for e_tilde, expected_abs_q in expected.items():
        primary = qumond_q(nu_rar, e_tilde)
        angular = qumond_q(nu_rar, e_tilde, angular_order=256)
        radial = qumond_q(nu_rar, e_tilde, radial_order=1024)
        bounds = qumond_q(nu_rar, e_tilde, log_v_bound=16.0)
        relative_error = abs(primary.q_abs - expected_abs_q) / expected_abs_q
        angular_change = abs(angular.q_abs - primary.q_abs) / primary.q_abs
        radial_change = abs(radial.q_abs - primary.q_abs) / primary.q_abs
        bounds_change = abs(bounds.q_abs - primary.q_abs) / primary.q_abs
        rows.append(
            {
                **asdict(primary),
                "expected_abs_q": expected_abs_q,
                "relative_error": relative_error,
                "angular_order_change": angular_change,
                "radial_order_change": radial_change,
                "log_bound_change": bounds_change,
                "pass": (
                    relative_error < 0.02
                    and angular_change < 0.005
                    and radial_change < 0.005
                    and bounds_change < 0.005
                ),
            }
        )

    asymptotic_y = np.array([1e-12, 1e12])
    simple = nu_n_family(asymptotic_y, 1.0)
    standard = nu_n_family(asymptotic_y, 2.0)
    asymptotes_pass = bool(
        np.isclose(simple[0] * math.sqrt(asymptotic_y[0]), 1.0, rtol=2e-6)
        and np.isclose(standard[0] * math.sqrt(asymptotic_y[0]), 1.0, rtol=2e-6)
        and abs(simple[1] - 1.0) < 2e-12
        and abs(standard[1] - 1.0) < 2e-12
    )
    return {
        "benchmark": "RAR q values from Desmond, Hees & Famaey (2024)",
        "rows": rows,
        "n_family_asymptotes_pass": asymptotes_pass,
        "all_pass": all(row["pass"] for row in rows) and asymptotes_pass,
    }


def target_scan(
    gammas: Sequence[float],
    *,
    a0: float,
    g_ext: float,
    angular_order: int,
    radial_order: int,
    log_v_bound: float,
) -> dict[str, object]:
    """Preparatory Solar-System-only target scan.

    This is not the registered joint test because a0 must ultimately be
    profiled from the frozen SPARC surface for every gamma.
    """
    rows = []
    mapping_y = np.logspace(-12, 12, 401)
    for gamma in gammas:
        nu = TanhLogNu(gamma)
        mapping_error = float(np.max(nu.mapping_relative_error(mapping_y)))
        result = qumond_q(
            nu,
            g_ext / a0,
            angular_order=angular_order,
            radial_order=radial_order,
            log_v_bound=log_v_bound,
        )
        q2 = q2_si(result.q, a0)
        integration_relative_error = (
            result.estimated_error / result.q_abs if result.q_abs else math.inf
        )
        rows.append(
            {
                "gamma": gamma,
                "a0_m_s2": a0,
                "g_ext_m_s2": g_ext,
                **asdict(result),
                "q2_s2": q2,
                "q2_abs_s2": abs(q2),
                "cassini_z": (q2 - 1.6e-27) / 1.8e-27,
                "max_mapping_relative_error": mapping_error,
                "mapping_pass": mapping_error < 1e-7,
                "integration_relative_error": integration_relative_error,
                "integration_pass": integration_relative_error < 0.005,
            }
        )
    return {
        "status": (
            "PREPARATORY SOLAR-SYSTEM-ONLY SCAN; not the registered joint "
            "test because a0 is not profiled from SPARC"
        ),
        "cassini_2026_mean_s2": 1.6e-27,
        "cassini_2026_sigma_s2": 1.8e-27,
        "rows": rows,
    }


def parse_gamma_list(raw: str) -> list[float]:
    values = [float(value) for value in raw.split(",") if value.strip()]
    if not values or any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("provide positive comma-separated gammas")
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--benchmark", action="store_true")
    mode.add_argument(
        "--target-gammas",
        type=parse_gamma_list,
        help="Preparatory scan, e.g. 0.489,0.5,1,2",
    )
    parser.add_argument("--a0", type=float, default=1.03e-10)
    parser.add_argument("--g-ext", type=float, default=2.32e-10)
    parser.add_argument("--angular-order", type=int, default=128)
    parser.add_argument("--radial-order", type=int, default=512)
    parser.add_argument("--log-v-bound", type=float, default=14.0)
    args = parser.parse_args()

    if args.benchmark:
        result = benchmark()
    else:
        result = target_scan(
            args.target_gammas,
            a0=args.a0,
            g_ext=args.g_ext,
            angular_order=args.angular_order,
            radial_order=args.radial_order,
            log_v_bound=args.log_v_bound,
        )
    print(json.dumps(result, indent=2))
    return 0 if result.get("all_pass", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
