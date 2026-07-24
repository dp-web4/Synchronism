#!/usr/bin/env python3
"""Execute the preregistered SPARC x Cassini set-intersection test."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
from pathlib import Path

import numpy as np
import scipy

from sparc_cassini_q2 import TanhLogNu, benchmark, mu_tanh_log, q2_si, qumond_q


ROOT = Path(__file__).resolve().parents[1]
CASSINI_2026 = {"mean_s2": 1.6e-27, "sigma_s2": 1.8e-27}
CASSINI_2014 = {"mean_s2": 3.0e-27, "sigma_s2": 3.0e-27}
G_EXT_VALUES = (2.00e-10, 2.32e-10, 2.48e-10)
PRIMARY_G_EXT = 2.32e-10
PROFILE_CONVENTIONS = (
    "delta_bic_profile",
    "delta_bic_fixed_gamma_vs_mcgaugh",
    "delta_bic_family_k2_vs_mcgaugh",
)
RADIAL_ORDERS = (512, 1024, 2048, 4096, 8192)
CONVERGENCE_TOLERANCE = 0.005


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_head() -> str:
    return subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def interval(measurement: dict[str, float]) -> tuple[float, float]:
    mean = measurement["mean_s2"]
    width = 1.96 * measurement["sigma_s2"]
    return mean - width, mean + width


def inside(value: float, bounds: tuple[float, float]) -> bool:
    return bounds[0] <= value <= bounds[1]


def grid_intervals(values: list[float], complete_grid: list[float]) -> list[list[float]]:
    """Compress surviving grid points without implying interpolated boundaries."""
    if not values:
        return []
    index = {value: position for position, value in enumerate(complete_grid)}
    result: list[list[float]] = []
    start = previous = values[0]
    for value in values[1:]:
        if index[value] != index[previous] + 1:
            result.append([start, previous])
            start = value
        previous = value
    result.append([start, previous])
    return result


def converged_quadrupole(
    nu: TanhLogNu,
    e_tilde: float,
) -> tuple[object, list[dict[str, float]], float]:
    """Require two consecutive radial convergences plus an angular check."""
    previous = None
    consecutive_passes = 0
    history: list[dict[str, float]] = []
    accepted = None
    for radial_order in RADIAL_ORDERS:
        result = qumond_q(
            nu,
            e_tilde,
            angular_order=128,
            radial_order=radial_order,
        )
        relative_change = (
            None
            if previous is None
            else abs(result.q - previous.q) / result.q_abs
        )
        history.append(
            {
                "radial_order": radial_order,
                "q": result.q,
                "change_from_previous": relative_change,
            }
        )
        if relative_change is not None and relative_change < CONVERGENCE_TOLERANCE:
            consecutive_passes += 1
        else:
            consecutive_passes = 0
        previous = result
        if consecutive_passes >= 2:
            accepted = result
            break
    if accepted is None:
        raise RuntimeError(
            f"radial quadrature did not converge through order {RADIAL_ORDERS[-1]}"
        )
    angular = qumond_q(
        nu,
        e_tilde,
        angular_order=256,
        radial_order=accepted.radial_order,
    )
    angular_change = abs(angular.q - accepted.q) / accepted.q_abs
    return accepted, history, angular_change


def execute(profile_path: Path) -> dict[str, object]:
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    profile_rows = profile["rows"]
    full_grid = [row["gamma"] for row in profile_rows]
    benchmark_result = benchmark()
    current_interval = interval(CASSINI_2026)
    legacy_interval = interval(CASSINI_2014)

    rows: list[dict[str, object]] = []
    for profile_row in profile_rows:
        gamma = profile_row["gamma"]
        a0 = profile_row["a0_m_s2"]
        sparc_relevant = any(
            profile_row[convention] <= 14.0
            for convention in PROFILE_CONVENTIONS
        )
        if not sparc_relevant:
            rows.append(
                {
                    **profile_row,
                    "cassini_evaluation_status": (
                        "not_required: excluded by SPARC at delta_BIC > 14 "
                        "under every registered convention"
                    ),
                    "max_mapping_relative_error": None,
                    "newtonian_residual_1_minus_mu_at_x30": float(
                        1.0 - mu_tanh_log(np.array([30.0]), gamma)[0]
                    ),
                    "solar_system": [],
                }
            )
            continue

        nu = TanhLogNu(gamma)
        mapping_error = float(
            np.max(nu.mapping_relative_error(np.logspace(-12, 12, 401)))
        )
        g_ext_rows = []
        for g_ext in G_EXT_VALUES:
            primary, radial_history, angular_change = converged_quadrupole(
                nu, g_ext / a0
            )
            q2 = q2_si(primary.q, a0)
            sensitivity: dict[str, object] = {
                "g_ext_m_s2": g_ext,
                "q": primary.q,
                "q2_s2": q2,
                "cassini_2026_z": (
                    q2 - CASSINI_2026["mean_s2"]
                )
                / CASSINI_2026["sigma_s2"],
                "cassini_2026_pass_95": inside(q2, current_interval),
                "cassini_2014_pass_95": inside(q2, legacy_interval),
                "integration_estimated_relative_error": (
                    primary.estimated_error / primary.q_abs
                ),
                "radial_convergence_history": radial_history,
                "angular_128_to_256_relative_change": angular_change,
            }
            g_ext_rows.append(sensitivity)

        rows.append(
            {
                **profile_row,
                "cassini_evaluation_status": "evaluated: inside union of SPARC sensitivity sets",
                "max_mapping_relative_error": mapping_error,
                "newtonian_residual_1_minus_mu_at_x30": float(
                    1.0 - mu_tanh_log(np.array([30.0]), gamma)[0]
                ),
                "solar_system": g_ext_rows,
            }
        )

    summary: dict[str, object] = {}
    for convention in PROFILE_CONVENTIONS:
        convention_summary: dict[str, object] = {}
        for threshold in (6.0, 10.0, 14.0):
            threshold_summary: dict[str, object] = {}
            sparc_values = [
                row["gamma"] for row in rows if row[convention] <= threshold
            ]
            for g_ext_index, g_ext in enumerate(G_EXT_VALUES):
                cassini_values = [
                    row["gamma"]
                    for row in rows
                    if row["solar_system"]
                    and row["solar_system"][g_ext_index]["cassini_2026_pass_95"]
                ]
                joint_values = sorted(set(sparc_values) & set(cassini_values))
                threshold_summary[f"g_ext_{g_ext:.2e}"] = {
                    "sparc_grid_intervals": grid_intervals(
                        sparc_values, full_grid
                    ),
                    "cassini_grid_intervals": grid_intervals(
                        cassini_values, full_grid
                    ),
                    "joint_grid_intervals": grid_intervals(
                        joint_values, full_grid
                    ),
                    "joint_grid_count": len(joint_values),
                }
            convention_summary[f"delta_bic_le_{threshold:g}"] = threshold_summary
        summary[convention] = convention_summary

    numerical_pass = all(
        row["max_mapping_relative_error"] is None
        or (
        row["max_mapping_relative_error"] < 1e-7
        and all(
            sensitivity["integration_estimated_relative_error"] < 0.005
            and sensitivity["angular_128_to_256_relative_change"] < 0.005
            for sensitivity in row["solar_system"]
        )
        )
        for row in rows
    )
    current_cassini_empty_all_g_ext = all(
        not any(
            row["solar_system"]
            and row["solar_system"][index]["cassini_2026_pass_95"]
            for row in rows
        )
        for index in range(len(G_EXT_VALUES))
    )
    legacy_cassini_empty_primary_g_ext = not any(
        row["solar_system"]
        and row["solar_system"][1]["cassini_2014_pass_95"]
        for row in rows
    )
    robust_empty = (
        benchmark_result["all_pass"]
        and numerical_pass
        and current_cassini_empty_all_g_ext
    )

    return {
        "artifact": "registered SPARC x Cassini tanh-log joint execution",
        "registration_commit": "9c77e7be8790b80e83e2167dc924691bb0ecffaf",
        "instrument_amendment_1_commit": "e2d25f74c0f41c1bd05b6130b6ea48a229a3ba51",
        "instrument_amendment_2": (
            "Research/preregistrations/sparc_cassini_tanhlog/"
            "INSTRUMENT_AMENDMENT_2.md; contained in execution_parent_commit"
        ),
        "execution_parent_commit": git_head(),
        "inputs": {
            "profile_path": str(profile_path.relative_to(ROOT)),
            "profile_sha256": sha256_file(profile_path),
            "profile_source_analysis": profile["source_analysis"],
            "profile_data": profile["data"],
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "cassini": {
            "current_2026": {
                **CASSINI_2026,
                "two_sided_95_interval_s2": list(current_interval),
            },
            "legacy_2014": {
                **CASSINI_2014,
                "two_sided_95_interval_s2": list(legacy_interval),
            },
        },
        "instrument_validation": {
            "benchmark": benchmark_result,
            "all_mapping_and_numerical_checks_pass": numerical_pass,
            "primary_angular_order": 128,
            "angular_check_order": 256,
            "adaptive_radial_orders": list(RADIAL_ORDERS),
            "convergence_rule": (
                "two consecutive radial changes below 0.5%, followed by "
                "angular 128-to-256 change below 0.5%"
            ),
            "short_circuit_rule": (
                "Cassini is not evaluated where SPARC delta_BIC exceeds 14 "
                "under all three recorded conventions; such a point cannot "
                "belong to any primary or mandatory-sensitivity intersection"
            ),
        },
        "primary_definition": {
            "sparc": "delta_bic_profile <= 10",
            "cassini": (
                "signed Q2 inside the current two-sided 95% interval at "
                "g_ext=2.32e-10 m s^-2"
            ),
            "ambiguity_guard": (
                "also report fixed-gamma-vs-McGaugh and family-k2-vs-McGaugh "
                "BIC conventions; the verdict may not depend on choosing among them"
            ),
        },
        "summary": summary,
        "outcome_checks": {
            "current_cassini_survival_set_empty_for_all_g_ext": (
                current_cassini_empty_all_g_ext
            ),
            "legacy_cassini_survival_set_empty_at_primary_g_ext": (
                legacy_cassini_empty_primary_g_ext
            ),
            "robust_empty_intersection": robust_empty,
        },
        "registered_outcome_branch": "A" if robust_empty else "B/C/D",
        "verdict": (
            "The registered scale-universal tanh-log family is incompatible "
            "with the joint SPARC and Cassini constraints as a QUMOND "
            "interpolation family at the registered resolution."
            if robust_empty
            else "No branch-A verdict: inspect failed checks and surviving intervals."
        ),
        "scope_caveat": (
            "This closes only the scale-universal QUMOND realization. It does "
            "not test modified inertia, engineering-only companders, multi-scale "
            "or system-dependent theories, dark/hybrid models, or the "
            "Synchronism umbrella ontology."
        ),
        "aqual_sensitivity": (
            "Not independently computed. Literature transfer says AQUAL is "
            "typically at least as constrained for the benchmark families; "
            "no AQUAL-specific verdict is claimed."
        ),
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("profile", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = execute(args.profile.resolve())
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0 if result["outcome_checks"]["robust_empty_intersection"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
