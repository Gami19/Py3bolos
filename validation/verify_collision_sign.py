"""Verify the corrected sign against independent integrals and legacy behavior.

Run from the repository root with ``python validation/verify_collision_sign.py``.
Use ``--skip-solver`` to run only independent integral/operator checks.
"""

from contextlib import contextmanager
import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path
import platform
import re
import sys
from unittest.mock import patch
import warnings

import numpy as np
import scipy
from scipy.integrate import quad, trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bolos import grid, parser, solver
from bolos.process import Process, int_linexp0

OUT = ROOT / "validation/collision_sign_corrected"
PRODUCTION_SCATTERINGS = Process.scatterings
SELECTED_IDS = [2, 4, 21, 25, 28, 31, 37, 43]


def write_csv(name, rows):
    """Write a nonempty list of result records to the audit directory."""
    with (OUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def legacy_scatterings(self, g, eps):
    """Reproduce the old handoff by undoing the production sign conversion."""
    return PRODUCTION_SCATTERINGS(self, -g, eps)


@contextmanager
def mode_context(mode):
    """Run corrected production or an explicitly named legacy reconstruction."""
    if mode == "legacy":
        with patch.object(Process, "scatterings", legacy_scatterings):
            yield
    else:
        yield


def integral_checks():
    """Check the helper contract against independent quadrature and polynomials."""
    rng = np.random.default_rng(20261105)
    cases = [(1.0, 2.0, 1.0, 1.0, -0.5, 1.5), (0.0, 1.0, 1.0, 1.0, 0.0, 0.5)]
    for slope in [0.0, -1e-12, 1e-12, -0.01, 0.01, -0.2, 0.2, -2.0, 2.0]:
        for _ in range(20):
            a = rng.uniform(0, 5)
            b = a + rng.uniform(0.001, 2)
            cases.append(
                (
                    a,
                    b,
                    rng.uniform(0.1, 2),
                    rng.uniform(0.1, 2),
                    slope,
                    rng.uniform(a, b),
                )
            )
    rows = []
    for a, b, u0, u1, slope, center in cases:
        sigma = lambda x: u0 + (u1 - u0) * (x - a) / (b - a)
        reference = quad(
            lambda x: sigma(x) * x * np.exp(slope * (center - x)),
            a,
            b,
            epsabs=1e-12,
            epsrel=1e-12,
        )[0]
        actual = float(int_linexp0(a, b, u0, u1, slope, center))
        relative = abs(actual / reference - 1)
        assert relative < 1e-9, (a, b, slope, actual, reference)
        if slope == 0:
            width = b - a
            polynomial = u0 * (a * width + width**2 / 2)
            polynomial += (u1 - u0) * (a * width / 2 + width**2 / 3)
            np.testing.assert_allclose(actual, polynomial, rtol=1e-12)
        rows.append(
            {
                "a": a,
                "b": b,
                "sigma_a": u0,
                "sigma_b": u1,
                "helper_g": slope,
                "center": center,
                "helper_value": actual,
                "independent_quad": reference,
                "relative_error": relative,
            }
        )
    # Broadcast/array execution must obey the same helper contract.
    arrays = np.asarray(cases).T
    np.testing.assert_allclose(
        int_linexp0(*arrays), [r["independent_quad"] for r in rows], rtol=1e-9
    )
    write_csv("helper_integrals.csv", rows)
    return {
        "cases": len(rows),
        "max_relative_error": max(r["relative_error"] for r in rows),
    }


def operator_checks():
    """Test rate and matrix action on distributions with exact logarithmic slopes."""
    rows = []
    for grid_name, cls in [
        ("linear", grid.LinearGrid),
        ("quadratic", grid.QuadraticGrid),
    ]:
        for n in [8, 16, 32, 64]:
            for slope in [-0.5, 0.0, 0.5]:
                for kind, threshold in [
                    ("EXCITATION", 0.7),
                    ("IONIZATION", 1.1),
                    ("ATTACHMENT", 0.0),
                ]:
                    s = solver.BoltzmannSolver(cls(0, 8, n))
                    process = s.add_process(
                        target="X",
                        kind=kind,
                        threshold=threshold,
                        data=[[0, 1], [1.3, 1.4], [3.2, 0.8], [8, 1.8]],
                    )
                    s.set_density("X", 0.7)
                    s._slope_diagnostics_seen = set()
                    f = np.exp(slope * s.cenergy)
                    np.testing.assert_allclose(s._g(f), slope, atol=2e-13)
                    process.set_grid_cache(s.grid)
                    expected = np.zeros(n)
                    exact_intervals = []
                    for index, (a, b) in enumerate(process.eps):
                        u0, u1 = process.sigma[index]
                        value = quad(
                            lambda x: (
                                (u0 + (u1 - u0) * (x - a) / (b - a))
                                * x
                                * np.exp(slope * x)
                            ),
                            a,
                            b,
                            epsabs=1e-11,
                            epsrel=1e-12,
                        )[0]
                        exact_intervals.append(value)
                        j, i = process.j[index], process.i[index]
                        expected[i] += (
                            0.7 * solver.GAMMA * process.in_factor * value
                        )
                        expected[j] -= 0.7 * solver.GAMMA * value
                    quad_rate = solver.GAMMA * sum(exact_intervals)
                    # Closed-form moments provide a second independent oracle.
                    lower, upper = process.eps[0, 0], process.eps[-1, 1]
                    knots = np.unique(
                        np.r_[
                            lower,
                            process.x[
                                (process.x > lower) & (process.x < upper)
                            ],
                            upper,
                        ]
                    )
                    analytic_integral = 0.0
                    for a, b in zip(knots[:-1], knots[1:]):
                        u0, u1 = process.interp([a, b])
                        c1 = (u1 - u0) / (b - a)
                        c0 = u0 - c1 * a
                        if slope == 0:
                            value = (
                                c0 * (b * b - a * a) / 2
                                + c1 * (b**3 - a**3) / 3
                            )
                        else:

                            def primitive(x):
                                return np.exp(slope * x) * (
                                    c0 * (x / slope - 1 / slope**2)
                                    + c1
                                    * (
                                        x * x / slope
                                        - 2 * x / slope**2
                                        + 2 / slope**3
                                    )
                                )

                            value = primitive(b) - primitive(a)
                        analytic_integral += value
                    exact_rate = solver.GAMMA * analytic_integral
                    np.testing.assert_allclose(
                        quad_rate, exact_rate, rtol=1e-12
                    )
                    for mode in ["legacy", "corrected"]:
                        with mode_context(mode):
                            rate = s.rate(f, process)
                            action = s._PQ(f).dot(f)
                        rate_error = abs(rate / exact_rate - 1)
                        action_error = np.linalg.norm(
                            action - expected, 1
                        ) / max(0.7 * exact_rate, np.linalg.norm(expected, 1))
                        # This checks indexing/multiplicity separately from reconstruction accuracy.
                        balance_error = abs(
                            action.sum() - 0.7 * (process.in_factor - 1) * rate
                        ) / (0.7 * rate)
                        assert balance_error < 1e-12
                        if mode == "corrected" or slope == 0:
                            assert rate_error < 1e-10 and action_error < 1e-10
                        rows.append(
                            {
                                "grid": grid_name,
                                "n": n,
                                "slope": slope,
                                "kind": kind,
                                "mode": mode,
                                "exact_rate": exact_rate,
                                "computed_rate": rate,
                                "independent_quad_rate": quad_rate,
                                "rate_relative_error": rate_error,
                                "matrix_L1_relative_error": action_error,
                                "particle_balance_relative_error": balance_error,
                            }
                        )
    write_csv("manufactured_operators.csv", rows)
    return {
        mode: {
            "cases": sum(r["mode"] == mode for r in rows),
            "max_rate_relative_error": max(
                r["rate_relative_error"] for r in rows if r["mode"] == mode
            ),
            "max_matrix_L1_relative_error": max(
                r["matrix_L1_relative_error"]
                for r in rows
                if r["mode"] == mode
            ),
        }
        for mode in ["legacy", "corrected"]
    }


def reference_cases():
    """Read reference values from the original BOLSIG+ export, not prior BOLOS results."""
    number = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?"
    text = (ROOT / "Phelps/N2-O2_80-20_2000Td.dat").read_text(
        encoding="utf-8-sig"
    )
    runs = []
    for match in re.finditer(
        r"^R(\d+)[ \t]*$([\s\S]*?)(?=^R\d+[ \t]*$|\Z)", text, re.M
    ):
        block = match.group(2)

        def scalar(label):
            value = re.search(
                r"^" + re.escape(label) + r"[ \t]+(" + number + r")[ \t]*$",
                block,
                re.M,
            )
            if value is None:
                raise ValueError(label)
            return float(value.group(1).replace("D", "E"))

        rates = {
            int(m.group(1)): float(m.group(3).replace("D", "E"))
            for m in re.finditer(
                r"^C(\d+)[ \t]+(.+?)[ \t]+(" + number + r")[ \t]*$",
                block,
                re.M,
            )
        }
        eedf = []
        for line in block.split("Energy (eV) EEDF (eV-3/2) Anisotropy", 1)[
            1
        ].splitlines():
            values = line.split()
            if len(values) == 3 and all(
                re.fullmatch(number, v) for v in values
            ):
                eedf.append([float(v.replace("D", "E")) for v in values])
            elif eedf:
                break
        runs.append(
            {
                "run": int(match.group(1)),
                "field": scalar("Electric field / N (Td)"),
                "temperature": scalar("Gas temperature (K)"),
                "N2": scalar("Mole fraction N2"),
                "O2": scalar("Mole fraction O2"),
                "mean_energy_eV": scalar("Mean energy (eV)"),
                "mobility_N": scalar("Mobility *N (1/m/V/s)"),
                "diffusion_N": scalar("Diffusion coefficient *N (1/m/s)"),
                "rates": rates,
                "eedf": np.asarray(eedf),
            }
        )
    selected = [
        min(runs, key=lambda r: abs(np.log(r["field"] / x)))
        for x in [1, 3, 10, 30, 100, 300, 1000, 2000]
    ]
    for r in selected:
        assert np.isclose(r["N2"], 0.8) and np.isclose(r["O2"], 0.2)
        assert r["temperature"] == 300 and len(r["rates"]) == 43
    return selected


def solver_checks():
    """Compare both runtime variants at fixed cutoffs and two resolutions."""
    with (ROOT / "Phelps/Phelps_N2-O2_20261026.txt").open(
        encoding="utf-8"
    ) as stream:
        definitions = parser.parse(stream)
    rows, rate_rows, eedf_rows, isolated_rows, diagnostic_rows = (
        [],
        [],
        [],
        [],
        [],
    )
    for ref in reference_cases():
        s = solver.BoltzmannSolver(grid.LinearGrid(0, 100, 500))
        processes = s.load_collisions(definitions)
        s.set_density("N2", 0.8)
        s.set_density("O2", 0.2)
        s.kT = 300 * solver.KB / solver.ELECTRONVOLT
        s.EN = ref["field"] * solver.TOWNSEND
        s.init()
        with mode_context("legacy"):
            bootstrap = s.converge(s.maxwell(2), maxn=1500, rtol=1e-6)
        bootstrap_grid = s.grid
        cutoff = max(20 * ref["mean_energy_eV"], 1.25 * ref["eedf"][-1, 0])
        for n in [600, 1200]:
            energy_grid = grid.QuadraticGrid(0, cutoff, n)
            guess = energy_grid.interpolate(bootstrap, bootstrap_grid)
            fixed_f = None
            for mode in ["legacy", "corrected"]:
                with (
                    mode_context(mode),
                    warnings.catch_warnings(record=True) as caught,
                ):
                    warnings.simplefilter("always")
                    s.grid = energy_grid
                    s.init()
                    f, iters, error = s.converge(
                        guess.copy(), maxn=1500, rtol=1e-7, full=True
                    )
                    assert (
                        np.all(np.isfinite(f)) and abs(s._norm(f) - 1) < 1e-10
                    )
                    metrics = {
                        "mean_energy_eV": s.mean_energy(f),
                        "mobility_N": s.mobility(f),
                        "diffusion_N": s.diffusion(f),
                    }
                    rates = [s.rate(f, p) for p in processes]
                    # Isolate coefficient reconstruction on an identical distribution.
                    if mode == "legacy":
                        fixed_f, fixed_rates = f.copy(), rates
                    else:
                        for i, p in enumerate(processes, 1):
                            changed = s.rate(fixed_f, p)
                            base = fixed_rates[i - 1]
                            isolated_rows.append(
                                {
                                    "field_Td": ref["field"],
                                    "n": n,
                                    "collision_ID": f"C{i}",
                                    "legacy_rate": base,
                                    "corrected_same_f_rate": changed,
                                    "change_percent": 100
                                    * (changed / base - 1)
                                    if base > 1e-25
                                    else None,
                                }
                            )
                diagnostic_rows.append(
                    {
                        "field_Td": ref["field"],
                        "n": n,
                        "mode": mode,
                        "iterations": iters,
                        "iteration_error": error,
                        "norm": s._norm(f),
                        "min_f": float(f.min()),
                        "runtime_warnings": sorted(
                            {str(w.message) for w in caught}
                        ),
                    }
                )
                for key, value in metrics.items():
                    rows.append(
                        {
                            "field_Td": ref["field"],
                            "n": n,
                            "mode": mode,
                            "quantity": key,
                            "value": value,
                            "BOLSIG": ref[key],
                            "signed_difference_percent": 100
                            * (value / ref[key] - 1),
                        }
                    )
                for i, value in enumerate(rates, 1):
                    reference = ref["rates"][i]
                    rate_rows.append(
                        {
                            "field_Td": ref["field"],
                            "n": n,
                            "mode": mode,
                            "collision_ID": f"C{i}",
                            "value": value,
                            "BOLSIG": reference,
                            "signed_difference_percent": 100
                            * (value / reference - 1)
                            if reference > 1e-25
                            else None,
                        }
                    )
                x, fr = ref["eedf"][:, :2].T
                fb = np.interp(x, energy_grid.c, f)
                l1 = trapezoid(abs(fb - fr) * np.sqrt(x), x) / trapezoid(
                    fr * np.sqrt(x), x
                )
                eedf_rows.append(
                    {
                        "field_Td": ref["field"],
                        "n": n,
                        "mode": mode,
                        "weighted_L1_percent": 100 * l1,
                    }
                )
                print(
                    f"{ref['field']:.6g} Td, n={n}, {mode}: {iters} iterations",
                    flush=True,
                )
    write_csv("transport.csv", rows)
    write_csv("rates_all_43.csv", rate_rows)
    write_csv("eedf_metrics.csv", eedf_rows)
    write_csv("same_distribution_rates.csv", isolated_rows)
    (OUT / "solver_diagnostics.json").write_text(
        json.dumps(diagnostic_rows, indent=2), encoding="utf-8"
    )
    summary = {}
    for n in [600, 1200]:
        for mode in ["legacy", "corrected"]:
            selected_rates = [
                abs(r["signed_difference_percent"])
                for r in rate_rows
                if r["n"] == n
                and r["mode"] == mode
                and int(r["collision_ID"][1:]) in SELECTED_IDS
                and r["signed_difference_percent"] is not None
            ]
            summary[f"{mode}_{n}"] = {
                "max_abs_transport_difference_percent": max(
                    abs(r["signed_difference_percent"])
                    for r in rows
                    if r["n"] == n and r["mode"] == mode
                ),
                "max_selected_rate_difference_percent": max(selected_rates),
                "max_eedf_L1_percent": max(
                    r["weighted_L1_percent"]
                    for r in eedf_rows
                    if r["n"] == n and r["mode"] == mode
                ),
            }
    return summary


def main():
    """Run the numerical audit and save versioned-source fingerprints."""
    arguments = argparse.ArgumentParser(description=__doc__)
    arguments.add_argument("--skip-solver", action="store_true")
    args = arguments.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    original_level = logging.root.manager.disable
    logging.disable(logging.WARNING)
    sources = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (ROOT / "bolos").glob("*.py")
    }
    summary = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256": sources,
        "input_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                ROOT / "Phelps/Phelps_N2-O2_20261026.txt",
                ROOT / "Phelps/N2-O2_80-20_2000Td.dat",
            ]
        },
        "definition": "corrected uses production; legacy reverses g at Process.scatterings entry",
    }
    try:
        summary["helper"] = integral_checks()
        summary["manufactured"] = operator_checks()
        if not args.skip_solver:
            summary["bolsig"] = solver_checks()
    finally:
        logging.disable(original_level)
    after = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (ROOT / "bolos").glob("*.py")
    }
    assert after == sources, "The audit must not edit package sources"
    summary["package_sources_unchanged"] = True
    (
        OUT / ("helper_summary.json" if args.skip_solver else "summary.json")
    ).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
