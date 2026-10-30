"""Reproducible numerical checks using the sibling Phelps data directory."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import types
import warnings

import numpy as np
import scipy
from scipy.integrate import quad

from bolos import grid, parser, solver
from bolos._compat import simps
from bolos.process import int_linexp0

ROOT = Path(__file__).resolve().parent


def original_solver(source):
    """Execute original equations with only Python 2 runtime compatibility shims.

    This is a reference run under Python 3, not an actual Python 2 execution.
    """
    process_module = types.ModuleType('process')
    exec(compile((source / 'bolos/process.py').read_text(), 'original/process.py', 'exec'),
         process_module.__dict__)
    target_module = types.ModuleType('target')
    previous = {key: sys.modules.get(key) for key in ('process', 'target')}
    try:
        sys.modules['process'] = process_module
        exec(compile((source / 'bolos/target.py').read_text(), 'original/target.py', 'exec'),
             target_module.__dict__)
        sys.modules['target'] = target_module
        module = types.ModuleType('original_solver')
        module.xrange, module.unicode = range, str
        exec(compile((source / 'bolos/solver.py').read_text(), 'original/solver.py', 'exec'),
             module.__dict__)
        module.integrate = types.SimpleNamespace(simps=simps)
        return module
    finally:
        for key, value in previous.items():
            if value is None:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = value


class LegacyDict(dict):
    def iteritems(self):
        return iter(self.items())


def solve(module, processes, en):
    coarse = grid.LinearGrid(0, 60., 500)
    b = module.BoltzmannSolver(coarse)
    if module is not solver:
        b.target = LegacyDict(b.target)
    loaded = b.load_collisions(processes)
    b.set_density('N2', 0.8)
    b.set_density('O2', 0.2)
    b.kT = 300 * module.KB / module.ELECTRONVOLT
    b.EN = en * module.TOWNSEND
    b.init()
    f = b.converge(b.maxwell(2.), maxn=1000, rtol=1e-5)
    b.grid = grid.QuadraticGrid(0, 15 * b.mean_energy(f), 400)
    b.init()
    f, iterations, error = b.converge(b.grid.interpolate(f, coarse),
                                     maxn=1000, rtol=1e-6, full=True)
    rates = np.array([b.rate(f, p) for p in loaded])
    # Test both the string-based and object-based public rate interfaces.
    p = next(p for p in loaded if p.product)
    np.testing.assert_allclose(b.rate(f, str(p.target_name) + ' -> ' + p.product),
                               b.rate(f, b.search(p.target_name, p.product)))
    assert np.isfinite(f).all() and np.min(f) >= -1e-12
    assert np.isfinite(rates).all() and np.min(rates) >= -1e-30
    np.testing.assert_allclose(b._norm(f), 1., rtol=1e-12)
    metrics = dict(mean_energy_eV=float(b.mean_energy(f)),
                   mobility_N=float(b.mobility(f)), diffusion_N=float(b.diffusion(f)))
    assert all(np.isfinite(value) and value > 0 for value in metrics.values())
    return f, rates, metrics, iterations, error


def read_bolsig(path):
    records, current = [], None
    fields = {'Mean energy (eV)': 'mean_energy_eV',
              'Mobility *N (1/m/V/s)': 'mobility_N',
              'Diffusion coefficient *N (1/m/s)': 'diffusion_N'}
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('Electric field / N (Td)'):
            current = {'EN_Td': float(line.split()[-1])}
            records.append(current)
        elif current is not None:
            for prefix, field in fields.items():
                if line.startswith(prefix):
                    current[field] = float(line.split()[-1])
    return records


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--data', type=Path, default=ROOT.parent / 'Phelps')
    args = cli.parse_args()
    # Independent collision-integral checks include flat and nearly flat EEDFs.
    collision_checks = 0
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for slope in (0., 1e-12, -1e-8, 0.01, -0.01, 1., -10.):
            actual = int_linexp0(np.array([1.]), np.array([1.3]),
                                 np.array([2.]), np.array([3.]),
                                 np.array([slope]), np.array([1.15]))[0]
            expected = quad(lambda x: (2+(x-1)/0.3)*x*np.exp(slope*(1.15-x)),
                            1., 1.3, epsabs=1e-12)[0]
            np.testing.assert_allclose(actual, expected, rtol=1e-10)
            collision_checks += 1
    # An independent per-panel integral checks legacy odd/even and uneven grids.
    for n in (2, 3, 4, 5, 200, 201):
        x = np.linspace(0.1, 2., n) ** 1.5
        y = np.exp(-x)
        def panels(start, stop):
            value = 0.
            for i in range(start, stop - 2, 2):
                polynomial = np.polyfit(x[i:i+3], y[i:i+3], 2)
                value += quad(lambda z: np.polyval(polynomial, z), x[i], x[i+2])[0]
            return value
        if n % 2:
            expected = panels(0, n)
        else:
            expected = 0.5 * (panels(0, n-1) + panels(1, n)
                              + (x[1]-x[0])*(y[1]+y[0])/2
                              + (x[-1]-x[-2])*(y[-1]+y[-2])/2)
        np.testing.assert_allclose(simps(y, x=x), expected, rtol=1e-12)

    datafile = args.data / 'Phelps_N2-O2_20261026.txt'
    with datafile.open(encoding='utf-8') as fp:
        processes = parser.parse(fp)
    assert len(processes) == 43
    reference = original_solver(ROOT.parent / 'bolos-master')
    bolsig = read_bolsig(args.data / 'N2-O2_80-20_2000Td.dat')
    report = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                  cross_section_sha256=hashlib.sha256(datafile.read_bytes()).hexdigest(),
                  process_count=len(processes), composition={'N2': 0.8, 'O2': 0.2},
                  gas_temperature_K=300, quadrature_checks=6,
                  collision_integral_checks=collision_checks, cases=[])
    # Stability fixes intentionally differ from the original in tiny tails and
    # flat-slope integrals. Compare meaningful outputs, not relative tail noise.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        for requested in (1, 10, 100, 1000, 2000):
            record = min(bolsig, key=lambda r: abs(r['EN_Td'] - requested))
            en = record['EN_Td']
            f, rates, metrics, iterations, error = solve(solver, processes, en)
            new_warnings = sorted({str(w.message) for w in caught})
            assert not new_warnings, new_warnings
            with warnings.catch_warnings(record=True) as legacy_warnings:
                warnings.simplefilter('always')
                old_f, old_rates, old_metrics, _, _ = solve(reference, processes, en)
            np.testing.assert_allclose(f, old_f, rtol=1e-6, atol=max(old_f)*1e-7)
            np.testing.assert_allclose(rates, old_rates, rtol=1e-4, atol=1e-28)
            np.testing.assert_allclose(list(metrics.values()), list(old_metrics.values()), rtol=1e-7)
            case = dict(EN_Td=en, iterations=iterations, convergence_error=float(error),
                        metrics=metrics, reference_max_abs_eedf_difference=float(np.max(abs(f-old_f))),
                        bolsig=record, bolsig_relative_difference={
                            k: (v-record[k])/record[k] for k, v in metrics.items()},
                        reference_metric_relative_difference={k: (v-old_metrics[k])/old_metrics[k] for k,v in metrics.items()},
                        reference_rate_max_relative_difference=float(np.max(abs(rates-old_rates)/np.maximum(abs(old_rates), 1e-28))),
                        runtime_warnings=new_warnings,
                        legacy_runtime_warnings=sorted({str(w.message) for w in legacy_warnings}),
                        rates_m3_per_s=rates.tolist())
            report['cases'].append(case)
            print(f'{en:g} Td: mean energy={metrics["mean_energy_eV"]:.8g} eV; reference OK', flush=True)
        report['runtime_warnings'] = sorted({str(w.message) for w in caught})
    output = ROOT / 'validation_results.json'
    output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(f'PASS: {output}')


if __name__ == '__main__':
    main()
