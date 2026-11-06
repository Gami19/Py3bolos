# BOLOS Python 3 Execution Guide

> In the public distribution, cross-section data and BOLSIG+ output are
> stored in the repository-level `Phelps/` directory. The included execution
> notebook is `run_basic.ipynb`. The original source and the comparison script
> `validate_port.py` are not included in the public distribution. Historical
> validation results are retained for reference. The `run_basic.py` mentioned
> in this guide is an example program that you save locally.

This guide explains how to calculate electron energy distribution functions
(EEDFs), reaction rate coefficients, and transport coefficients using the
provided Phelps N2/O2 cross sections.

Work in the root of the cloned repository, the directory that contains
`pyproject.toml`, `bolos/`, and `Phelps/`. The procedure is the same on
Windows, macOS, and Linux. Only virtual-environment creation and activation
differ.

Run an existing sample or the notebook first. Then read Section 5 to
understand the role of each step in a complete calculation.

## Contents

1. [Files and execution methods](#1-files-and-execution-methods)
2. [Check the Python environment](#2-check-the-python-environment)
3. [Run the existing samples](#3-run-the-existing-samples)
4. [Run the notebook](#4-run-the-notebook)
5. [Basic execution program](#5-basic-execution-program)
6. [Meaning of each program section](#6-meaning-of-each-program-section)
7. [Command-line arguments](#7-command-line-arguments)
8. [Saving results and changing the field](#8-saving-results-and-changing-the-field)
9. [Points to check in the results](#9-points-to-check-in-the-results)
10. [Errors and warnings](#10-errors-and-warnings)
11. [Re-running validation and related files](#11-re-running-validation-and-related-files)

## 1. Files and execution methods

The relevant layout is:

```text
Py3bolos/                               Repository root
├─ .venv/                               Local Python virtual environment
├─ bolos/                               Computational library
├─ Phelps/
│  ├─ Phelps_N2-O2_20261026.txt         Input cross sections
│  └─ N2-O2_80-20_2000Td.dat            BOLSIG+ output
├─ samples/
│  ├─ single.py                         Single-condition calculation
│  └─ bolsig.py                         Comparison with tabular BOLSIG+ output
├─ run_basic.ipynb                      Execution notebook
└─ EXECUTION_GUIDE_EN.md                This guide
```

`bolos/solver.py` is the computational library. Normally, do not run it
directly; import it from an execution program such as `single.py`.

| Method | Inputs and conditions | Purpose |
|---|---|---|
| `samples/single.py` | Input file, `--en`, and `--temp` on the command line | Quickly calculate one condition |
| `run_basic.ipynb` | Variables changed in notebook cells | Inspect calculations, plots, and tables |
| Your own program | Python variables or command-line arguments | Add parameter scans and result storage |

`Phelps_N2-O2_20261026.txt` contains collision cross sections.
`N2-O2_80-20_2000Td.dat` is calculated BOLSIG+ output and is not used as
cross-section input.

## 2. Check the Python environment

Run the following commands in a terminal at the repository root. They are
different from code entered at the Python `>>>` prompt. Use Python 3.10 or
later. On macOS and Linux the executable is often `python3`; on Windows,
`py -3` selects a Python 3 installation.

### 2.1 Create and activate a virtual environment

macOS and Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Windows (PowerShell):

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

The prompt should start with `(.venv)`. Check the Python executable:

```bash
python -c "import sys; print(sys.executable)"
```

The output should point inside this repository's `.venv`: `.venv/bin/python`
on macOS and Linux, or `.venv\Scripts\python.exe` on Windows.

If Windows blocks `Activate.ps1` because of the execution policy, skip
activation and use that executable in place of `python`:

```powershell
.\.venv\Scripts\python.exe -c "import sys; print(sys.executable)"
```

### 2.2 Install BOLOS and check the import

From the repository root, with the virtual environment active:

```bash
python -m pip install -U pip
python -m pip install -e ".[notebook]"
python -c "from bolos import solver; print(solver.__file__)"
```

`-e` installs the project in editable mode. `[notebook]` also installs
Matplotlib and ipykernel. The quotes keep the shell from treating
`.[notebook]` as a filename pattern. Using `python -m pip` installs into the
Python environment used by the command.

Confirm that `solver.__file__` points to `bolos/solver.py` inside this
repository. Use `".[samples]"` for the plotting samples, `".[validation]"`
for the BOLSIG+ comparison notebook, or `python -m pip install -e .` for
calculations alone.

## 3. Run the existing samples

After activating the virtual environment, run the remaining commands from the
repository root. The path separator is `/`. The same form works in Windows
PowerShell.

### 3.1 Run one calculation

```bash
python samples/single.py Phelps/Phelps_N2-O2_20261026.txt --en=100 --temp=300
```

On Windows without activation, use `.\.venv\Scripts\python.exe` instead of
`python`.

| Argument | Meaning |
|---|---|
| `python` | Python used to run the program |
| `single.py` | Program to run |
| Path to `Phelps_N2-O2_20261026.txt` | Required `input` cross-section file |
| `--en=100` | Reduced electric field E/N = 100 Td |
| `--temp=300` | Gas temperature T = 300 K |

Do not type the literal word `input`; replace it with the actual file path.
The path is supplied when the program starts. These commands assume the
terminal's current directory is the repository root. Adjust the paths if you
run them elsewhere, and quote paths that contain spaces.

### 3.2 Change the conditions

Calculate at 200 Td and 300 K:

```bash
python samples/single.py Phelps/Phelps_N2-O2_20261026.txt --en=200 --temp=300
```

The N2/O2 fractions in `single.py` are fixed at N2 = 0.8 and O2 = 0.2.
Edit the corresponding lines or use the notebook to change the mixture.

Display all command-line options with:

```bash
python samples/single.py --help
```

Add `--debug` to display detailed iteration logs:

```bash
python samples/single.py Phelps/Phelps_N2-O2_20261026.txt --en=100 --temp=300 --debug
```

### 3.3 Editor run buttons

If an editor runs only `python single.py`, the required cross-section path is
missing and argparse reports:

```text
the following arguments are required: input
```

Run the command from Section 3.1 in a terminal first. The Section 5 example
sets the input path in the program, so it does not require a file argument.

### 3.4 Run `bolsig.py`

`samples/bolsig.py` requires two arguments: `input` and `bolsigdata`. It
expects a numeric table readable by `numpy.loadtxt`; the implementation uses
the second column for the field and columns four and later for reaction rates.

The supplied `N2-O2_80-20_2000Td.dat` contains explanatory lines, conditions,
and EEDF blocks, so it is not in the format expected by this sample. The
original block-format comparison script is local-only and is not included.

## 4. Run the notebook

1. Open `run_basic.ipynb` at the repository root.
2. Select the Python environment where BOLOS is installed as the kernel.
3. Set the notebook working directory to the repository root. Check that
   `Path.cwd()` points to the directory containing `bolos/` and `Phelps/`.
4. Run the cells from top to bottom.

The terminal environment and notebook kernel are selected independently.
The input path should have the following form:

```python
script_dir = Path.cwd()
input_file = script_dir / "Phelps" / "Phelps_N2-O2_20261026.txt"
```

The main notebook variables are:

```python
en_td_values = np.linspace(200.0, 1000.0, 9)
gas_temperature_k = 300.0
fractions = {"N2": 0.6, "O2": 0.4}
```

Make the fractions sum to 1. After changing conditions, re-run that cell and
all dependent cells. Notebook settings may differ from the single-condition
example and historical validation records.

## 5. Basic execution program

A BOLOS program follows this sequence:

```text
Import the library
    ↓
Set the input file and conditions
    ↓
Create the energy grid and solver
    ↓
Load cross sections
    ↓
Set composition, temperature, and field; call init()
    ↓
Create an initial EEDF and call converge()
    ↓
Optionally rebuild the grid and calculate again
    ↓
Print mean energy, transport coefficients, and rates
```

Save the following complete example as `run_basic.py` directly under the
repository root. The file is an example and is not included by default.

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as co

from bolos import grid, parser, solver


def main():
    script_dir = Path(__file__).resolve().parent
    input_file = script_dir / "Phelps" / "Phelps_N2-O2_20261026.txt"

    en_td = 100.0
    gas_temperature_k = 300.0
    fractions = {"N2": 0.8, "O2": 0.2}

    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    if (
        any(value < 0 for value in fractions.values())
        or not np.isclose(sum(fractions.values()), 1.0)
    ):
        raise ValueError("Fractions must be nonnegative and sum to 1.")

    coarse_grid = grid.LinearGrid(0.0, 60.0, 500)
    bsolver = solver.BoltzmannSolver(coarse_grid)

    with input_file.open(encoding="utf-8") as fp:
        collisions = parser.parse(fp)
    bsolver.load_collisions(collisions)

    for species, fraction in fractions.items():
        bsolver.set_density(species, fraction)
    bsolver.kT = gas_temperature_k * co.k / co.eV
    bsolver.EN = en_td * solver.TOWNSEND
    bsolver.init()

    f_initial = bsolver.maxwell(2.0)
    f_coarse = bsolver.converge(f_initial, maxn=1000, rtol=1e-5)

    coarse_mean = bsolver.mean_energy(f_coarse)
    refined_grid = grid.QuadraticGrid(0.0, 15.0 * coarse_mean, 400)
    bsolver.grid = refined_grid
    bsolver.init()

    f_guess = refined_grid.interpolate(f_coarse, coarse_grid)
    eedf, iterations, error = bsolver.converge(
        f_guess, maxn=1000, rtol=1e-6, full=True
    )

    if not np.isfinite(eedf).all() or np.min(eedf) < -1e-12:
        raise ValueError("EEDF contains non-finite or negative values.")

    print(f"E/N = {en_td:g} Td, T = {gas_temperature_k:g} K")
    print(f"Iterations = {iterations}, final error = {error:.3e}")
    print(f"Mean energy = {bsolver.mean_energy(eedf):.8g} eV")
    print(f"Mobility * N = {bsolver.mobility(eedf):.8e} 1/(m V s)")
    print(f"Diffusion * N = {bsolver.diffusion(eedf):.8e} 1/(m s)")

    print("\nReaction rate coefficients:")
    for target, process in bsolver.iter_all():
        rate = bsolver.rate(eedf, process)
        print(f"{str(process):45s} {rate:.8e} m^3/s")

    positive = eedf > 0
    plt.semilogy(bsolver.grid.c[positive], eedf[positive])
    plt.xlabel("Electron energy (eV)")
    plt.ylabel("EEDF (eV$^{-3/2}$)")
    plt.title(f"N2/O2, E/N = {en_td:g} Td")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
```

Run it from the repository root:

```bash
python run_basic.py
```

The input file is set inside the program, so no cross-section path is needed
on the command line. `__file__` refers to the program's own location, making
the input path independent of the terminal's current directory.

## 6. Meaning of each program section

### 6.1 Imports and `main()`

```python
from bolos import parser, solver, grid
```

- `parser` reads BOLSIG+-compatible cross-section files.
- `grid` creates the energy grid.
- `solver` solves for the EEDF and calculates rates and transport values.

Putting the workflow in `main()` makes the entry point explicit:

```python
if __name__ == "__main__":
    main()
```

This calls `main()` when the file is run directly, but not merely when it is
imported.

### 6.2 Energy grid

```python
coarse_grid = grid.LinearGrid(0.0, 60.0, 500)
```

This divides 0--60 eV into 500 cells. Boundaries are in `coarse_grid.b`,
centers are in `coarse_grid.c`, and the EEDF array corresponds to cell
centers. The appropriate upper energy and resolution depend on the condition.

### 6.3 Loading cross sections

```python
with input_file.open(encoding="utf-8") as fp:
    collisions = parser.parse(fp)
bsolver.load_collisions(collisions)
```

`parser.parse()` returns process dictionaries. `load_collisions()` registers
them as usable processes and targets. The supplied file contains 43 processes;
target names such as `N2` and `O2` must match the names used for composition.

### 6.4 Setting the mixture

```python
bsolver.set_density("N2", 0.8)
bsolver.set_density("O2", 0.2)
```

In this usage, `density` means molar fraction, not particle number density
in m^-3. The following is equivalent:

```python
bsolver.target["N2"].density = 0.8
bsolver.target["O2"].density = 0.2
```

### 6.5 Temperature and field units

```python
bsolver.kT = gas_temperature_k * co.k / co.eV
bsolver.EN = en_td * solver.TOWNSEND
```

| Variable | Meaning and unit |
|---|---|
| `gas_temperature_k` | Input gas temperature [K] |
| `bsolver.kT` | Gas temperature converted to energy [eV] |
| `en_td` | Input reduced field [Td] |
| `solver.TOWNSEND` | `1e-21`, conversion from Td to V m^2 |
| `bsolver.EN` | Internal reduced field [V m^2] |

For example, 100 Td is assigned as
`bsolver.EN = 100.0 * solver.TOWNSEND`.

### 6.6 When to call `init()`

```python
bsolver.init()
```

This prepares coefficients based on the mixture, temperature, field, and grid.
Call it again after changing any of those conditions, including after
assigning a new grid to `bsolver.grid`.

### 6.7 Initial distribution and iteration

```python
f_initial = bsolver.maxwell(2.0)
```

Here 2.0 is the initial Maxwell distribution parameter kT_e = 2 eV. It is
separate from the 300 K gas temperature and does not constrain the final
electron temperature or mean energy.

```python
eedf, iterations, error = bsolver.converge(
    f_guess, maxn=1000, rtol=1e-6, full=True
)
```

| Item | Meaning |
|---|---|
| `f_guess` | Initial distribution |
| `maxn` | Maximum iteration count |
| `rtol` | Convergence threshold |
| `full=True` | Also return iteration count and final error |
| `eedf` | Calculated EEDF |
| `iterations` | Actual iteration count |
| `error` | Final distribution difference |

The implementation evaluates the error as `_norm(abs(f_old - f_new))`.
Reducing `rtol` does not automatically reduce grid or cross-section errors.
If convergence is not reached, `ConvergenceError` is raised.

### 6.8 Why rebuild the grid?

The example first solves on a provisional range, then uses the mean energy:

```python
refined_grid = grid.QuadraticGrid(0.0, 15.0 * coarse_mean, 400)
```

The quadratic grid provides more resolution at low energy for this case.
When changing grids, interpolate the old EEDF:

```python
f_guess = refined_grid.interpolate(f_coarse, coarse_grid)
```

Passing an EEDF array unchanged to a grid with different cells breaks the
energy-to-distribution correspondence.

### 6.9 Meaning of output values

| Call | Result | Unit |
|---|---|---|
| `mean_energy(eedf)` | Mean electron energy | eV |
| `mobility(eedf)` | Reduced mobility, μN | 1/(m V s) |
| `diffusion(eedf)` | Reduced diffusion, DN | 1/(m s) |
| `rate(eedf, process)` | Rate coefficient, k | m^3/s |

Mobility and diffusion are multiplied by the gas number density N. If N is
known separately, divide these values by N to obtain the usual μ and D.
`rate(..., weighted=False)` does not multiply by the target molar fraction;
`weighted=True` multiplies by that fraction, not by a particle number density.

You can select a process by signature:

```python
k_ionization = bsolver.rate(eedf, "N2 -> N2^+")
```

If several processes have the same target and product, string lookup uses the
first match. To handle every process, pass process objects from `iter_all()`.

## 7. Command-line arguments

`single.py` uses `argparse`. Its basic pattern is:

```python
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("input", help="Cross-section file")
argparser.add_argument("--en", type=float, default=100.0)
argparser.add_argument("--temp", type=float, default=300.0)
args = argparser.parse_args()

input_file = args.input
en_td = args.en
gas_temperature_k = args.temp
```

For:

```text
python single.py cs.txt --en=200 --temp=350
```

the values are:

```python
args.input == "cs.txt"
args.en == 200.0
args.temp == 350.0
```

Arguments without `--` are positional and required here. `--en` and `--temp`
are optional and use their defaults when omitted. See `samples/single.py` for
the complete `parse_args()` and `main()` integration.

## 8. Saving results and changing the field

### 8.1 Save an EEDF as text

Add this after calculating `eedf` in Section 5:

```python
output_dir = script_dir / "results"
output_dir.mkdir(exist_ok=True)

np.savetxt(
    output_dir / "eedf.txt",
    np.column_stack((bsolver.grid.c, eedf)),
    header="energy_eV eedf_eV_minus_3_over_2",
)
```

The columns are energy [eV] and EEDF [eV^-3/2]. Existing files are
overwritten; include the field in the filename if you need separate results.

### 8.2 Save a plot

Add this before `plt.show()`:

```python
plt.savefig(output_dir / f"eedf_{en_td:g}Td.png", dpi=200)
```

### 8.3 Calculate several fields

After loading the cross sections and setting the mixture and temperature:

```python
for en_td in [50.0, 100.0, 200.0]:
    bsolver.EN = en_td * solver.TOWNSEND
    bsolver.init()
    initial = bsolver.maxwell(2.0)
    eedf = bsolver.converge(initial, maxn=1000, rtol=1e-6)
    print(en_td, bsolver.mean_energy(eedf))
```

Check the energy range and resolution separately for each field. For a broad
field scan, the two-stage grid calculation from Section 5 can be used for
each condition.

## 9. Points to check in the results

### 9.1 Reference values

The saved notebook result for N2:O2 = 80:20, 300 K, and 100 Td was:

| Quantity | Saved result |
|---|---:|
| Mean electron energy | 2.6272569 eV |
| μN | 1.16865409 x 10^24 1/(m V s) |
| DN | 2.56912809 x 10^24 1/(m s) |

These values depend on the notebook grid and convergence settings. `single.py`
uses the original settings, while the notebook and complete example use a
15-times-mean-energy upper limit and 400 cells, so their results can differ.

### 9.2 Convergence and grid dependence

Convergence of the iteration is not the same as convergence with respect to
the grid. Compare outputs after:

- Increasing the final grid cell count.
- Changing the upper energy limit.
- Tightening the convergence threshold.
- Changing the initial distribution.

The high-energy EEDF tail can affect high-threshold reaction rates even when
its population is small. Check the rates of interest, not only mean energy.

### 9.3 Inspect the EEDF

```python
print("Finite values only:", np.isfinite(eedf).all())
print("Minimum:", np.min(eedf))
print("BOLOS normalization:", bsolver._norm(eedf))
```

`_norm()` is an internal method, but it is useful for this validation. Its
value should be approximately 1. `np.sum(eedf)` alone does not check the
energy-weighted normalization.

## 10. Errors and warnings

### `ModuleNotFoundError: No module named 'bolos'`

Check the Python executable:

```bash
python -c "import sys; print(sys.executable)"
```

If the output is not inside this repository's `.venv`, activate the
environment from Section 2.1. If that environment still cannot import BOLOS,
install the project from the repository root:

```bash
python -m pip install -e ".[notebook]"
```

### `the following arguments are required: input`

Add the cross-section path after `single.py`, as shown in Section 3.1.

### `the following arguments are required: input, bolsigdata`

`bolsig.py` was run without its two required arguments. Use `single.py` for a
single calculation, or the notebook for the supplied BOLSIG+ comparison data.

### `FileNotFoundError`

Check the working directory and input path.

macOS and Linux:

```bash
pwd
test -f Phelps/Phelps_N2-O2_20261026.txt && echo yes
```

Windows (PowerShell):

```powershell
Get-Location
Test-Path .\Phelps\Phelps_N2-O2_20261026.txt
```

If the file is reported missing, change to the repository root or correct the
path.

### `KeyError: 'N2'` or `KeyError: 'O2'`

Check that the loaded cross sections contain those exact target names,
including capitalization and whitespace:

```python
print(list(bsolver.target))
```

### `ConvergenceError`

The convergence criterion was not met before the iteration limit. Check field
and temperature units, composition, grid/distribution correspondence, and
whether `init()` was called after changing conditions. Then adjust the grid,
initial distribution, or iteration limit. Increasing `maxn` alone is not
always sufficient.

### Negative cross section after subtracting inelastic processes

The supplied N2 data can produce negative values during the
EFFECTIVE-to-ELASTIC conversion. BOLOS clips them to zero and logs a warning.
This is a data-processing warning, not a Python syntax error, and it should be
considered when interpreting results.

### `RuntimeWarning` messages

Some unused boundary divisions use zero-flux conditions and are intentionally
not evaluated. The stable local-coordinate integration also avoids replacing
zero-over-zero expressions with arbitrary zeros.

`EEDF tail slopes undefined` is different. The legacy zero-contribution
fallback is retained, and the affected energy range and discrete population
fraction are logged. Even a small fraction can affect high-threshold rates;
check grid and cutoff sensitivity before ignoring it.

Messages such as `Energy boundary ... (safe)` and `is outside ... (safe)` are
informational. Enable them with:

```python
import logging

logging.getLogger().setLevel(logging.INFO)
```

Do not assume that every warning means calculation failure. Check convergence,
finite and nonnegative final EEDF values, normalization, and output
coefficients. Results containing `NaN` or `inf` must not be used.

### No plot appears or the program does not exit

`plt.show()` may open a window and wait until it is closed. For file-only
output, call `plt.savefig()` followed by `plt.close()`. Notebook plots usually
appear as cell output.

## 11. Re-running validation and related files

`VALIDATION.md` and `validation_results.json` are historical validation
records. The local comparison script `validate_port.py` is excluded by
`.gitignore`, and the original source is not included in the public
repository. Therefore, the original-source comparison cannot be reproduced
from the public files alone.

Historical validation covered 43 Phelps processes and N2:O2 = 80:20 at 300 K
for approximately 1, 10, 100, 1000, and 2000 Td. It uses a compatibility
port on Python 3, not a direct execution comparison in Python 2. The original
source is not required for normal calculations.

Related files:

- [Basic calculation notebook](run_basic.ipynb)
- [Single-condition sample](samples/single.py)
- [Porting notes](PORTING_NOTES.md)
- [Historical validation explanation](VALIDATION.md)
- [Detailed historical results](validation_results.json)
- [Japanese execution guide](EXECUTION_GUIDE_js.md)
