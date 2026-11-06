# Py3BOLOS — BOLOS Translated to Python 3

**This project is a Python 3 translation and port of the original Python 2 [BOLOS (bolos)](https://github.com/aluque/bolos) code developed by Alejandro Luque.** It is a derivative project, not an official release by the original author. Here, "translation" refers to converting the source code for Python 3 compatibility, rather than simply translating the documentation into another language.

BOLOS is a Python library for solving the electron Boltzmann equation in nonthermal plasmas. It reads collision cross sections in a BOLSIG+-compatible format and calculates electron energy distributions, mean electron energy, electron mobility, diffusion coefficients, and reaction rate coefficients.

## Original Project and Changes

- Original project: [aluque/bolos](https://github.com/aluque/bolos)
- Original author: Alejandro Luque (Instituto de Astrofísica de Andalucía, IAA-CSIC)
- Original documentation: [BOLOS documentation](https://bolos.readthedocs.io/)
- License: The original project's LGPLv2 designation is retained; see below.

The main changes include Python 3 syntax, string and iterator handling, explicit relative imports, integration compatibility for recent SciPy versions, and numerical stability improvements. This project also includes Japanese and English execution guides, a notebook, and historical validation records. It retains the original physical model and numerical approach, but identical results under all conditions are not guaranteed. See [PORTING_NOTES.md](PORTING_NOTES.md) for details.

## Requirements and Installation

`setup.py` specifies Python 3.10 or later. NumPy and SciPy are required for computation. Not every combination of Python and dependency versions has been tested.

Run the commands from the directory containing this README. Create a virtual environment there. `.venv/` is listed in `.gitignore`.

macOS and Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[notebook]"
```

Windows (PowerShell):

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -U pip
python -m pip install -e ".[notebook]"
```

If Windows blocks `Activate.ps1`, use `.\.venv\Scripts\python.exe` in place of `python`. The quotes around `".[notebook]"` keep the shell from treating that extra as a filename pattern.

Choose the extra for the tools you need:

- `".[samples]"` adds Matplotlib for the plotting samples.
- `".[notebook]"` adds Matplotlib and ipykernel.
- `".[validation]"` adds the dependencies used by `validate_bolsig.ipynb`.
- `python -m pip install -e .` installs the computational library alone.

To open `run_basic.ipynb`, install JupyterLab separately or use a compatible editor. Select this virtual environment as the notebook kernel and use the repository root as the notebook's working directory. Step-by-step instructions are in [EXECUTION_GUIDE_EN.md](EXECUTION_GUIDE_EN.md) and [EXECUTION_GUIDE_js.md](EXECUTION_GUIDE_js.md).

**The import name is `bolos`, as in the original project.** A clone's directory name can differ from the package name. Check which source file is being imported with:

```console
python -c "from bolos import solver; print(solver.__file__)"
```

The output should point to `bolos/solver.py` inside this project. If you move or rename the project directory, run `python -m pip install -e .` again from its new location to update the editable installation.

If ty reports `Cannot resolve imported module bolos`, first check the runtime import using the command above. If the runtime import succeeds but static analysis fails, configure ty to use the same virtual environment as your interpreter. When this repository is the workspace root, set `root = ["."]` under `[environment]` in `ty.toml`.

## Basic Usage

Run the following example from the repository root. It uses the BOLSIG+-compatible N2/O2 collision cross-section file in `Phelps/`. Replace the input path if you use a different dataset.

```python
import scipy.constants as constants
from bolos import grid, parser, solver

with open("Phelps/Phelps_N2-O2_20261026.txt", encoding="utf-8") as stream:
    collisions = parser.parse(stream)

energy_grid = grid.LinearGrid(0.0, 100.0, 1000)
boltzmann = solver.BoltzmannSolver(energy_grid)
boltzmann.load_collisions(collisions)
boltzmann.set_density("N2", 0.8)
boltzmann.set_density("O2", 0.2)
boltzmann.kT = 300.0 * constants.k / constants.eV
boltzmann.EN = 100.0 * solver.TOWNSEND  # E/N = 100 Td
boltzmann.init()

distribution = boltzmann.converge(
    boltzmann.maxwell(2.0), maxn=1000, rtol=1e-6
)
print("Mean energy [eV]:", boltzmann.mean_energy(distribution))
print("Mobility * N [1/(m V s)]:", boltzmann.mobility(distribution))
print("Diffusion coefficient * N [1/(m s)]:", boltzmann.diffusion(distribution))
```

Adjust the gas composition, gas temperature, reduced electric field, and energy grid for your application. Check convergence and sensitivity to the grid resolution and energy range.

## Files and Validation

| File or Directory | Description |
| --- | --- |
| `bolos/` | Python 3 computational library |
| `Phelps/` | N2/O2 collision cross sections and BOLSIG+ comparison output |
| `samples/` | Computation and plotting examples |
| `run_basic.ipynb` | Calculation notebook |
| [validate_bolsig.ipynb](validate_bolsig.ipynb) | Executed EEDF, reaction-rate, and transport comparison against the included BOLSIG+ export |
| `validation/bolsig_comparison_corrected/` | Current comparison tables, figures, and provenance after the reconstruction-sign correction |
| `validation/bolsig_comparison/` | Archived comparison results before the reconstruction-sign correction |
| `validation/collision_sign_corrected/` | Independent integral checks and corrected-versus-legacy solver comparisons |
| [EXECUTION_GUIDE_EN.md](EXECUTION_GUIDE_EN.md) | Detailed English execution guide |
| [EXECUTION_GUIDE_js.md](EXECUTION_GUIDE_js.md) | Detailed Japanese execution guide |
| [PORTING_NOTES.md](PORTING_NOTES.md) | Porting changes and limitations |
| [VALIDATION.md](VALIDATION.md) / [validation_results.json](validation_results.json) | Existing validation records |
| `docs/` | Sphinx documentation based on the original project |

The historical validation covers a limited set of cases, including five electric-field conditions at approximately 1–2000 Td for N2:O2 = 80:20 at 300 K. It is not a direct comparison against execution in a Python 2 environment. `validation_results.json` records earlier validation runs and is not regenerated by the new notebook.

For a reproducible comparison using the files included here, run [validate_bolsig.ipynb](validate_bolsig.ipynb). It compares eight exported field conditions against BOLSIG+ 07/2024, including EEDFs, eight representative reaction rates, mean energy, reduced mobility, reduced diffusion, and their derived characteristic energy. It also checks mesh and cutoff sensitivity and exports all 43 rates. Install its dependencies with `python -m pip install -e ".[validation]"`. The reference is the supplied BOLSIG+ output; the notebook does not launch BOLSIG+ itself.

Some historical documentation and validation records refer to `run_phelps.ipynb`, which is not included in this distribution. Use `run_basic.ipynb` for the included notebook example.

The current branch corrects the local exponential reconstruction sign in collision integrals. See [the sign audit](docs/COLLISION_SIGN_VALIDATION.md) for the derivation and numerical impact. To refresh the comparison notebook and explicitly export its tables, run `python validation/refresh_bolsig_notebook.py --export`. Its existing commented export cell is preserved; the runner exports into `validation/bolsig_comparison_corrected/` separately from archived results.

**The original BOLOS source and the local comparison script `validate_port.py` are not included in this distribution.** The comparison script requires a separate copy of the original source and is excluded from version control via `.gitignore`. The notebook and standard calculation examples do not require the original source. Historical validation records are retained for reference.

The repository uses the following layout:

```text
Py3bolos/
├── bolos/
├── Phelps/
│   ├── Phelps_N2-O2_20261026.txt
│   └── N2-O2_80-20_2000Td.dat
├── samples/
└── run_basic.ipynb
```

Cross-section data and BOLSIG+ output in `Phelps/` are separate from the BOLOS source license. Their use and redistribution are subject to the terms of their respective providers.

## License and Acknowledgments

This project retains the **LGPLv2** designation stated in the original BOLOS README and `setup.py`. The license referenced by the original project is the **GNU Library General Public License, version 2 (June 1991)**. Its complete official English text is included in [LICENSE.txt](LICENSE.txt).

Retain the original author's copyright and attribution notices, and follow the license terms when modifying or redistributing the software. This derivative project does not imply endorsement or support by the original author. The software is provided without warranty, as specified in the license.

We acknowledge Alejandro Luque, the original author of BOLOS, and the researchers and developers whose work underlies the software. The original project cites:

G. J. M. Hagelaar and L. C. Pitchford, “Solving the Boltzmann equation to obtain electron transport coefficients and rate coefficients for fluid models,” *Plasma Sources Science and Technology*, **14** (2005), 722–733.
