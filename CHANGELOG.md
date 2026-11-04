# Changelog

## 1.0.0

### BOLOS Python 3: Modernization, reliability fixes, and documentation

Version 1.0.0 establishes the current Python 3 derivative as its release
baseline. It includes code modernization, numerical robustness improvements,
and documentation of the computational workflow.

See [RELEASE_NOTES.md](RELEASE_NOTES.md) for the publication-ready v1.0.0 notes,
including the BOLSIG+ comparison results and compatibility information.

### Improvements

- Ported the computational library to Python 3 with explicit relative imports
  and updated iterator and string handling.
- Adopted `scipy.integrate.simpson` and removed the temporary legacy
  integration compatibility layer. Historical results using the previous
  integration rule are not expected to match exactly.
- Improved collision-integral stability for flat and nearly flat distributions
  and avoided divisions at unused boundary faces.
- Removed redundant matrix construction and duplicated processing.
- Simplified reaction-rate calculations using direct accumulation.
- Improved type annotations, logging, and property definitions.
- Standardized computational-package formatting to PEP 8 and added Ruff checks.

### Fixes

- Fixed handling of empty reaction sets, including transport calculations
  without ionization or attachment.
- Fixed shared-list references after EFFECTIVE-to-ELASTIC conversion.
- Corrected package exports.
- Preserved two-dimensional cross-section tables when parsing a single row.
- Added validation of cross-section table shape, finite and nonnegative
  values, and strictly increasing energies.

### Documentation

- Added NumPy-style docstrings and executable examples throughout `bolos/`.
- Added an [English source-code guide](docs/CODE_GUIDE.md) covering all six
  computational modules, their responsibilities, and execution flow.
- Documented units, cache behavior, numerical processing, calculation examples,
  and grid refinement.
- Added a [documentation index](docs/README.md).
- Added an executed [BOLSIG+ comparison notebook](validate_bolsig.ipynb) and
  exported tables, figures, and provenance in `validation/bolsig_comparison/`.

### Validation

- The documented regression suite passed 102 tests, including executable
  docstring examples.
- Ruff and pycodestyle checks passed during the documentation and style update.
- Comparisons before and after the computational refactoring covered
  1, 10, 100, 1000, and 2000 Td for an N2/O2 mixture at 300 K.
- EEDF values were unchanged in those comparisons. Transport coefficients and
  all 43 reaction-rate coefficients agreed within the documented tolerances.
- The source-code guide's calculation and grid-refinement examples were
  executed successfully.
- The BOLSIG+ comparison notebook completed all eight code cells and compares
  eight field conditions across 1–2000 Td. Maximum discrepancies were 0.501%
  for population-weighted EEDF L1, 0.536% for mean energy/reduced mobility/
  reduced diffusion, and 13.32% for selected rates above the documented floor.

These checks record validation performed during development. The before/after
refactoring comparisons assess consistency with the preceding implementation.
The separate notebook compares BOLOS with an existing BOLSIG+ export; it does
not execute BOLSIG+ itself. See [VALIDATION.md](VALIDATION.md) and
[PORTING_NOTES.md](PORTING_NOTES.md) for the conditions and limitations.

### Requirements and compatibility

- Python 3.10 or later.
- NumPy and SciPy 1.11 or later.
- The package import name remains `bolos`.
- Cross-section tables with nonfinite values, duplicate energies, or unsorted
  energies now raise `ValueError`; energies must be strictly increasing.
- Small negative EEDF values observed at low reduced fields in the preceding
  implementation remain unchanged by the refactoring. Check convergence,
  nonnegativity, and sensitivity to the grid and energy range for each case.

### Acknowledgments

This project is an independent Python 3 derivative of BOLOS by Alejandro
Luque. The original attribution and LGPLv2 license designation are retained.
This release does not imply endorsement by the original author.
