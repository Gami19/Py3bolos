"""Energy grids, coordinate mappings, and EEDF interpolation.

Grid stores cell geometry. MappedGrid builds it from an invertible mapping;
AutomaticGrid builds it from an EEDF. Energies are in eV. Short attribute
names are retained for compatibility with the solver.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import interp1d

FloatArray: TypeAlias = NDArray[np.float64]
EnergyValue: TypeAlias = float | FloatArray


class Grid:
    """Cell geometry shared by all energy grids.

    Parameters
    ----------
    boundaries : array-like
        Finite, nonnegative, strictly increasing cell boundaries in eV.

    Attributes
    ----------
    b : ndarray
        Cell boundaries, with length ``n + 1``.
    c : ndarray
        Cell centers, with length ``n``.
    d : ndarray
        Cell widths, with length ``n``.
    d32 : ndarray
        Differences of boundary energies raised to 3/2, for EEDF integrals.

    Use a concrete grid such as ``LinearGrid(x0, x1, n)`` for the usual API.
    To add a coordinate mapping, subclass MappedGrid and implement f and finv.
    Older custom subclasses of Grid using the (x0, x1, n) constructor should
    likewise inherit from MappedGrid instead.
    """

    b: FloatArray
    c: FloatArray
    d: FloatArray
    d32: FloatArray
    n: int
    x0: float
    x1: float
    delta: float

    def __init__(self, boundaries: ArrayLike) -> None:
        """Create a grid from strictly increasing energy boundaries.

        Args:
            boundaries: One-dimensional sequence of finite, nonnegative
                boundary energies. At least two boundaries are required.

        Raises:
            ValueError: If the boundaries are not one-dimensional, contain
                fewer than two values, contain non-finite or negative values,
                or are not strictly increasing.
        """
        self.b = np.array(boundaries, dtype=np.float64, copy=True)
        if self.b.ndim != 1 or self.b.size < 2:
            raise ValueError("A grid requires at least two 1-D boundaries.")
        if not np.all(np.isfinite(self.b)) or np.any(self.b < 0):
            raise ValueError("Energy boundaries must be finite and nonnegative.")
        if np.any(np.diff(self.b) <= 0):
            raise ValueError("Energy boundaries must be strictly increasing.")

        self.n = self.b.size - 1
        self.x0 = float(self.b[0])
        self.x1 = float(self.b[-1])
        self.delta = self.x1 - self.x0
        self.c = 0.5 * (self.b[1:] + self.b[:-1])
        self.d = np.diff(self.b)
        self.d32 = self.b[1:] ** 1.5 - self.b[:-1] ** 1.5

    def interpolate(self, f: ArrayLike, other: Grid) -> FloatArray:
        """Interpolate an EEDF from other onto this grid's centers.

        Endpoint samples extend to the source boundaries; outside values are
        zero. Rebuild on every call because the source EEDF or grid can change.
        """
        values = np.asarray(f, dtype=np.float64)
        if values.shape != (other.n,):
            raise ValueError("The EEDF must have one value per source cell.")
        interpolator = interp1d(
            np.r_[other.x0, other.c, other.x1],
            np.r_[values[0], values, values[-1]],
            bounds_error=False,
            fill_value=0,
        )
        return np.asarray(interpolator(self.c), dtype=np.float64)

    def cell(self, x: float) -> int:
        """Return the cell containing x, including x1 in the final cell.

        Internal boundaries belong to the cell on their right. Out-of-range
        and nonfinite energies raise ValueError. No analytic mapping is needed.
        """
        if not np.isfinite(x) or not self.x0 <= x <= self.x1:
            raise ValueError("Energy is outside the grid.")
        index = int(np.searchsorted(self.b, x, side="right")) - 1
        return max(0, min(index, self.n - 1))


class MappedGrid(Grid, ABC):
    """A grid uniformly spaced in a transformed energy coordinate.

    f maps energy to the uniform coordinate; finv maps it back. Both must
    support scalar floats and NumPy arrays. Subclasses with mapping parameters
    must set them before calling this initializer.
    """

    fx0: float
    fx1: float
    df: float

    def __init__(self, x0: float, x1: float, n: int) -> None:
        if not np.isfinite(x0) or not np.isfinite(x1) or not 0 <= x0 < x1:
            raise ValueError("Require finite energies with 0 <= x0 < x1.")
        if n < 1:
            raise ValueError("The number of cells must be positive.")

        # Mapping functions may use the energy range (e.g. GeometricGrid).
        self.x0, self.x1 = x0, x1
        self.delta = x1 - x0
        self.fx0 = float(self.f(x0))
        self.fx1 = float(self.f(x1))
        mapped_boundaries = np.linspace(self.fx0, self.fx1, n + 1)
        self.df = float(mapped_boundaries[1] - mapped_boundaries[0])
        boundaries = np.asarray(self.finv(mapped_boundaries), dtype=np.float64)
        # Preserve existing grid values, except tiny negative roundoff at zero.
        if x0 == 0 and -np.finfo(np.float64).eps * max(1.0, x1) <= boundaries[0] < 0:
            boundaries[0] = 0
        super().__init__(boundaries)

        # Preserve the requested range rather than roundoff in finv(f(x)).
        self.x0, self.x1 = x0, x1
        self.delta = x1 - x0

    @abstractmethod
    def f(self, x: EnergyValue) -> EnergyValue:
        """Map energy to the coordinate in which cells have equal widths."""
        raise NotImplementedError

    @abstractmethod
    def finv(self, w: EnergyValue) -> EnergyValue:
        """Map the uniform coordinate back to energy."""
        raise NotImplementedError


class LinearGrid(MappedGrid):
    """A grid with constant cell widths."""

    def f(self, x: EnergyValue) -> EnergyValue:
        return x

    def finv(self, w: EnergyValue) -> EnergyValue:
        return w


class QuadraticGrid(MappedGrid):
    """A grid whose cell widths increase linearly with cell index."""

    def f(self, x: EnergyValue) -> EnergyValue:
        return np.sqrt(x - self.x0)

    def finv(self, w: EnergyValue) -> EnergyValue:
        return w**2 + self.x0


class GeometricGrid(MappedGrid):
    """A grid with each cell r times as wide as the preceding one.

    r must be positive and different from one; use LinearGrid for equal widths.
    """

    def __init__(self, x0: float, x1: float, n: int, r: float = 1.1) -> None:
        if not np.isfinite(r) or r <= 0 or r == 1:
            raise ValueError("The width ratio r must be positive and not one.")
        self.r = r
        self.logr = float(np.log(r))
        self.rn_minus_1 = float(np.exp(n * self.logr) - 1)
        super().__init__(x0, x1, n)

    def f(self, x: EnergyValue) -> EnergyValue:
        return np.log(1 + (x - self.x0) * self.rn_minus_1 / self.delta) / self.logr

    def finv(self, w: EnergyValue) -> EnergyValue:
        return self.x0 + self.delta * (np.exp(w * self.logr) - 1) / self.rn_minus_1


class LogGrid(MappedGrid):
    """A pseudo-logarithmic grid with positive offset s to avoid log(0).

    Spacing is approximately logarithmic when x - x0 is much larger than s.
    """

    def __init__(self, x0: float, x1: float, n: int, s: float = 10.0) -> None:
        if not np.isfinite(s) or s <= 0:
            raise ValueError("The logarithmic offset s must be positive.")
        self.s = s
        super().__init__(x0, x1, n)

    def f(self, x: EnergyValue) -> EnergyValue:
        return np.log(self.s + (x - self.x0))

    def finv(self, w: EnergyValue) -> EnergyValue:
        return np.exp(w) - self.s + self.x0


class AutomaticGrid(Grid):
    """Redistribute boundaries to give roughly equal population per cell.

    f0 is an EEDF on grid. The cell count is preserved. delta is retained for
    call compatibility but, as in the original implementation, is unused.
    """

    def __init__(self, grid: Grid, f0: ArrayLike, delta: float = 1e-4) -> None:
        values = np.asarray(f0, dtype=np.float64)
        if values.shape != (grid.n,):
            raise ValueError("The EEDF must have one value per source cell.")
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("The EEDF must be finite and nonnegative.")

        cumulative = np.r_[0.0, np.cumsum(grid.d32 * values)]
        if not np.isfinite(cumulative[-1]) or cumulative[-1] <= 0:
            raise ValueError("The EEDF must have a finite, positive population.")
        cumulative /= cumulative[-1]

        # Zero-population cells give repeated cumulative values. Use the last
        # boundary of each plateau to make the inverse interpolation unique.
        keep = np.r_[np.diff(cumulative) > 0, True]
        inverse_population = interp1d(cumulative[keep], grid.b[keep])
        fractions = np.linspace(0.0, 1.0, grid.n + 1)
        boundaries = np.asarray(inverse_population(fractions), dtype=np.float64)
        boundaries[0], boundaries[-1] = grid.b[0], grid.b[-1]
        super().__init__(boundaries)


_GRID_CLASSES: dict[str, Callable[..., MappedGrid]] = {
    "linear": LinearGrid,
    "lin": LinearGrid,
    "quadratic": QuadraticGrid,
    "quad": QuadraticGrid,
    "geometric": GeometricGrid,
    "geo": GeometricGrid,
    "logarithmic": LogGrid,
    "log": LogGrid,
}


def mkgrid(
    kind: str, x0: float, x1: float, n: int, *args: float, **kwargs: float
) -> MappedGrid:
    """Build a mapped grid by name, forwarding optional mapping parameters.

    Names are linear/lin, quadratic/quad, geometric/geo, logarithmic/log.
    """
    return _GRID_CLASSES[kind](x0, x1, n, *args, **kwargs)
