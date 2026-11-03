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

    b: FloatArray # energy of cell boundaries
    c: FloatArray # energy of cell centers
    d: FloatArray # cell widths
    d32: FloatArray # differences of boundary energies raised to 3/2
    n: int # number of cells
    x0: float # minimum energy
    x1: float # maximum energy
    delta: float # x1 - x0

    def __init__(self, boundaries: ArrayLike) -> None:
        """Create a grid from strictly increasing energy boundaries.

        Args:
            boundaries: One-dimensional sequence of finite, nonnegative
                boundary energies. At least two boundaries are required.

        Raises:
            ValueError:
            1. If the boundaries are not one-dimensional,
            2. If the boundaries contain fewer than two values, contain non-finite or negative values,
            3. If the boundaries are not strictly increasing.
        """
        self.b = np.array(boundaries, dtype=np.float64, copy=True)
        if self.b.ndim != 1 or self.b.size < 2:
            raise ValueError("A grid requires at least two 1-D boundaries.")
        if not np.all(np.isfinite(self.b)) or np.any(self.b < 0):
            raise ValueError("Energy boundaries must be finite and nonnegative.")
        if np.any(np.diff(self.b) <= 0):
            raise ValueError("Energy boundaries must be strictly increasing.")

        self.n = self.b.size - 1 # number of cells is one less than number of boundaries
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

        Args:
            f (ArrayLike): EEDF values on the source grid's centers.
            other (Grid): The source grid.
        Returns:
            FloatArray: Interpolated EEDF values on this grid's centers.
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
        """ Search for the cell containing a given energy.
        Return the cell containing x, including x1 in the final cell.

        Internal boundaries belong to the cell on their right. Out-of-range
        and nonfinite energies raise ValueError. No analytic mapping is needed.

        Args:
            x (float): Energy in eV.
        Returns:
            int: Index of the cell containing x, in the range [0, n-1].
        Raises:
            ValueError: If x is not finite or outside the grid's energy range.
        """
        if not np.isfinite(x) or not self.x0 <= x <= self.x1:
            raise ValueError("Energy is outside the grid.")
        index = int(np.searchsorted(self.b, x, side="right")) - 1
        # Clamp to the valid range [0, n-1] to handle x == x1.
        return max(0, min(index, self.n - 1))


class MappedGrid(Grid, ABC):
    """A grid uniformly spaced in a transformed energy coordinate.

    f maps energy to the uniform coordinate; finv maps it back. Both must
    support scalar floats and NumPy arrays. Subclasses with mapping parameters
    must set them before calling this initializer.
    """

    # Transformed energy boundaries
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
    """A grid with constant cell widths.

    Cell boundaries:
        x_i = x0 + i * delta / n
        for i = 0, ..., n.

    Cell widths:
        dx_i = delta / n
        for i = 0, ..., n-1.

    where delta = x1 - x0.
    """

    def f(self, x: EnergyValue) -> EnergyValue:
        return x

    def finv(self, w: EnergyValue) -> EnergyValue:
        return w


class QuadraticGrid(MappedGrid):
    """A grid whose cell widths increase linearly with cell index.

    Cell boundaries:
        x_i = x0 + (i/n)**2 * delta
        for i = 0, ..., n.

    Cell widths:
        dx_i = (2*i + 1) * delta / n**2
        for i = 0, ..., n-1.

    where delta = x1 - x0.
    """

    def f(self, x: EnergyValue) -> EnergyValue:
        return np.sqrt(x - self.x0)

    def finv(self, w: EnergyValue) -> EnergyValue:
        return w**2 + self.x0


class GeometricGrid(MappedGrid):
    """A grid with geometrically varying cell widths.

    Each cell is r times as wide as the preceding cell.

    Cell widths:
        dx_i = delta * (r - 1) * r**i / (r**n - 1)
        for i = 0, ..., n-1.

    Cell boundaries:
        x_i = x0 + delta * (r**i - 1) / (r**n - 1)
        for i = 0, ..., n.

    where delta = x1 - x0.

    The width ratio r must be positive and different from one.
    For r > 1, cell widths increase; for 0 < r < 1, they decrease.
    Use LinearGrid for equal cell widths.
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

    Coordinate mapping:
        w = ln(s + x - x0)

    Cell boundaries:
        x_i = x0 + s * ((1 + delta/s)**(i/n) - 1)
        for i = 0, ..., n.

    Cell widths:
        dx_i = s * (1 + delta/s)**(i/n) * ((1 + delta/s)**(1/n) - 1)
        for i = 0, ..., n-1.

    where delta = x1 - x0 and s > 0.

    Cell widths increase geometrically with the ratio
        r = (1 + delta/s)**(1/n).

    A larger s produces a grid closer to uniform spacing.
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
    """Redistribute grid boundaries to approximately equalize EEDF population.

    Cell weights:
        p_i = f0_i * (b_{i+1}**1.5 - b_i**1.5)
        for i = 0, ..., n-1.

    Normalized cumulative population:
        P_j = sum(p_i, i=0..j-1) / sum(p_i, i=0..n-1)
        for j = 0, ..., n.

    New cell boundaries:
        x_j = inverse_interp(P, b, j/n)
        for j = 0, ..., n.

    where b is the original grid boundary array and f0 is
    the EEDF evaluated on the original grid.

    The number of cells and the original energy range are preserved.
    Zero-population plateaus are handled during inverse interpolation.

    The parameter delta is retained for backward compatibility
    but is not used.
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
