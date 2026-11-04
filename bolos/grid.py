r"""Energy grids, coordinate mappings, and EEDF interpolation.

Grid stores cell geometry. MappedGrid builds it from an invertible
mapping;
AutomaticGrid builds it from an EEDF. Energies are in eV. Short
attribute
names are retained for compatibility with the solver.

Examples
--------
>>> from bolos.grid import LinearGrid
>>> LinearGrid(0, 10, 5).n
5
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray: TypeAlias = NDArray[np.float64]
EnergyValue: TypeAlias = float | FloatArray


class Grid:
    r"""Cell geometry shared by all energy grids.

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
        Differences of boundary energies raised to 3/2, for EEDF
        integrals.

    Notes
    -----
    Use a concrete grid such as ``LinearGrid(x0, x1, n)`` for the usual
    API.
    To add a coordinate mapping, subclass MappedGrid and implement f and
    finv.
    Older custom subclasses of Grid using the (x0, x1, n) constructor
    should
    likewise inherit from MappedGrid instead.

    Examples
    --------
    >>> energy_grid = Grid([0, 1, 3])
    >>> energy_grid.n
    2
    """

    b: FloatArray  # energy of cell boundaries
    c: FloatArray  # energy of cell centers
    d: FloatArray  # cell widths
    d32: FloatArray  # differences of boundary energies raised to 3/2
    n: int  # number of cells
    x0: float  # minimum energy
    x1: float  # maximum energy
    delta: float  # x1 - x0

    def __init__(self, boundaries: ArrayLike) -> None:
        r"""Create a grid from strictly increasing energy boundaries.

        Parameters
        ----------
        boundaries : array-like of shape (n + 1,)
            Finite, nonnegative cell boundaries in eV. At least two are
            required.

        Raises
        ------
        ValueError
            If boundaries are not one-dimensional, finite, nonnegative,
            and
            strictly increasing, or fewer than two boundaries are
            supplied.

        Examples
        --------
        >>> grid = Grid([0, 1, 3])
        >>> grid.c.tolist()
        [0.5, 2.0]
        """
        self.b = np.array(boundaries, dtype=np.float64, copy=True)
        if self.b.ndim != 1 or self.b.size < 2:
            raise ValueError("A grid requires at least two 1-D boundaries.")
        if not np.all(np.isfinite(self.b)) or np.any(self.b < 0):
            raise ValueError(
                "Energy boundaries must be finite and nonnegative."
            )
        if np.any(np.diff(self.b) <= 0):
            raise ValueError("Energy boundaries must be strictly increasing.")

        self.n = self.b.size - 1  # number of cells
        self.x0 = float(self.b[0])
        self.x1 = float(self.b[-1])
        self.delta = self.x1 - self.x0
        self.c = 0.5 * (self.b[1:] + self.b[:-1])
        self.d = np.diff(self.b)
        self.d32 = self.b[1:] ** 1.5 - self.b[:-1] ** 1.5

    def interpolate(self, f: ArrayLike, other: Grid) -> FloatArray:
        r"""Interpolate an EEDF onto this grid's centers.

        Parameters
        ----------
        f : array-like of shape (other.n,)
            Source EEDF values on the other grid's cell centers.
        other : Grid
            Source grid.

        Returns
        -------
        numpy.ndarray of shape (n,)
            Linearly interpolated values. Endpoint samples extend to the
            source
            boundaries, and values outside those boundaries are zero.

        Raises
        ------
        ValueError
            If the EEDF does not have one value per source cell.

        Examples
        --------
        >>> source = LinearGrid(0, 4, 2)
        >>> destination = LinearGrid(0, 4, 4)
        >>> destination.interpolate([1, 3], source).tolist()
        [1.0, 1.5, 2.5, 3.0]
        """
        values = np.asarray(f, dtype=np.float64)
        if values.shape != (other.n,):
            raise ValueError("The EEDF must have one value per source cell.")
        return np.interp(
            self.c,
            np.r_[other.x0, other.c, other.x1],
            np.r_[values[0], values, values[-1]],
            left=0.0,
            right=0.0,
        )

    def cell(self, x: float) -> int:
        r"""Return the cell index containing an energy.

        Parameters
        ----------
        x : float
            Energy in eV. The final boundary belongs to the last cell.

        Returns
        -------
        int
            Cell index between zero and ``n - 1``. Internal boundaries
            belong
            to the cell on their right.

        Raises
        ------
        ValueError
            If the energy is nonfinite or outside the grid's range.

        Examples
        --------
        >>> grid = Grid([0, 1, 3])
        >>> grid.cell(1), grid.cell(3)
        (1, 1)
        """
        if not np.isfinite(x) or not self.x0 <= x <= self.x1:
            raise ValueError("Energy is outside the grid.")
        index = int(np.searchsorted(self.b, x, side="right")) - 1
        # Clamp to the valid range [0, n-1] to handle x == x1.
        return max(0, min(index, self.n - 1))


class MappedGrid(Grid, ABC):
    r"""A grid uniformly spaced in a transformed energy coordinate.

    Parameters
    ----------
    x0, x1 : float
        Finite minimum and maximum energies in eV, with 0 <= x0 < x1.
    n : int
        Positive number of cells.

    Notes
    -----
    f maps energy to the uniform coordinate; finv maps it back. Both
    must
    support scalar floats and NumPy arrays. Subclasses with mapping
    parameters
    must set them before calling this initializer.

    Examples
    --------
    >>> energy_grid = LinearGrid(0, 10, 5)
    >>> energy_grid.n
    5
    """

    # Transformed energy boundaries
    fx0: float
    fx1: float
    df: float

    def __init__(self, x0: float, x1: float, n: int) -> None:
        r"""Initialize a grid spaced uniformly in mapped coordinates.

        Parameters
        ----------
        x0, x1 : float
            Finite minimum and maximum energies in eV, with 0 <= x0 <
            x1.
        n : int
            Positive number of cells.

        Raises
        ------
        ValueError
            If the range, cell count, mapping parameter, or resulting
            boundaries
            are invalid.

        Examples
        --------
        >>> grid = LinearGrid(0, 10, 5)
        >>> grid.n
        5
        """
        if not np.isfinite(x0) or not np.isfinite(x1) or not 0 <= x0 < x1:
            raise ValueError("Require finite energies with 0 <= x0 < x1.")
        if n < 1:
            raise ValueError("The number of cells must be positive.")

        # Mapping functions may use the energy range (e.g.
        # GeometricGrid).
        self.x0, self.x1 = x0, x1
        self.delta = x1 - x0
        self.fx0 = float(self.f(x0))
        self.fx1 = float(self.f(x1))
        mapped_boundaries = np.linspace(self.fx0, self.fx1, n + 1)
        self.df = float(mapped_boundaries[1] - mapped_boundaries[0])
        boundaries = np.asarray(self.finv(mapped_boundaries), dtype=np.float64)
        # Preserve existing grid values, except tiny negative roundoff
        # at zero.
        if (
            x0 == 0
            and -np.finfo(np.float64).eps * max(1.0, x1) <= boundaries[0] < 0
        ):
            boundaries[0] = 0
        super().__init__(boundaries)

        # Preserve the requested range rather than roundoff in
        # finv(f(x)).
        self.x0, self.x1 = x0, x1
        self.delta = x1 - x0

    @abstractmethod
    def f(self, x: EnergyValue) -> EnergyValue:
        r"""Map energy into the uniform coordinate.

        Parameters
        ----------
        x : float or numpy.ndarray
            Energy in eV.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = LinearGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        raise NotImplementedError

    @abstractmethod
    def finv(self, w: EnergyValue) -> EnergyValue:
        r"""Map the uniform coordinate back to energy.

        Parameters
        ----------
        w : float or numpy.ndarray
            Uniform mapping coordinate.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = LinearGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        raise NotImplementedError


class LinearGrid(MappedGrid):
    r"""A grid with constant cell widths.

    Parameters
    ----------
    x0, x1 : float
        Finite minimum and maximum energies in eV, with 0 <= x0 < x1.
    n : int
        Positive number of cells.

    Notes
    -----
    Cell boundaries:
        x_i = x0 + i * delta / n
        for i = 0, ..., n.

    Cell widths:
        dx_i = delta / n
        for i = 0, ..., n-1.

    where delta = x1 - x0.

    Examples
    --------
    >>> energy_grid = LinearGrid(0, 10, 5)
    >>> energy_grid.n
    5
    """

    def f(self, x: EnergyValue) -> EnergyValue:
        r"""Map energy into the uniform coordinate.

        Parameters
        ----------
        x : float or numpy.ndarray
            Energy in eV.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = LinearGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return x

    def finv(self, w: EnergyValue) -> EnergyValue:
        r"""Map the uniform coordinate back to energy.

        Parameters
        ----------
        w : float or numpy.ndarray
            Uniform mapping coordinate.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = LinearGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return w


class QuadraticGrid(MappedGrid):
    r"""A grid whose cell widths increase linearly with cell index.

    Parameters
    ----------
    x0, x1 : float
        Finite minimum and maximum energies in eV, with 0 <= x0 < x1.
    n : int
        Positive number of cells.

    Notes
    -----
    Cell boundaries:
        x_i = x0 + (i/n)**2 * delta
        for i = 0, ..., n.

    Cell widths:
        dx_i = (2*i + 1) * delta / n**2
        for i = 0, ..., n-1.

    where delta = x1 - x0.

    Examples
    --------
    >>> energy_grid = QuadraticGrid(0, 10, 5)
    >>> energy_grid.n
    5
    """

    def f(self, x: EnergyValue) -> EnergyValue:
        r"""Map energy into the uniform coordinate.

        Parameters
        ----------
        x : float or numpy.ndarray
            Energy in eV.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = QuadraticGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return np.sqrt(x - self.x0)

    def finv(self, w: EnergyValue) -> EnergyValue:
        r"""Map the uniform coordinate back to energy.

        Parameters
        ----------
        w : float or numpy.ndarray
            Uniform mapping coordinate.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = QuadraticGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return w**2 + self.x0


class GeometricGrid(MappedGrid):
    r"""A grid with geometrically varying cell widths.

    Parameters
    ----------
    x0, x1 : float
        Finite minimum and maximum energies in eV, with 0 <= x0 < x1.
    n : int
        Positive number of cells.
    r : float, optional
        Positive cell-width ratio, different from one. Default is 1.1.

    Notes
    -----
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

    Examples
    --------
    >>> energy_grid = GeometricGrid(0, 10, 5)
    >>> energy_grid.n
    5
    """

    def __init__(self, x0: float, x1: float, n: int, r: float = 1.1) -> None:
        r"""Initialize a grid spaced uniformly in mapped coordinates.

        Parameters
        ----------
        x0, x1 : float
            Finite minimum and maximum energies in eV, with 0 <= x0 <
            x1.
        n : int
            Positive number of cells.
        r : float, optional
            Positive cell-width ratio, different from one. Default is
            1.1.

        Raises
        ------
        ValueError
            If the range, cell count, mapping parameter, or resulting
            boundaries
            are invalid.

        Examples
        --------
        >>> grid = GeometricGrid(0, 10, 5)
        >>> grid.n
        5
        """
        if not np.isfinite(r) or r <= 0 or r == 1:
            raise ValueError("The width ratio r must be positive and not one.")
        self.r = r
        self.logr = float(np.log(r))
        self.rn_minus_1 = float(np.exp(n * self.logr) - 1)
        super().__init__(x0, x1, n)

    def f(self, x: EnergyValue) -> EnergyValue:
        r"""Map energy into the uniform coordinate.

        Parameters
        ----------
        x : float or numpy.ndarray
            Energy in eV.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = GeometricGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return (
            np.log(1 + (x - self.x0) * self.rn_minus_1 / self.delta)
            / self.logr
        )

    def finv(self, w: EnergyValue) -> EnergyValue:
        r"""Map the uniform coordinate back to energy.

        Parameters
        ----------
        w : float or numpy.ndarray
            Uniform mapping coordinate.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = GeometricGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return (
            self.x0
            + self.delta * (np.exp(w * self.logr) - 1) / self.rn_minus_1
        )


class LogGrid(MappedGrid):
    r"""Define a pseudo-logarithmic grid with a positive offset.

    Parameters
    ----------
    x0, x1 : float
        Finite minimum and maximum energies in eV, with 0 <= x0 < x1.
    n : int
        Positive number of cells.
    s : float, optional
        Positive logarithmic energy offset in eV. Default is 10.

    Notes
    -----
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

    Examples
    --------
    >>> energy_grid = LogGrid(0, 10, 5)
    >>> energy_grid.n
    5
    """

    def __init__(self, x0: float, x1: float, n: int, s: float = 10.0) -> None:
        r"""Initialize a grid spaced uniformly in mapped coordinates.

        Parameters
        ----------
        x0, x1 : float
            Finite minimum and maximum energies in eV, with 0 <= x0 <
            x1.
        n : int
            Positive number of cells.
        s : float, optional
            Positive logarithmic energy offset in eV. Default is 10.

        Raises
        ------
        ValueError
            If the range, cell count, mapping parameter, or resulting
            boundaries
            are invalid.

        Examples
        --------
        >>> grid = LogGrid(0, 10, 5)
        >>> grid.n
        5
        """
        if not np.isfinite(s) or s <= 0:
            raise ValueError("The logarithmic offset s must be positive.")
        self.s = s
        super().__init__(x0, x1, n)

    def f(self, x: EnergyValue) -> EnergyValue:
        r"""Map energy into the uniform coordinate.

        Parameters
        ----------
        x : float or numpy.ndarray
            Energy in eV.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = LogGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return np.log(self.s + (x - self.x0))

    def finv(self, w: EnergyValue) -> EnergyValue:
        r"""Map the uniform coordinate back to energy.

        Parameters
        ----------
        w : float or numpy.ndarray
            Uniform mapping coordinate.

        Returns
        -------
        float or numpy.ndarray
            Transformed values with the input shape.

        Examples
        --------
        >>> grid = LogGrid(0, 10, 5)
        >>> values = np.array([0.0, 1.0, 10.0])
        >>> np.allclose(grid.finv(grid.f(values)), values)
        True
        """
        return np.exp(w) - self.s + self.x0


class AutomaticGrid(Grid):
    r"""Redistribute boundaries to equalize EEDF population.

    Parameters
    ----------
    grid : Grid
        Source energy grid.
    f0 : array-like
        Nonnegative EEDF with positive population.
    delta : float, optional
        Unused compatibility parameter.

    Notes
    -----
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

    Examples
    --------
    >>> source = LinearGrid(0, 10, 5)
    >>> energy_grid = AutomaticGrid(source, np.ones(5))
    >>> energy_grid.n
    5
    """

    def __init__(self, grid: Grid, f0: ArrayLike, delta: float = 1e-4) -> None:
        r"""Construct a grid with nearly equal population per cell.

        Parameters
        ----------
        grid : Grid
            Source grid, whose cell count and energy range are
            preserved.
        f0 : array-like of shape (grid.n,)
            Finite, nonnegative EEDF with a positive total population.
        delta : float, optional
            Unused parameter retained for compatibility. Default is
            1e-4.

        Raises
        ------
        ValueError
            If the EEDF shape, values, or total population are invalid.

        Examples
        --------
        >>> source = LinearGrid(0, 10, 5)
        >>> redistributed = AutomaticGrid(source, np.ones(5))
        >>> redistributed.n
        5
        """
        values = np.asarray(f0, dtype=np.float64)
        if values.shape != (grid.n,):
            raise ValueError("The EEDF must have one value per source cell.")
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("The EEDF must be finite and nonnegative.")

        cumulative = np.r_[0.0, np.cumsum(grid.d32 * values)]
        if not np.isfinite(cumulative[-1]) or cumulative[-1] <= 0:
            raise ValueError(
                "The EEDF must have a finite, positive population."
            )
        cumulative /= cumulative[-1]

        # Zero-population cells give repeated cumulative values. Use the
        # last
        # boundary of each plateau to make the inverse interpolation
        # unique.
        keep = np.r_[np.diff(cumulative) > 0, True]
        fractions = np.linspace(0.0, 1.0, grid.n + 1)
        boundaries = np.interp(fractions, cumulative[keep], grid.b[keep])
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
    r"""Build a mapped grid by name.

    Parameters
    ----------
    kind : str
        Grid name: linear/lin, quadratic/quad, geometric/geo, or
        logarithmic/log.
    x0, x1 : float
        Minimum and maximum energies in eV.
    n : int
        Positive number of cells.
    *args : float
        Additional positional arguments for the grid constructor.
    **kwargs : float
        Additional keyword arguments, such as ``r`` or ``s``.

    Returns
    -------
    MappedGrid
        Newly constructed grid.

    Raises
    ------
    KeyError
        If the grid name is unknown.
    ValueError
        If the grid constructor rejects its arguments.

    Examples
    --------
    >>> mkgrid("lin", 0, 10, 5).b.tolist()
    [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
    """
    return _GRID_CLASSES[kind](x0, x1, n, *args, **kwargs)
