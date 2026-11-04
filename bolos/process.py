r"""Collision processes and numerical helpers for cross-section data.

Examples
--------
>>> from bolos.process import Process
>>> collision = Process(
...     target="Ar",
...     kind="ELASTIC",
...     data=[[0, 1e-20], [10, 1e-20]],
... )
>>> collision.x.tolist()
[0.0, 10.0]
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.interpolate import interp1d

from .grid import FloatArray, Grid

_QUAD_X, _QUAD_W = np.polynomial.legendre.leggauss(8)


class Process:
    r"""Represent a collision process and its cross section.

    Parameters
    ----------
    target : str or None, optional
        Name of the collision target.
    kind : str or None, optional
        Process type, such as ``'ELASTIC'`` or ``'EXCITATION'``.
    data : array-like
        Two-column table containing energy in eV and cross section
        values.
    comment : str, optional
        Free-form description of the process.
    mass_ratio : float or None, optional
        Projectile-to-target mass ratio.
    product : str or None, optional
        Product species for processes that create one.
    threshold : float, optional
        Energy shift in eV used for inelastic collisions.
    weight_ratio : float or None, optional
        Statistical weight ratio associated with the process.

    Attributes
    ----------
    x, y : numpy.ndarray
        Energy and cross-section columns from ``data``.
    interp : callable
        Linear cross-section interpolator with constant endpoint
        extrapolation.
    isnull : bool
        Whether this is a regular process (always ``False`` here).

    Examples
    --------
    >>> process = Process(
    ...     target="Ar",
    ...     kind="ELASTIC",
    ...     data=[[0.0, 1.0e-20], [10.0, 2.0e-20]],
    ... )
    >>> process.x.tolist()
    [0.0, 10.0]
    """

    # IN_FACTOR is the number of electrons in the final state for an
    # inelastic electron collision.
    # ex. EXCITATION: M + e -> M* + e', IN_FACTOR = 1, SHIFT_FACTOR = 1
    # ex. IONIZATION: M + e -> M+ + 2e', IN_FACTOR = 2, SHIFT_FACTOR = 2
    IN_FACTOR = {
        "EXCITATION": 1,
        "IONIZATION": 2,
        "ATTACHMENT": 0,
        "ELASTIC": 1,
        "MOMENTUM": 1,
        "EFFECTIVE": 1,
    }

    # SHIFT_FACTOR is the energy shift factor for inelastic electron
    # collisions.
    # This factor represents the energy distribution of the outgoing
    # electron(s) after the collision.
    # Before energy = SHIFT_FACTOR * After energy + threshold
    SHIFT_FACTOR = {
        "EXCITATION": 1,
        "IONIZATION": 2,
        "ATTACHMENT": 1,
        "ELASTIC": 1,
        "MOMENTUM": 1,
        "EFFECTIVE": 1,
    }

    def __init__(
        self,
        target: str | None = None,
        kind: str | None = None,
        data: ArrayLike | None = None,
        comment: str = "",
        mass_ratio: float | None = None,
        product: str | None = None,
        threshold: float = 0,
        weight_ratio: float | None = None,
    ) -> None:
        r"""Initialize a process from a validated cross-section table.

        Parameters
        ----------
        target : str or None, optional
            Name of the collision target.
        kind : str or None, optional
            Process type, such as ``'ELASTIC'`` or ``'EXCITATION'``.
        data : array-like
            Two-column table containing energy in eV and cross section
            values.
        comment : str, optional
            Free-form description of the process.
        mass_ratio : float or None, optional
            Projectile-to-target mass ratio.
        product : str or None, optional
            Product species for processes that create one.
        threshold : float, optional
            Energy shift in eV used for inelastic collisions.
        weight_ratio : float or None, optional
            Statistical weight ratio associated with the process.

        Raises
        ------
        ValueError
            If the table is empty, not two-column, nonfinite, negative,
            or
            its energy nodes are not strictly increasing.

        Examples
        --------
        >>> collision = Process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> collision.data.shape
        (2, 2)
        """
        self.target_name = target
        # Set by Target.add_process() when this process is registered.
        self.target = None
        self.kind = kind
        self.data = np.array(data, dtype=np.float64, copy=True)
        if (
            self.data.ndim != 2
            or self.data.shape[1] != 2
            or not len(self.data)
        ):
            raise ValueError(
                "Cross-section data must be a nonempty two-column table."
            )
        if not np.all(np.isfinite(self.data)):
            raise ValueError(
                "Cross-section data must contain only finite values."
            )
        if np.any(self.data < 0):
            raise ValueError(
                "Energies and cross sections must be nonnegative."
            )
        if np.any(np.diff(self.data[:, 0]) <= 0):
            raise ValueError(
                "Cross-section energies must be strictly increasing."
            )
        self.x = self.data[:, 0]
        self.y = self.data[:, 1]
        self.comment = comment
        self.mass_ratio = mass_ratio
        self.product = product
        self.threshold = threshold
        self.weight_ratio = weight_ratio
        self.interp = padinterp(self.data)
        self.isnull = False
        self.in_factor = self.IN_FACTOR.get(self.kind)
        self.shift_factor = self.SHIFT_FACTOR.get(self.kind)
        self.cached_grid = None

    def scatterings(self, g: FloatArray, eps: FloatArray) -> FloatArray:
        r"""Calculate scattering integrals for the cached energy grid.

        Parameters
        ----------
        g : array-like
            Exponential slopes for the energy cells.
        eps : array-like
            Energy values at which the integrals are evaluated.

        Returns
        -------
        numpy.ndarray
            Scattering contribution for each selected energy cell. An
            empty
            float array is returned when no cached intervals overlap the
            grid.

        Notes
        -----
        Call :meth:`set_grid_cache` before this method. The cache must
        correspond to the same grid used to produce ``g`` and ``eps``.
        "gj" and "epsj" are the slope and energy values for the cached
        intervals.
        ".eps" and ".sigma" are the energy and cross-section values for
        the cached intervals.

        Examples
        --------
        >>> process = Process(
        ...     kind="ELASTIC",
        ...     data=[[0.0, 1.0], [2.0, 1.0]],
        ... )
        >>> process.set_grid_cache(
        ...     type("Grid", (), {"b": np.array([0.0, 1.0, 2.0])})()
        ... )
        >>> np.allclose(
        ...     process.scatterings(np.zeros(2), np.ones(2)),
        ...     [0.5, 1.5],
        ... )
        True
        """
        # self.j is the index of the grid cell that to be integrated
        # over.
        if len(self.j) == 0:
            # Avoid indexing g[self.j] when no interval overlaps the
            # grid.
            return np.array([], dtype=np.float64)

        gj = g[self.j]
        epsj = eps[self.j]
        # Integrate the linear interpolated cross section times the
        # exponential over the cached intervals.
        r = int_linexp0(
            self.eps[:, 0],
            self.eps[:, 1],
            self.sigma[:, 0],
            self.sigma[:, 1],
            gj,
            epsj,
        )
        return r

    def set_grid_cache(self, grid: Grid) -> None:
        r"""Cache intersections between grid cells and shifted cells.

        Parameters
        ----------
        grid : object
            Energy grid with an increasing ``b`` array of cell
            boundaries.

        Returns
        -------
        None
            The calculated intervals are stored on ``self``.

        Notes
        -----
        Reusing the identical grid object does not recalculate the
        cache.
        The cached arrays ``i``, ``j``, ``sigma``, and ``eps`` describe
        overlapping intervals and interpolated cross sections.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> collision = Process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> collision.set_grid_cache(LinearGrid(0, 10, 5))
        >>> collision.eps.shape[1]
        2
        """
        # Create an array with matching
        # rows ([i], [j], [eps1, eps2], [sigma1, sigma2])
        # containing the overlap between shifted cell i and cell j.
        # There may
        # be multiple rows for a given i, j if an interpolation point
        # for
        # sigma falls inside the interval.

        if self.cached_grid is grid:
            return

        self.cached_grid = grid

        eps1 = self.shift_factor * grid.b + self.threshold
        eps1 = np.maximum(eps1, grid.b[0] + 1e-9)
        eps1 = np.minimum(eps1, grid.b[-1] - 1e-9)

        fltb = np.logical_and(grid.b >= eps1[0], grid.b <= eps1[-1])
        fltx = np.logical_and(self.x >= eps1[0], self.x <= eps1[-1])
        nodes = np.unique(np.r_[eps1, grid.b[fltb], self.x[fltx]])

        sigma0 = self.interp(nodes)
        self.j = np.searchsorted(grid.b, nodes[1:]) - 1
        self.i = np.searchsorted(eps1, nodes[1:]) - 1
        self.sigma = np.c_[sigma0[:-1], sigma0[1:]]
        self.eps = np.c_[nodes[:-1], nodes[1:]]

    def __str__(self):
        r"""Return a concise human-readable description of the process.

        Returns
        -------
        str
            Process kind, target name, and optional product.

        Examples
        --------
        >>> collision = Process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> str(collision)
        '{ELASTIC: Ar }'
        """
        return "{%s: %s %s}" % (
            self.kind,
            self.target_name,
            "-> " + self.product if self.product else "",
        )


class NullProcess(Process):
    r"""Represent a process with an identically zero cross section.

    Parameters
    ----------
    target : str
        Name of the collision target.
    kind : str
        Process type.

    Examples
    --------
    >>> process = NullProcess("Ar", "ELASTIC")
    >>> process.isnull
    True
    >>> process.interp(np.array([0.0, 1.0])).tolist()
    [0.0, 0.0]
    """

    def __init__(self, target: str, kind: str) -> None:
        r"""Initialize a zero-cross-section process.

        Parameters
        ----------
        target : str
            Name of the collision target.
        kind : str
            Process type.

        Examples
        --------
        >>> collision = NullProcess("Ar", "ELASTIC")
        >>> collision.isnull
        True
        """
        self.data = np.empty((0, 2))
        self.interp = lambda x: np.zeros_like(x)
        self.target_name = target
        self.kind = kind
        self.isnull = True

        self.comment = None
        self.mass_ratio = None
        self.product = None
        self.threshold = None
        self.weight_ratio = None

        self.x = np.array([])
        self.y = np.array([])

    def __str__(self):
        r"""Return the fixed display name for a null process.

        Returns
        -------
        str
            ``"{NULL}"``.

        Examples
        --------
        >>> str(NullProcess("Ar", "ELASTIC"))
        '{NULL}'
        """
        return "{NULL}"


def padinterp(data: FloatArray) -> interp1d:
    r"""Create a linear interpolator with constant extrapolation.

    Parameters
    ----------
    data : numpy.ndarray
        Two-column array containing nonnegative energies and cross
        sections.
        The first column must be sorted in ascending order.

    Returns
    -------
    scipy.interpolate.interp1d
        Linear interpolator. Within 0 to 1e8 eV, values outside the
        tabulated energies are extended from the respective endpoint.
        Values outside the padded range raise ValueError.

    Examples
    --------
    >>> interpolator = padinterp(
    ...     np.array([[1.0, 2.0], [3.0, 4.0]])
    ... )
    >>> interpolator(np.array([0.0, 2.0, 4.0])).tolist()
    [2.0, 3.0, 4.0]
    """
    if data[0, 0] > 0:
        x = np.r_[0.0, data[:, 0], 1e8]
        y = np.r_[data[0, 1], data[:, 1], data[-1, 1]]
    else:
        x = np.r_[data[:, 0], 1e8]
        y = np.r_[data[:, 1], data[-1, 1]]

    return interp1d(x, y, kind="linear")


def int_linexp0(
    a: ArrayLike,
    b: ArrayLike,
    u0: ArrayLike,
    u1: ArrayLike,
    g: ArrayLike,
    x0: ArrayLike,
) -> FloatArray:
    r"""Integrate a linearly interpolated function times an exponential.

    The integral is

    .. math::

       \int_a^b u(x) \exp(g (x_0 - x)) x\,dx,

    where ``u`` is linear and has endpoint values ``u0`` and ``u1``.
    Inputs are broadcast according to NumPy rules.

    Parameters
    ----------
    a, b : array-like
        Lower and upper integration bounds.
    u0, u1 : array-like
        Values of the linear function at ``a`` and ``b``.
    g : array-like
        Exponential slopes. Non-finite slopes contribute zero.
    x0 : array-like
        Reference energy in the exponential.

    Returns
    -------
    numpy.ndarray
        Broadcasted integral values.

    Raises
    ------
    FloatingPointError
        If a finite slope produces a non-finite integral.

    Examples
    --------
    >>> int_linexp0(0.0, 1.0, 1.0, 1.0, 0.0, 0.0)
    array(0.5)
    """
    # Since u(x) is linear, we calculate separately the coefficients
    # of degree 0 and 1 which, after multiplying by the x in the
    # integrand
    # correspond to 1 and 2

    # The expressions involve the following exponentials that are
    # problematic:
    # expa = np.exp(g * (-a + x0))
    # expb = np.exp(g * (-b + x0))
    # The problems come with small g: in that case, the exp() rounds to
    # 1
    # and neglects the order 1 and 2 terms that are required to cancel
    # the
    # 1/g**2 and 1/g**3 below.  The solution is to rewrite the
    # expressions
    # as functions of expm1(x) = exp(x) - 1, which is guaranteed to be
    # accurate
    # even for small x.
    a, b, u0, u1, g, x0 = np.broadcast_arrays(a, b, u0, u1, g, x0)
    finite = np.isfinite(g)
    small = finite & (
        np.abs(g) * np.maximum(np.abs(a - x0), np.abs(b - x0)) < 0.1
    )
    # Do not evaluate singular formulae for a flat distribution.
    # Nonfinite
    # slopes retain the legacy zero contribution; the solver diagnoses
    # them.
    work_g = np.where(
        finite & ~small,
        g,
        1 / (1 + np.maximum(np.abs(a - x0), np.abs(b - x0))),
    )
    expm1a = np.expm1(work_g * (-a + x0))
    expm1b = np.expm1(work_g * (-b + x0))

    ag = a * work_g
    bg = b * work_g

    ag1 = ag + 1
    bg1 = bg + 1

    g2 = work_g * work_g
    g3 = g2 * work_g

    # These are the expressions as functions of expa/expb
    # A1 = (  expa * ag1
    #        - expb * bg1) / g2

    # A2 = (expa * (2 * ag1 + ag * ag) -
    #       expb * (2 * bg1 + bg * bg)) / g3

    A1 = (expm1a * ag1 + ag - expm1b * bg1 - bg) / g2

    A2 = (
        expm1a * (2 * ag1 + ag * ag)
        + ag * (ag + 2)
        - expm1b * (2 * bg1 + bg * bg)
        - bg * (bg + 2)
    ) / g3

    # The factors multiplying each coefficient can be obtained by
    # the interpolation formula of u(x) = c0 + c1 * x
    c0 = (a * u1 - b * u0) / (a - b)
    c1 = (u0 - u1) / (a - b)

    r = c0 * A1 + c1 * A2

    if np.any(small):
        # Integrate in a local coordinate to avoid cancellation at g ->
        # 0.
        # For g=0 the integrand is a quadratic polynomial, integrated
        # exactly.
        t = (_QUAD_X + 1) / 2
        width = (b - a)[small, None]
        x = a[small, None] + width * t
        u = u0[small, None] * (1 - t) + u1[small, None] * t
        values = (
            width[:, 0]
            / 2
            * np.sum(
                _QUAD_W
                * u
                * x
                * np.exp(g[small, None] * (x0[small, None] - x)),
                axis=1,
            )
        )
        r = np.array(r, copy=True)
        r[small] = values
    r = np.where(finite, r, 0.0)
    if not np.all(np.isfinite(r)):
        raise FloatingPointError(
            "Nonfinite collision integral for a finite EEDF slope; "
            "check the EEDF and energy grid."
        )
    return r
