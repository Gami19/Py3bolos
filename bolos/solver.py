r"""Solve the electron Boltzmann equation.

The :class:`BoltzmannSolver` class stores collision processes and
physical
conditions, then solves for the electron energy distribution function
(EEDF).

Examples
--------
>>> from bolos.grid import LinearGrid
>>> from bolos.solver import BoltzmannSolver
>>> BoltzmannSolver(LinearGrid(0, 10, 5)).n
5
"""

from __future__ import annotations

__docformat__ = "restructuredtext en"

import logging
from collections.abc import Iterator
from itertools import chain
from math import sqrt

import numpy as np

# Units in this module will be SI units, except energies, which are
# expressed
# in eV.
# The scipy.constants contains the recommended CODATA for all physical
# constants in SI units.
import scipy.constants as co
from scipy import sparse
from scipy.integrate import simpson
from scipy.sparse.linalg import spsolve

from .grid import FloatArray, Grid
from .process import Process
from .target import Target

logger = logging.getLogger(__name__)

GAMMA = sqrt(2 * co.elementary_charge / co.electron_mass)
TOWNSEND = 1e-21
KB = co.k
ELECTRONVOLT = co.eV


class ConvergenceError(Exception):
    r"""Signal that the EEDF iteration exceeded its iteration limit.

    Examples
    --------
    >>> error = ConvergenceError("Iteration limit reached")
    >>> isinstance(error, Exception)
    True
    """


class BoltzmannSolver:
    r"""Solve the Boltzmann equation for electrons in a gas.

    Parameters
    ----------
    grid : grid.Grid
        Energy grid on which the EEDF is evaluated.

    Attributes
    ----------
    benergy, cenergy, denergy : numpy.ndarray
        Cell boundaries, centers, and widths, respectively.
    denergy32 : numpy.ndarray
        Differences of boundary energies raised to the power 3/2.
    n : int
        Number of energy cells.
    kT : float
        Gas temperature in eV. Must be set by the user.
    EN : float
        Reduced electric field in SI units. Must be set by the user.
    target : dict
        Mapping from target names to :class:`target.Target` instances.

    Examples
    --------
    >>> from bolos.grid import LinearGrid
    >>> bsolver = BoltzmannSolver(LinearGrid(0, 60.0, 400))
    >>> bsolver.grid.n
    400
    """

    def __init__(self, grid: Grid) -> None:
        r"""Initialize a solver with an energy grid.

        Parameters
        ----------
        grid : grid.Grid
            Energy grid on which the EEDF is evaluated.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> bsolver.n
        30
        """
        self.EN = None
        self.kT = None
        self.grid = grid
        self.target: dict[str, Target] = {}

    @property
    def grid(self) -> Grid:
        r"""Return the solver's energy grid.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> bsolver.grid.n
        30
        """
        return self._grid

    @grid.setter
    def grid(self, grid: Grid) -> None:
        r"""Set the energy grid and update derived grid attributes.

        Parameters
        ----------
        grid : grid.Grid
            New energy grid.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> bsolver.grid = LinearGrid(0, 20, 40)
        >>> bsolver.n
        40
        """
        self._grid = grid
        # These are cell boundary values at i - 1/2
        self.benergy = self.grid.b
        # these are cell centers
        self.cenergy = self.grid.c
        # And these are the deltas
        self.denergy = self.grid.d
        # This is useful when integrating the growth term.
        self.denergy32 = grid.d32.copy()
        self.n = grid.n

    def set_density(self, species: str, density: float) -> None:
        r"""Set the molar fraction of a target species.

        Parameters
        ----------
        species : str
           Target species name.
        density : float
           New molar fraction.

        Raises
        ------
        KeyError
            If ``species`` has not been added to the solver.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> from bolos.target import Target
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 10))
        >>> bsolver.target["N2"] = Target("N2")
        >>> bsolver.set_density("N2", 0.8)
        >>> bsolver.target["N2"].density
        0.8
        """
        target = self.target.get(species)
        if target is None:
            raise KeyError(f"Unknown target species: {species!r}")

        target.density = density

    def load_collisions(self, dict_processes):
        r"""Load collision processes from dictionaries.

        Parameters
        ----------
        dict_processes : iterable of dict
            Process definitions accepted by :meth:`add_process`.

        Returns
        -------
        list of process.Process
            Added process instances.

        See Also
        --------
        add_process

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> definitions = [
        ...     dict(
        ...         target="Ar",
        ...         kind="ELASTIC",
        ...         mass_ratio=1e-5,
        ...         data=[[0, 1e-20], [10, 1e-20]],
        ...     )
        ... ]
        >>> len(bsolver.load_collisions(definitions))
        1
        """
        plist = [self.add_process(**p) for p in dict_processes]
        # We make sure that all targets have their elastic cross-
        # sections
        # in the form of ELASTIC cross sections (not EFFECTIVE /
        # MOMENTUM)
        for item in self.target.values():
            item.ensure_elastic()
        return plist

    def add_process(self, **kwargs):
        r"""Add a collision process to the solver.

        Parameters
        ----------
        **kwargs
            Keyword arguments passed to :class:`process.Process`,
            including
            ``target``, ``kind``, ``data``, and process metadata.

        Returns
        -------
        process.Process
            The newly added process.

        Examples
        --------
        >>> import numpy as np
        >>> from bolos import solver, grid
        >>> energy_grid = grid.LinearGrid(0, 60.0, 400)
        >>> bsolver = solver.BoltzmannSolver(energy_grid)
        >>> # This is an example cross-section that decays exponentially
        >>> energy = np.linspace(0, 10)
        >>> cross_section = 1e-20 * np.exp(-energy)
        >>> process = bsolver.add_process(
        ...     kind="EXCITATION",
        ...     target="Kriptonite",
        ...     mass_ratio=1e-5,
        ...     threshold=10,
        ...     data=np.c_[energy, cross_section],
        ... )

        See Also
        --------
        load_collisions
        """
        proc = Process(**kwargs)
        if proc.target_name is None:
            raise ValueError("Process target name must be specified.")
        try:
            target = self.target[proc.target_name]
        except KeyError:
            target = Target(proc.target_name)
            self.target[proc.target_name] = target

        target.add_process(proc)
        return proc

    def search(self, signature, product=None, first=True):
        r"""Find one or more processes matching a signature.

        Parameters
        ----------
        signature : str
            Process signature in the form ``"TARGET -> PRODUCT"``.
        product : str, optional
            If provided, interpret ``signature`` as a target name and
            search
            for this product.
        first : bool, optional
            If ``True``, return only the first match; otherwise return
            all
            matches.

        Returns
        -------
        process.Process or list of process.Process
            First matching process or all matching processes.

        Raises
        ------
        KeyError
            If the target or product is not found.


        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="EXCITATION",
        ...     product="Ar*",
        ...     threshold=1,
        ...     data=[[0, 0], [1, 0], [10, 1e-20]],
        ... )
        >>> bsolver.search("Ar -> Ar*") is collision
        True
        """
        if product is not None:
            matches = self.target[signature].by_product[product]
            if not matches:
                raise KeyError("Process %s not found" % signature)

            return matches[0] if first else matches

        t, p = [x.strip() for x in signature.split("->")]
        return self.search(t, p, first=first)

    def _iter_processes(
        self, *groups: str
    ) -> Iterator[tuple[Target, Process]]:
        r"""Iterate over selected groups for active targets.

        Parameters
        ----------
        *groups : str
            Names of Target attributes containing process lists.

        Yields
        ------
        tuple of (Target, Process)
            Active target and collision process, preserving group order.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 10))
        >>> list(bsolver._iter_processes("elastic"))
        []
        """
        for target in self.target.values():
            if target.density > 0:
                processes = (getattr(target, group) for group in groups)
                for process in chain.from_iterable(processes):
                    yield target, process

    def iter_elastic(self) -> Iterator[tuple[Target, Process]]:
        r"""Iterate over active elastic processes.

        Yields
        ------
        tuple
            A ``(target, process)`` pair.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> list(bsolver.iter_elastic())
        []
        """
        yield from self._iter_processes("elastic")

    def iter_inelastic(self) -> Iterator[tuple[Target, Process]]:
        r"""Iterate over active inelastic processes.

        Yields
        ------
        tuple
            A ``(target, process)`` pair.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> list(bsolver.iter_inelastic())
        []
        """
        yield from self._iter_processes("inelastic")

    def iter_growth(self) -> Iterator[tuple[Target, Process]]:
        r"""Iterate over active ionization and attachment processes.

        Yields
        ------
        tuple
            A ``(target, process)`` pair.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> list(bsolver.iter_growth())
        []
        """
        yield from self._iter_processes("ionization", "attachment")

    def iter_all(self) -> Iterator[tuple[Target, Process]]:
        r"""Iterate over all active elastic and inelastic processes.

        Yields
        ------
        tuple
            A ``(target, process)`` pair.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> list(bsolver.iter_all())
        []
        """
        yield from self.iter_elastic()
        yield from self.iter_inelastic()

    def iter_momentum(self) -> Iterator[tuple[Target, Process]]:
        r"""Iterate over processes used for momentum calculations.

        Yields
        ------
        tuple
            A ``(target, process)`` pair from :meth:`iter_all`.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> list(bsolver.iter_momentum())
        []
        """
        return self.iter_all()

    def init(self) -> None:
        r"""Initialize cached collision and transport coefficients.

        Call this method after changing target densities, gas
        temperature,
        electric field, or the energy grid and before iterating the
        EEDF.

        The expensive grid-dependent collision calculations are cached
        on
        each process and reused until the grid object changes.

        Notes
        -----
        Reinitialization is required whenever physical conditions
        change.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bsolver.sigma_m.shape
        (31,)
        """
        self.sigma_eps = np.zeros_like(self.benergy)
        self.sigma_m = np.zeros_like(self.benergy)
        electric_field = self.EN
        gas_temperature = self.kT
        if electric_field is None:
            raise ValueError(
                "Reduced electric field EN must be set before init()."
            )
        if gas_temperature is None:
            raise ValueError("Gas temperature kT must be set before init().")

        for target, process in self.iter_elastic():
            if target.mass_ratio is None:
                raise ValueError("Mass ratio for target '%s' is not set." % target.name)
            s = target.density * process.interp(self.benergy)
            self.sigma_eps += 2 * target.mass_ratio * s
            self.sigma_m += s
            process.set_grid_cache(self.grid)

        for target, process in self.iter_inelastic():
            self.sigma_m += target.density * process.interp(self.benergy)
            process.set_grid_cache(self.grid)

        self.W = -GAMMA * self.benergy**2 * self.sigma_eps

        # Coefficient of sigma_tilde.
        self.DA = GAMMA / 3.0 * electric_field**2 * self.benergy

        # Independent term.
        self.DB = GAMMA * gas_temperature * self.benergy**2 * self.sigma_eps

        self._slope_diagnostics_seen = set()
        for target in self.target.values():
            affected = target.clipped_elastic_range
            if target.density > 0 and affected is not None:
                overlap = (
                    self.grid.x1 > affected[0] and self.grid.x0 < affected[1]
                )
                log = logger.warning if overlap else logger.info
                log(
                    "Target %s: clipped ELASTIC range %.6g..%.6g eV %s "
                    "grid %.6g..%.6g eV. %s",
                    target.name,
                    *affected,
                    "overlaps" if overlap else "is outside",
                    self.grid.x0,
                    self.grid.x1,
                    "Check cross-section extrapolation."
                    if overlap
                    else "No direct effect on this grid (safe).",
                )
        logger.info(
            "Energy boundary faces use zero flux: unused boundary divisions "
            "are skipped (safe)."
        )

        logger.info("Solver successfully initialized/updated")

    def maxwell(self, kT: float) -> FloatArray:
        r"""Return a normalized Maxwell-Boltzmann EEDF.

        Parameters
        ----------
        kT : float
            Electron temperature in eV.

        Returns
        -------
        numpy.ndarray
            Normalized Maxwell-Boltzmann EEDF at the cell centers.

        Notes
        -----
        This is useful as an initial guess for :meth:`converge`.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> bsolver.maxwell(2.0).shape
        (30,)
        """
        return (
            2
            * np.sqrt(1 / np.pi)
            * kT ** (-3.0 / 2.0)
            * np.exp(-self.cenergy / kT)
        )

    def iterate(self, f0, delta=1e14):
        r"""Perform one normalized EEDF iteration.

        Parameters
        ----------
        f0 : array-like
            Previous EEDF.
        delta : float, optional
            Iteration parameter. Larger values can accelerate
            convergence but
            may introduce instability.

        Returns
        -------
        numpy.ndarray
            Normalized updated EEDF.

        Notes
        -----
        This is a low-level routine not intended for normal uses.  The
        standard entry point for the iterative solution of the EEDF is
        the :func:`BoltzmannSolver.converge` method.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> updated = bsolver.iterate(distribution)
        >>> bool(np.isfinite(updated).all())
        True
        """
        A, Q = self._linsystem(f0)
        f1 = spsolve(sparse.eye(self.n) + delta * A - delta * Q, f0)
        return self._normalized(f1)

    def converge(
        self, f0, maxn=100, rtol=1e-5, delta0=1e14, m=4.0, full=False, **kwargs
    ):
        r"""Iterate the EEDF until convergence is reached.

        Parameters
        ----------
        f0 : array-like
            Initial EEDF.
        maxn : int, optional
            Maximum number of iterations.
        rtol : float, optional
            Stop when the weighted L1 difference is below this
            tolerance.
        delta0 : float, optional
            Initial iteration parameter, adapted between iterations.
        m : float, optional
            Target error reduction factor for Richardson extrapolation.
        full : bool, optional
            If ``True``, return convergence information in addition to
            the
            final EEDF.

        **kwargs : object
            Additional arguments forwarded to ``iterate``. The default
            implementation accepts only the explicitly documented
            arguments.

        Returns
        -------
        numpy.ndarray or tuple
            Final EEDF, or ``(EEDF, iterations, error)`` when ``full``
            is
            ``True``.

        Notes
        -----
        If convergence is not achieved after ``maxn`` iterations, an
        exception
        of type ``ConvergenceError`` is raised.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> final = bsolver.converge(distribution, maxn=1000)
        >>> bool(np.isclose(bsolver._norm(final), 1.0))
        True
        """
        err0 = err1 = 0
        delta = delta0

        for i in range(maxn):
            # If we have already two error estimations we use Richardson
            # extrapolation to obtain a new delta and speed up
            # convergence.
            if 0 < err1 < err0:
                # Log extrapolation attempting to reduce the error a
                # factor m
                delta = delta * np.log(m) / (np.log(err0) - np.log(err1))

            f1 = self.iterate(f0, delta=delta, **kwargs)
            err0 = err1
            err1 = self._norm(abs(f0 - f1))

            logger.debug(
                "After iteration %3d, err = %g (target: %g)",
                i + 1,
                err1,
                rtol,
            )
            if err1 < rtol:
                logger.info(
                    "Convergence achieved after %d iterations. err = %g",
                    i + 1,
                    err1,
                )
                if full:
                    return f1, i + 1, err1

                return f1
            f0 = f1

        logger.error("Convergence failed")

        raise ConvergenceError()

    def _linsystem(self, F):
        r"""Build the transport and collision matrices for an EEDF.

        Parameters
        ----------
        F : array-like
            Current EEDF.

        Returns
        -------
        tuple
            Transport matrix ``A`` and collision matrix ``Q``.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> transport, collisions = bsolver._linsystem(distribution)
        >>> transport.shape, collisions.nnz
        ((30, 30), 0)
        """
        Q = self._PQ(F)

        nu = np.sum(Q.dot(F))

        # Only interior faces enter the zero-flux boundary
        # discretization.
        # Use finite placeholders at the two unused boundary faces.
        sigma_tilde = self.sigma_m.copy()
        sigma_tilde[1:-1] += nu / np.sqrt(self.benergy[1:-1]) / GAMMA

        # The R (G) term, which we add to A.
        G = 2 * self.denergy32 * nu / 3

        A = self._scharf_gummel(sigma_tilde, G)

        return A, Q

    def _norm(self, f: FloatArray) -> float:
        r"""Calculate the EEDF normalization integral.

        Parameters
        ----------
        f : array-like
            EEDF values at cell centers.

        Returns
        -------
        float
            Integral of ``f * sqrt(E)`` over energy.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> bool(bsolver._norm(bsolver.maxwell(2.0)) > 0)
        True
        """
        return simpson(f * np.sqrt(self.cenergy), x=self.cenergy)

    def _normalized(self, f: FloatArray) -> FloatArray:
        r"""Normalize an EEDF using :meth:`_norm`.

        Parameters
        ----------
        f : array-like
            EEDF values at cell centers.

        Returns
        -------
        numpy.ndarray
            Normalized EEDF.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> f = bsolver._normalized(bsolver.maxwell(2.0))
        >>> bool(np.isclose(bsolver._norm(f), 1.0))
        True
        """
        N = self._norm(f)
        return f / N

    def _scharf_gummel(self, sigma_tilde, G=0):
        r"""Build the Scharfetter-Gummel transport matrix.

        Parameters
        ----------
        sigma_tilde : array-like
            Effective momentum-transfer cross section.
        G : float or array-like, optional
            Growth contribution added to the matrix diagonal.

        Returns
        -------
        scipy.sparse.dia_matrix
            Tridiagonal discretization matrix.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bsolver._scharf_gummel(bsolver.sigma_m).shape
        (30, 30)
        """
        D = np.ones_like(self.DA)
        D[1:-1] = self.DA[1:-1] / sigma_tilde[1:-1] + self.DB[1:-1]

        # Boundary values of z are unused because of the zero-flux
        # conditions.
        z = self.W * np.r_[np.nan, np.diff(self.cenergy), np.nan] / D

        a0 = self.W / (1 - np.exp(-z))
        a1 = self.W / (1 - np.exp(z))

        diags = np.zeros((3, self.n))

        # No flux at the energy = 0 boundary
        diags[0, 0] = a0[1]

        diags[0, 1:] = a0[2:] - a1[1:-1]
        diags[1, :] = a1[:-1]
        diags[2, :] = -a0[1:]

        # zero flux b.c.
        diags[2, -2] = -a0[-2]
        diags[0, -1] = -a1[-2]

        diags[0, :] += G

        A = sparse.dia_matrix((diags, [0, 1, -1]), shape=(self.n, self.n))

        return A

    def _g(self, F0):
        r"""Estimate the logarithmic EEDF slope in each cell.

        Parameters
        ----------
        F0 : array-like
            EEDF values at cell centers.

        Returns
        -------
        numpy.ndarray
            Approximation to ``d(log(F0)) / dE``.

        Raises
        ------
        FloatingPointError
            If ``F0`` contains non-finite values.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bool(np.allclose(bsolver._g(np.ones(30)), 0))
        True
        """
        Fp = np.r_[F0[0], F0, F0[-1]]
        cenergyp = np.r_[self.cenergy[0], self.cenergy, self.cenergy[-1]]
        # Zero-valued tails can produce 0/0 or log(0). Diagnose their
        # impact
        # explicitly instead of emitting repeated low-level NumPy
        # warnings.
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            g = np.log(Fp[2:] / Fp[:-2]) / (cenergyp[2:] - cenergyp[:-2])
        if not np.all(np.isfinite(F0)):
            raise FloatingPointError(
                "EEDF must be finite before constructing collision slopes."
            )
        if np.any(F0 < 0) and "negative" not in self._slope_diagnostics_seen:
            self._slope_diagnostics_seen.add("negative")
            negative_mass = np.sum(
                self.denergy32 * np.maximum(-F0, 0)
            ) / np.sum(self.denergy32 * np.abs(F0))
            log = logger.info if negative_mass < 1e-12 else logger.warning
            log(
                "Negative EEDF values during iteration: relative negative "
                "population %.3g; verify the final EEDF and rates.",
                negative_mass,
            )
        bad = ~np.isfinite(g)
        if np.any(bad):
            occupied = bad & (F0 > 0)
            label = "occupied" if np.any(occupied) else "empty"
            if label not in self._slope_diagnostics_seen:
                self._slope_diagnostics_seen.add(label)
                if label == "empty":
                    logger.info(
                        "Nonfinite slopes only in nonpositive-EEDF cells: "
                        "legacy zero-integral fallback; check final EEDF "
                        "nonnegativity."
                    )
                else:
                    mass = np.sum(
                        self.denergy32[occupied] * F0[occupied]
                    ) / np.sum(self.denergy32 * F0)
                    logger.warning(
                        "EEDF tail slopes undefined in %d occupied cells "
                        "(energy %.6g..%.6g eV; discrete population fraction "
                        "%.3g). Legacy zero-integral fallback retained; "
                        "check grid/cutoff convergence before ignoring.",
                        np.count_nonzero(occupied),
                        self.cenergy[occupied].min(),
                        self.cenergy[occupied].max(),
                        mass,
                    )

        return g

    def _PQ(self, F0, reactions=None):
        r"""Build the sparse collision transition matrix.

        Parameters
        ----------
        F0 : array-like
            Current EEDF used to estimate collision slopes.
        reactions : iterable of tuple, optional
            ``(target, process)`` pairs. Defaults to active inelastic
            processes.

        Returns
        -------
        scipy.sparse.coo_matrix
            Collision matrix mapping the EEDF to collision source terms.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bsolver._PQ(distribution).nnz
        0
        """
        g = self._g(F0)
        if reactions is None:
            reactions = self.iter_inelastic()

        data = []
        rows = []
        cols = []
        for t, k in reactions:
            r = t.density * GAMMA * k.scatterings(g, self.cenergy)
            in_factor = k.in_factor
            if in_factor is None:
                raise ValueError("Process %s has no in_factor defined." % k)

            data.extend([in_factor * r, -r])
            rows.extend([k.i, k.j])
            cols.extend([k.j, k.j])

        if not data:
            return sparse.coo_matrix((self.n, self.n))

        data, rows, cols = (np.hstack(x) for x in (data, rows, cols))
        PQ = sparse.coo_matrix(
            (data, (rows, cols)),
            shape=(self.n, self.n),
        )

        return PQ

    def rate(
        self, F0: FloatArray, k: Process | str, weighted: bool = False
    ) -> float:
        r"""Calculate a collision-process rate from an EEDF.

        Parameters
        ----------
        F0 : array-like
            EEDF values at cell centers.
        k : process.Process or str
            Process instance or signature accepted by :meth:`search`.
        weighted : bool, optional
            If ``True``, multiply the rate by the target density.

        Returns
        -------
        float
            Rate of the process for the given EEDF.

        See Also
        --------
        search

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bool(bsolver.rate(distribution, collision) > 0)
        True
        """
        g = self._g(F0)

        if isinstance(k, str):
            k = self.search(k)

        k.set_grid_cache(self.grid)
        r = k.scatterings(g, self.cenergy)
        rate = GAMMA * np.dot(F0[k.j], r)
        if weighted:
            target = k.target
            if target is None:
                raise ValueError(
                    "Process must be attached to a target before "
                    "calculating a weighted rate."
                )
            rate *= target.density

        return rate

    def mobility(self, F0: FloatArray) -> float:
        r"""Calculate the reduced electron mobility from an EEDF.

        Parameter
        ----------
        F0 : array-like
            EEDF values at cell centers.

        Returns
        -------
        float
            Reduced mobility (mobility times gas density) in SI units.

        See Also
        --------
        diffusion

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bool(np.isfinite(bsolver.mobility(distribution)))
        True
        """
        DF0 = np.r_[0.0, np.diff(F0) / np.diff(self.cenergy), 0.0]
        Q = self._PQ(F0, reactions=self.iter_growth())

        nu = np.sum(Q.dot(F0)) / GAMMA
        # The lower endpoint is explicitly zero in the mobility
        # quadrature.
        sigma_tilde = self.sigma_m[1:] + nu / np.sqrt(self.benergy[1:])
        y = np.zeros_like(self.benergy)
        y[1:] = DF0[1:] * self.benergy[1:] / sigma_tilde

        return -(GAMMA / 3) * simpson(y, x=self.benergy)

    def diffusion(self, F0: FloatArray) -> float:
        r"""Calculate the reduced electron diffusion coefficient.

        Parameters
        ----------
        F0 : array-like
            EEDF values at cell centers.

        Returns
        -------
        float
            Reduced diffusion coefficient in SI units.

        See Also
        --------
        mobility

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> collision = bsolver.add_process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     mass_ratio=1e-5,
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> bsolver.set_density("Ar", 1.0)
        >>> bsolver.EN = 100 * TOWNSEND
        >>> bsolver.kT = 0.026
        >>> bsolver.init()
        >>> distribution = bsolver.maxwell(2.0)
        >>> bool(bsolver.diffusion(distribution) > 0)
        True
        """
        Q = self._PQ(F0, reactions=self.iter_growth())
        nu = np.sum(Q.dot(F0)) / GAMMA
        sigma_m = np.zeros_like(self.cenergy)

        for target, process in self.iter_momentum():
            s = target.density * process.interp(self.cenergy)
            sigma_m += s

        sigma_tilde = sigma_m + nu / np.sqrt(self.cenergy)
        y = F0 * self.cenergy / sigma_tilde
        return (GAMMA / 3) * simpson(y, x=self.cenergy)

    def mean_energy(self, F0: FloatArray) -> float:
        r"""Calculate the mean electron energy.

        Parameters
        ----------
        F0 : array-like
            EEDF values at cell centers.

        Returns
        -------
        float
            Mean electron energy in eV.

        Examples
        --------
        >>> from bolos.grid import LinearGrid
        >>> bsolver = BoltzmannSolver(LinearGrid(0, 10, 30))
        >>> bool(bsolver.mean_energy(bsolver.maxwell(2.0)) > 0)
        True
        """
        de52 = np.diff(self.benergy**2.5)
        return np.sum(0.4 * F0 * de52)
