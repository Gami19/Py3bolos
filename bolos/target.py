r"""Group collision processes by target species and reaction kind.

Examples
--------
>>> from bolos.target import Target
>>> Target("Ar").name
'Ar'
"""

from __future__ import annotations

import logging
from collections import defaultdict

import numpy as np

from .process import Process

logger = logging.getLogger(__name__)


class Target:
    r"""Store collision processes and composition for one species.

    Parameters
    ----------
    name : str
        Species name used to identify the collision target.

    Attributes
    ----------
    name : str
        Target species name.
    density : float
        Molar fraction, initially zero.
    mass_ratio : float or None
        Electron-to-target mass ratio supplied by a registered process.
    elastic, effective : list of Process
        Elastic and effective collision processes.
    attachment, ionization, excitation : list of Process
        Inelastic collision processes.
    weighted_elastic : list of Process
        Compatibility storage for weighted elastic processes.
    kind : dict of str to list of Process
        Kind lookup. MOMENTUM and EFFECTIVE share the same list.
    by_product : defaultdict of list
        Processes indexed by product name.
    clipped_elastic_range : tuple of float or None
        Energy range in eV affected by clipping a derived elastic table.

    Examples
    --------
    >>> from bolos.target import Target
    >>> target = Target("Ar")
    >>> target.density = 1.0
    >>> target.name, target.density
    ('Ar', 1.0)
    """

    def __init__(self, name: str) -> None:
        r"""Initialize an empty target with zero molar fraction.

        Parameters
        ----------
        name : str
            Species name.

        Examples
        --------
        >>> Target("Ar").density
        0.0
        """
        self.name = name
        self.mass_ratio = None
        self.density = 0.0
        self.clipped_elastic_range = None

        # Lists of all processes pertaining this target
        self.elastic = []
        self.effective = []
        self.attachment = []
        self.ionization = []
        self.excitation = []
        self.weighted_elastic = []

        self.kind = {
            "ELASTIC": self.elastic,
            "EFFECTIVE": self.effective,
            "MOMENTUM": self.effective,
            "ATTACHMENT": self.attachment,
            "IONIZATION": self.ionization,
            "EXCITATION": self.excitation,
            "WEIGHTED_ELASTIC": self.weighted_elastic,
        }

        self.by_product = defaultdict(list)

        logger.debug("Target %s created.", self)

    def add_process(self, process: Process) -> None:
        r"""Register a process in its kind and product lookup tables.

        Parameters
        ----------
        process : Process
            Collision process to associate with this target.

        Raises
        ------
        KeyError
            If the process kind is not supported by this target.
        ValueError
            If the process supplies a conflicting mass ratio.

        Notes
        -----
        The process's ``target`` attribute is set to this instance.

        Examples
        --------
        >>> from bolos.process import Process
        >>> target = Target("Ar")
        >>> collision = Process(
        ...     target="Ar",
        ...     kind="ELASTIC",
        ...     data=[[0, 1e-20], [10, 1e-20]],
        ... )
        >>> target.add_process(collision)
        >>> target.elastic[0] is collision
        True
        """
        kind = self.kind[process.kind]
        kind.append(process)

        if process.mass_ratio is not None:
            logger.debug("Mass ratio (=%g) for %s", process.mass_ratio, self)

            if (
                self.mass_ratio is not None
                and self.mass_ratio != process.mass_ratio
            ):
                raise ValueError(
                    "More than one mass ratio for target '%s'" % self.name
                )

            self.mass_ratio = process.mass_ratio

        process.target = self

        self.by_product[process.product].append(process)

        logger.debug("Process %s added to target %s", process, self)

    def ensure_elastic(self) -> None:
        r"""Derive an elastic table from an effective cross section.

        Raises
        ------
        ValueError
            If both elastic and effective tables exist, or multiple
            effective
            tables are registered.

        Notes
        -----
        Inelastic cross sections are subtracted on the effective table's
        energy
        nodes. Negative results are clipped to zero and logged. If no
        elastic
        or effective process exists, a warning is logged and no table is
        added.

        Examples
        --------
        >>> from bolos.process import Process
        >>> target = Target("Ar")
        >>> target.add_process(
        ...     Process(
        ...         target="Ar",
        ...         kind="EFFECTIVE",
        ...         mass_ratio=1e-5,
        ...         data=[[0, 1e-20], [10, 1e-20]],
        ...     )
        ... )
        >>> target.ensure_elastic()
        >>> len(target.elastic), len(target.effective)
        (1, 0)
        """
        if self.elastic and self.effective:
            raise ValueError(
                "In target '%s': EFFECTIVE/MOMENTUM and ELASTIC"
                "cross-sections are incompatible." % self
            )

        if self.elastic:
            return

        if len(self.effective) > 1:
            raise ValueError(
                "In target '%s': Can't handle more that 1 "
                "EFFECTIVE/MOMENTUM for a given target" % self
            )

        if not self.effective:
            logger.warning(
                "Target %s has no ELASTIC or EFFECTIVE cross sections", self
            )
            return

        newdata = self.effective[0].data.copy()
        for p in self.inelastic:
            newdata[:, 1] -= p.interp(newdata[:, 0])

        if np.amin(newdata[:, 1]) < 0:
            bad = np.flatnonzero(newdata[:, 1] < 0)
            self.clipped_elastic_range = (
                newdata[max(0, bad[0] - 1), 0],
                newdata[min(len(newdata) - 1, bad[-1] + 1), 0],
            )
            logger.warning(
                (
                    "Target %s: EFFECTIVE minus INELASTIC is negative at "
                    "%d table nodes (%.6g..%.6g eV); clipped to zero. "
                    "Interpolated ELASTIC may be affected within "
                    "%.6g..%.6g eV. This can reflect constant extrapolation "
                    "beyond individual tables, not necessarily invalid "
                    "source data."
                ),
                self.name,
                len(bad),
                newdata[bad[0], 0],
                newdata[bad[-1], 0],
                *self.clipped_elastic_range,
            )
            newdata[:, 1] = np.maximum(newdata[:, 1], 0)

        newelastic = Process(
            target=self.name,
            kind="ELASTIC",
            data=newdata,
            mass_ratio=self.effective[0].mass_ratio,
            comment="Calculated from EFFECTIVE cross sections",
        )

        logger.debug("EFFECTIVE -> ELASTIC for target %s", self)
        self.add_process(newelastic)

        # Remove the EFFECTIVE processes.
        self.effective.clear()

    @property
    def inelastic(self) -> list[Process]:
        r"""Return attachment, ionization, and excitation processes.

        Returns
        -------
        list of Process
            A new list of the selected processes, in group order.

        Examples
        --------
        >>> target = Target("Ar")
        >>> target.inelastic
        []
        """
        return self.attachment + self.ionization + self.excitation

    @property
    def everything(self) -> list[Process]:
        r"""Return elastic and all inelastic processes.

        Returns
        -------
        list of Process
            A new list of the selected processes, in group order.

        Examples
        --------
        >>> target = Target("Ar")
        >>> target.everything
        []
        """
        return (
            self.elastic + self.attachment + self.ionization + self.excitation
        )

    def __repr__(self) -> str:
        r"""Return a representation containing the target name.

        Returns
        -------
        str
            Constructor-like representation.

        Examples
        --------
        >>> print(repr(Target("Ar")))
        Target('Ar')
        """
        return "Target(%s)" % repr(self.name)

    def __str__(self) -> str:
        r"""Return the target species name.

        Returns
        -------
        str
            Text representation of the target.

        Examples
        --------
        >>> str(Target("Ar"))
        'Ar'
        """
        return self.name
