r"""Parse BOLSIG+-compatible cross-section files independently of BOLOS.

Most users only need :func:`parse`. Processes are returned as
dictionaries
so that this module can also be reused outside BOLOS.

Examples
--------
>>> from bolos.parser import parse
>>> from io import StringIO
>>> parse(StringIO("Unrecognized header\n"))
[]
"""

import logging
import re
from collections.abc import Iterator
from typing import TypedDict

import numpy as np

logger = logging.getLogger(__name__)


class ProcessDict(TypedDict, total=False):
    r"""Describe parsed fields accepted by the collision constructors.

    Notes
    -----
    This TypedDict describes ordinary dictionaries at runtime. Fields
    vary
    with the process kind; parsing always supplies kind, target,
    comment,
    and data. The data table remains two-dimensional for a single row.

    Examples
    --------
    >>> definition = ProcessDict(kind="ELASTIC", target="Ar")
    >>> definition["target"]
    'Ar'
    """

    kind: str
    target: str
    product: str
    comment: str
    data: list[list[float]]
    mass_ratio: float
    threshold: float
    weight_ratio: float


def parse(fp: Iterator[str]) -> list[ProcessDict]:
    r"""Parse processes from a LXCat cross-section file.

    Parameters
    ----------
    fp : iterator of str
        Text stream positioned before the processes to read. Each
        process starts with a keyword recognized by `KEYWORDS`.

    Returns
    -------
    processes : list of dict
        Process dictionaries in file order. Each contains ``kind``,
        ``target``, ``comment``, and ``data``, plus fields supplied by
        the reader for that process type.

    Raises
    ------
    StopIteration
        A process block ends before its target or required argument.
    ValueError
        A numeric argument, reaction, or data table cannot be parsed.
    IndexError
        A required numeric argument is missing from its line.

    Notes
    -----
    Unrecognized lines outside process blocks are ignored. The stream is
    consumed but is not closed. To construct a BOLOS process from a
    returned dictionary ``d``, use ``bolos.process.Process(**d)`` after
    importing ``bolos.process``.

    Examples
    --------
    >>> from io import StringIO
    >>> text = "ATTACHMENT\nO2 -> O-\nExample\n-----\n"
    >>> text += "0 0\n1 2e-20\n-----\n"
    >>> processes = parse(StringIO(text))
    >>> processes[0]["kind"], processes[0]["target"]
    ('ATTACHMENT', 'O2')
    >>> processes[0]["data"]
    [[0.0, 0.0], [1.0, 2e-20]]
    """
    processes: list[ProcessDict] = []
    for line in fp:
        key = line.strip()
        reader = KEYWORDS.get(key)
        if reader is None:
            continue
        logger.debug("New process of type '%s'", key)
        process = reader(fp)
        process["kind"] = key
        processes.append(process)

    logger.info("Parsing complete. %d processes read.", len(processes))

    return processes


# BOLSIG+ separators begin with at least five dashes.
RE_SEP = re.compile("-----+")


def _read_until_sep(fp: Iterator[str]) -> list[str]:
    r"""Read stripped lines up to a separator or the end of the stream.

    Parameters
    ----------
    fp : iterator of str
        Text stream to consume from its current position.

    Returns
    -------
    lines : list of str
        Lines with leading and trailing whitespace removed. The
        separator is consumed and excluded from the result. Blank lines
        are retained.

    Notes
    -----
    A separator is any line whose stripped text starts with at least
    five dashes. Reaching the end of the stream raises no exception.

    Examples
    --------
    >>> from io import StringIO
    >>> stream = StringIO("  comment  \n\n-----\nremaining\n")
    >>> _read_until_sep(stream)
    ['comment', '']
    >>> next(stream).strip()
    'remaining'
    """
    lines = []
    for line in fp:
        if RE_SEP.match(line.strip()):
            break
        lines.append(line.strip())

    return lines


def _read_block(
    fp: Iterator[str], has_arg: bool = True
) -> tuple[str, str | None, str, list[list[float]]]:
    r"""Read a process block after its keyword line.

    Parameters
    ----------
    fp : iterator of str
        Text stream positioned at the target line. The block contains a
        target, an optional argument line, comments terminated by a
        separator, and a numeric table terminated by another separator
        or the end of the stream.
    has_arg : bool, optional
        Whether to read an argument after the target. Default is True.

    Returns
    -------
    target : str
        Stripped target or reaction line.
    arg : str or None
        Stripped argument line, or None when `has_arg` is False.
    comment : str
        Stripped comment lines joined by newline characters.
    data : list of list of float
        Numeric table retaining two dimensions, including for a single
        row.

    Raises
    ------
    StopIteration
        The target or required argument line is absent.
    ValueError
        The data table cannot be parsed by `numpy.loadtxt`.

    Examples
    --------
    >>> from io import StringIO
    >>> text = "Ar\n1e-5\nExample\n-----\n0 0\n1 2e-20\n-----\n"
    >>> target, arg, comment, data = _read_block(StringIO(text))
    >>> target, arg, comment
    ('Ar', '1e-5', 'Example')
    >>> data
    [[0.0, 0.0], [1.0, 2e-20]]
    """
    target = next(fp).strip()
    if has_arg:
        arg = next(fp).strip()
    else:
        arg = None

    comment = "\n".join(_read_until_sep(fp))

    logger.debug("Read process '%s'", target)
    data = np.loadtxt(_read_until_sep(fp), ndmin=2).tolist()

    return target, arg, comment, data


# Specialized readers return dictionaries containing process attributes.
def _read_momentum(fp: Iterator[str]) -> ProcessDict:
    r"""Read a MOMENTUM, ELASTIC, or EFFECTIVE process block.

    Parameters
    ----------
    fp : iterator of str
        Text stream positioned at the target line after the keyword. The
        next line supplies the mass ratio as its first value, followed
        by comments and the numeric table.

    Returns
    -------
    process : dict
        Dictionary containing ``target`` (str), ``mass_ratio`` (float),
        ``comment`` (str), and ``data`` (as described in `_read_block`).
        The process keyword is not included.

    Raises
    ------
    StopIteration
        The target or argument line is absent.
    ValueError
        The mass ratio or numeric table cannot be parsed.
    IndexError
        The argument line contains no mass ratio.

    Examples
    --------
    >>> from io import StringIO
    >>> text = "Ar\n1e-5\nExample\n-----\n0 0\n1 2e-20\n-----\n"
    >>> process = _read_momentum(StringIO(text))
    >>> process["target"], process["mass_ratio"]
    ('Ar', 1e-05)
    """
    target, arg, comment, data = _read_block(fp, has_arg=True)
    mass_ratio = float(arg.split()[0])
    d = dict(
        target=target,
        mass_ratio=mass_ratio,
        comment=comment,
        data=data,
    )

    return d


RE_ARROW = re.compile("<?->")


def _read_excitation(fp: Iterator[str]) -> ProcessDict:
    r"""Read an EXCITATION or IONIZATION process block.

    Parameters
    ----------
    fp : iterator of str
        Text stream positioned at a reaction line containing ``->`` or
        ``<->``. The argument line gives the threshold energy and, for a
        reversible reaction, the statistical weight ratio. Comments and
        the numeric table follow.

    Returns
    -------
    process : dict
        Dictionary containing ``target`` and ``product`` (str),
        ``comment`` (str), ``data`` (as described in `_read_block`), and
        ``threshold`` (float). A reversible reaction also includes
        ``weight_ratio`` (float). The process keyword is not included.

    Raises
    ------
    StopIteration
        The reaction or argument line is absent.
    ValueError
        The reaction does not split into exactly two sides, or a numeric
        argument or table cannot be parsed.
    IndexError
        A required numeric argument is missing.

    Notes
    -----
    The reversible arrow ``<->`` must be surrounded by whitespace for
    the weight ratio to be read.

    Examples
    --------
    >>> from io import StringIO
    >>> text = "Ar <-> Ar*\n11.5 3\nExample\n-----\n"
    >>> text += "12 0\n13 2e-20\n-----\n"
    >>> process = _read_excitation(StringIO(text))
    >>> process["target"], process["product"]
    ('Ar', 'Ar*')
    >>> process["threshold"], process["weight_ratio"]
    (11.5, 3.0)
    """
    target, arg, comment, data = _read_block(fp, has_arg=True)
    lhs, rhs = [s.strip() for s in RE_ARROW.split(target)]

    d = dict(
        target=lhs,
        product=rhs,
        comment=comment,
        data=data,
    )

    if "<->" in target.split():
        threshold, weight_ratio = float(arg.split()[0]), float(arg.split()[1])
        d["weight_ratio"] = weight_ratio
    else:
        threshold = float(arg.split()[0])

    d["threshold"] = threshold
    return d


def _read_attachment(fp: Iterator[str]) -> ProcessDict:
    r"""Read an ATTACHMENT process block without an argument line.

    Parameters
    ----------
    fp : iterator of str
        Text stream positioned at the target or reaction line after the
        keyword. Comments and the numeric table follow immediately.

    Returns
    -------
    process : dict
        Dictionary containing ``target`` (str), ``comment`` (str),
        ``data`` (as described in `_read_block`), and ``threshold``
        (always 0.0). If the reaction splits into exactly two sides at
        ``->`` or ``<->``, ``product`` (str) is also included. Otherwise
        the whole target line is retained. The keyword is not included.

    Raises
    ------
    StopIteration
        The target line is absent.
    ValueError
        The numeric table cannot be parsed.

    Examples
    --------
    >>> from io import StringIO
    >>> text = "O2 -> O-\nExample\n-----\n0 0\n1 2e-20\n-----\n"
    >>> process = _read_attachment(StringIO(text))
    >>> (
    ...     process["target"],
    ...     process["product"],
    ...     process["threshold"],
    ... )
    ('O2', 'O-', 0.0)
    >>> text = "O2\nExample\n-----\n0 0\n1 2e-20\n-----\n"
    >>> "product" in _read_attachment(StringIO(text))
    False
    """
    target, arg, comment, data = _read_block(fp, has_arg=False)

    d = dict(
        comment=comment,
        data=data,
        threshold=0.0,
    )
    lr = [s.strip() for s in RE_ARROW.split(target)]

    if len(lr) == 2:
        d["target"] = lr[0]
        d["product"] = lr[1]
    else:
        d["target"] = target

    return d


KEYWORDS = {
    "MOMENTUM": _read_momentum,
    "ELASTIC": _read_momentum,
    "EFFECTIVE": _read_momentum,
    "EXCITATION": _read_excitation,
    "IONIZATION": _read_excitation,
    "ATTACHMENT": _read_attachment,
}
