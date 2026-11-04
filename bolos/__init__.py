# Usually these are the only modules that are useful to the user.
r"""Solve electron Boltzmann equations using tabulated cross sections.

Examples
--------
>>> from bolos import grid, parser, process, solver
>>> grid.LinearGrid(0, 10, 5).n
5
"""

__all__ = ["solver", "grid", "parser", "process"]
