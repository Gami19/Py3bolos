"""Run package documentation examples with the regular regression suite."""

import doctest
import importlib


def load_tests(loader, tests, pattern):
    """Include executable examples from every BOLOS module.

    Parameters
    ----------
    loader : unittest.TestLoader
        Loader supplied by unittest discovery.
    tests : unittest.TestSuite
        Suite to extend with documentation examples.
    pattern : str or None
        Discovery filename pattern, unused by this hook.

    Returns
    -------
    unittest.TestSuite
        Suite including independent module, class, and method examples.
    """
    for name in (
        "bolos",
        "bolos.grid",
        "bolos.parser",
        "bolos.process",
        "bolos.solver",
        "bolos.target",
    ):
        module = importlib.import_module(name)
        tests.addTests(doctest.DocTestSuite(module))
    # DocTestFinder discovers a property's getter, but not its setter.
    module = importlib.import_module("bolos.solver")
    setter = module.BoltzmannSolver.grid.fset
    example = doctest.DocTestParser().get_doctest(
        setter.__doc__, vars(module), "BoltzmannSolver.grid.setter", None, 0
    )
    tests.addTest(doctest.DocTestCase(example))
    return tests
