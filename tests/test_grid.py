"""Grid geometry and regression tests. Run with unittest discover -s tests."""

import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from bolos.grid import (
    AutomaticGrid,
    GeometricGrid,
    Grid,
    LinearGrid,
    LogGrid,
    MappedGrid,
    QuadraticGrid,
    mkgrid,
)


class GridTests(unittest.TestCase):
    def test_mappings_and_shared_geometry(self):
        for cls in (LinearGrid, QuadraticGrid, GeometricGrid, LogGrid):
            with self.subTest(grid=cls.__name__):
                grid = cls(0, 60, 100)
                self.assertEqual(grid.b.shape, (101,))
                self.assertEqual(grid.c.shape, (100,))
                self.assertTrue(np.all(grid.d > 0))
                assert_allclose(grid.finv(grid.f(grid.c)), grid.c, atol=1e-12)
                assert_allclose(np.diff(grid.f(grid.b)), grid.df, atol=1e-12)
                assert_allclose(np.sum(grid.d), 60)
                assert_allclose(np.sum(grid.d32), 60**1.5)

    def test_missing_mapping_cannot_be_instantiated(self):
        class MissingInverse(MappedGrid):
            def f(self, x):
                return x

        with self.assertRaises(TypeError):
            MissingInverse(0, 1, 10)

    def test_cell_lookup_for_all_grid_types(self):
        source = LinearGrid(0, 60, 20)
        automatic = AutomaticGrid(source, np.exp(-source.c / 2))
        for grid in (source, QuadraticGrid(0, 60, 20), automatic):
            with self.subTest(grid=type(grid).__name__):
                self.assertEqual(grid.cell(grid.x0), 0)
                self.assertEqual(grid.cell(grid.x1), grid.n - 1)
                for i, center in enumerate(grid.c):
                    self.assertEqual(grid.cell(float(center)), i)
                for i, boundary in enumerate(grid.b[1:-1], start=1):
                    self.assertEqual(grid.cell(float(boundary)), i)
                for x in (-1, 61, np.nan, np.inf):
                    with self.assertRaises(ValueError):
                        grid.cell(x)

    def test_interpolation_refreshes_distribution_and_source(self):
        destination = LinearGrid(0, 6, 6)
        source = LinearGrid(1, 5, 4)
        assert_allclose(destination.interpolate(np.ones(4), source), [0, 1, 1, 1, 1, 0])
        assert_allclose(destination.interpolate(np.full(4, 2), source), [0, 2, 2, 2, 2, 0])
        wider_source = LinearGrid(0, 6, 3)
        assert_allclose(destination.interpolate(np.full(3, 3), wider_source), 3)

    def test_interpolation_of_linear_function(self):
        source = LinearGrid(0, 10, 10)
        destination = QuadraticGrid(1, 9, 30)
        assert_allclose(destination.interpolate(2 * source.c + 1, source), 2 * destination.c + 1)

    def test_automatic_grid_population_and_reuse(self):
        source = LinearGrid(0, 10, 10)
        f0 = np.exp(-source.c)
        automatic = AutomaticGrid(source, f0)
        cumulative = np.r_[0, np.cumsum(source.d32 * f0)]
        cumulative /= cumulative[-1]
        assert_allclose(np.interp(automatic.b, source.b, cumulative), np.linspace(0, 1, 11))
        self.assertEqual(automatic.d32.shape, (10,))
        # A second redistribution needs d32, missing in the old implementation.
        repeated = AutomaticGrid(automatic, np.ones(10))
        self.assertTrue(np.all(repeated.d > 0))

    def test_automatic_grid_with_empty_cells(self):
        source = LinearGrid(0, 6, 6)
        automatic = AutomaticGrid(source, [0, 1, 0, 1, 0, 0])
        self.assertEqual(automatic.x0, source.x0)
        self.assertEqual(automatic.x1, source.x1)
        self.assertTrue(np.all(automatic.d > 0))

    def test_invalid_inputs(self):
        for boundaries in ([0], [0, 0, 1], [-1, 0], [0, np.nan], [[0, 1]]):
            with self.assertRaises(ValueError):
                Grid(boundaries)
        source = LinearGrid(0, 3, 3)
        for f0 in ([0, 0, 0], [-1, 1, 1], [1, np.nan, 1], [1, 2]):
            with self.assertRaises(ValueError):
                AutomaticGrid(source, f0)
        with self.assertRaises(ValueError):
            source.interpolate([1, 2], source)
        for args in ((0, 1, 0), (1, 0, 10), (-1, 1, 10)):
            with self.assertRaises(ValueError):
                LinearGrid(*args)
        with self.assertRaises(ValueError):
            GeometricGrid(0, 10, 10, r=1)
        with self.assertRaises(ValueError):
            LogGrid(0, 10, 10, s=0)

    def test_factory_aliases_and_parameters(self):
        for names, cls in ((('linear', 'lin'), LinearGrid),
                           (('quadratic', 'quad'), QuadraticGrid),
                           (('geometric', 'geo'), GeometricGrid),
                           (('logarithmic', 'log'), LogGrid)):
            for name in names:
                self.assertIsInstance(mkgrid(name, 0, 10, 10), cls)
        assert_array_equal(mkgrid('geo', 0, 10, 10, 1.2).b, GeometricGrid(0, 10, 10, r=1.2).b)
        assert_array_equal(mkgrid('log', 0, 10, 10, s=0.1).b, LogGrid(0, 10, 10, s=0.1).b)


if __name__ == '__main__':
    unittest.main()
