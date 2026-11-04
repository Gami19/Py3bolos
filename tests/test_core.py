"""Regression checks for collision loading, empty reactions, and rate reduction."""

from io import StringIO
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from scipy import sparse

from bolos import grid, parser, solver
from bolos.process import Process, padinterp
from bolos.target import Target


class CoreTests(unittest.TestCase):
    def make_solver(self):
        result = solver.BoltzmannSolver(grid.LinearGrid(0, 10, 30))
        result.add_process(target="Ar", kind="ELASTIC", mass_ratio=1e-5,
                           data=[[0, 1e-20], [10, 1e-20]])
        result.set_density("Ar", 1)
        result.EN = 100 * solver.TOWNSEND
        result.kT = 0.026
        result.init()
        return result

    def test_empty_reactions_and_no_growth_transport(self):
        result = self.make_solver()
        f = result.maxwell(2)
        for reactions in (None, [], iter(())):
            matrix = result._PQ(f, reactions=reactions)
            self.assertEqual(matrix.shape, (30, 30))
            self.assertEqual(matrix.nnz, 0)
        self.assertTrue(np.isfinite(result.mobility(f)))
        self.assertGreater(result.diffusion(f), 0)
        final = result.converge(f, maxn=1000)
        self.assertTrue(np.all(np.isfinite(final)))
        assert_allclose(result._norm(final), 1)

    def test_effective_aliases_remain_synchronized(self):
        target = Target("Ar")
        target.add_process(Process(target="Ar", kind="EFFECTIVE", mass_ratio=1e-5,
                                   data=[[0, 1e-20], [10, 1e-20]]))
        effective = target.effective
        target.ensure_elastic()
        self.assertIs(target.effective, effective)
        for key in ("EFFECTIVE", "MOMENTUM"):
            self.assertIs(target.kind[key], effective)
            self.assertEqual(target.kind[key], [])
        target.ensure_elastic()  # Repeated loading must not report a conflict.
        self.assertEqual(len(target.elastic), 1)

    def test_public_exports(self):
        namespace = {}
        exec("from bolos import *", namespace)
        self.assertIs(namespace["parser"], parser)

    def test_single_row_table(self):
        text = "ELASTIC\nAr\n1e-5\nExample\n-----\n0 1e-20\n-----\n"
        definition = parser.parse(StringIO(text))[0]
        self.assertEqual(definition["data"], [[0, 1e-20]])
        process = Process(**definition)
        assert_allclose(process.interp([0, 1, 10]), 1e-20)

    def test_invalid_cross_sections(self):
        for data in ([], [0, 1], [[0, 1, 2]], [[0, np.nan]], [[0, np.inf]],
                     [[-1, 1]], [[0, -1]], [[0, 1], [0, 2]], [[1, 1], [0, 2]]):
            with self.subTest(data=data), self.assertRaises(ValueError):
                Process(target="Ar", kind="ELASTIC", data=data)

    def test_parser_reader_errors_are_not_silenced(self):
        with patch.dict(parser.KEYWORDS, {"BROKEN": lambda fp: {}["missing"]}):
            with self.assertRaises(KeyError):
                parser.parse(StringIO("Ignored header\nBROKEN\n"))

    def test_rate_matches_sparse_reduction_with_repeated_indices(self):
        result = self.make_solver()
        process = result.add_process(
            target="Ar", kind="EXCITATION", product="Ar*", threshold=1,
            data=[[0, 0], [1, 0], [1.1, 1e-20], [1.2, 2e-20], [10, 1e-20]],
        )
        result.init()
        f = result.maxwell(2)
        process.set_grid_cache(result.grid)
        self.assertGreater(len(process.j), len(np.unique(process.j)))
        r = process.scatterings(result._g(f), result.cenergy)
        reference = sparse.coo_matrix(
            (solver.GAMMA * r, (process.j, np.zeros(len(r), dtype=int))),
            shape=(result.n, 1),
        ).toarray().ravel()
        expected = f.dot(reference)
        assert_allclose(result.rate(f, process), expected, rtol=1e-14, atol=0)
        result.set_density("Ar", 0.25)
        assert_allclose(result.rate(f, "Ar -> Ar*", weighted=True), expected * 0.25)

    def test_padded_interpolation_keeps_range_behavior(self):
        interpolate = padinterp(np.array([[1., 2.], [3., 4.]]))
        assert_allclose(interpolate([0, 1, 2, 3, 1e8]), [2, 2, 3, 4, 4])
        for energy in (-1, 1e8 + 1):
            with self.assertRaises(ValueError):
                interpolate(energy)


if __name__ == "__main__":
    unittest.main()
