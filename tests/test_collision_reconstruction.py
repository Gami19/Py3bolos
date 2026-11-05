"""Check local exponential reconstruction against independent integrals."""

import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad

from bolos import grid, solver
from bolos.process import Process, int_linexp0


class CollisionReconstructionTests(unittest.TestCase):
    """Keep helper and distribution-slope conventions consistent."""

    def test_helper_keeps_its_reverse_coordinate_contract(self):
        """Verify the helper independently for both signs and zero slope."""
        for slope in (-0.5, -1e-12, 0.0, 1e-12, 0.5):
            with self.subTest(slope=slope):
                expected = quad(
                    lambda energy: (
                        energy
                        * (0.5 + energy)
                        * np.exp(slope * (1.5 - energy))
                    ),
                    1.0,
                    2.0,
                    epsabs=1e-12,
                    epsrel=1e-12,
                )[0]
                actual = int_linexp0(1, 2, 1.5, 2.5, slope, 1.5)
                assert_allclose(actual, expected, rtol=1e-12)

    def test_logarithmic_slope_retains_derivative_definition(self):
        """Verify the derivative sign on uniform and nonuniform grids."""
        for cls in (grid.LinearGrid, grid.QuadraticGrid):
            s = solver.BoltzmannSolver(cls(0, 8, 16))
            for slope in (-0.5, 0.0, 0.5):
                with self.subTest(grid=cls.__name__, slope=slope):
                    f = np.exp(slope * s.cenergy)
                    assert_allclose(s._g(f), slope, atol=1e-13)

    def test_scattering_interval_uses_forward_energy_reconstruction(self):
        """Detect a sign reversal even when no solver iteration is run."""
        energy_grid = grid.LinearGrid(0, 4, 4)
        collision = Process(
            target="X", kind="EXCITATION", data=[[0, 1], [4, 1]]
        )
        collision.set_grid_cache(energy_grid)
        index = np.flatnonzero(
            (collision.j == 1) & (collision.eps[:, 0] == 1)
        )[0]
        original = Process.scatterings
        for slope in (-0.5, 0.5):
            with self.subTest(slope=slope):
                slopes = np.full(4, slope)
                expected = quad(
                    lambda energy: energy * np.exp(slope * (energy - 1.5)),
                    1,
                    2,
                )[0]
                actual = collision.scatterings(slopes, energy_grid.c)[index]
                assert_allclose(actual, expected, rtol=1e-12)

                # The old handoff must fail this independent oracle.
                def reversed_scatterings(process, g, eps):
                    return original(process, -g, eps)

                with patch.object(
                    Process, "scatterings", reversed_scatterings
                ):
                    reversed_value = collision.scatterings(
                        slopes, energy_grid.c
                    )[index]
                self.assertGreater(abs(reversed_value / expected - 1), 0.01)

    def test_rate_and_collision_matrix_match_manufactured_distribution(self):
        """Check source/destination contributions for three reaction kinds."""
        for cls in (grid.LinearGrid, grid.QuadraticGrid):
            for kind, threshold in (
                ("EXCITATION", 0.7),
                ("IONIZATION", 1.1),
                ("ATTACHMENT", 0.0),
            ):
                s = solver.BoltzmannSolver(cls(0, 8, 16))
                collision = s.add_process(
                    target="X",
                    kind=kind,
                    threshold=threshold,
                    data=[[0, 1], [1.3, 1.4], [3.2, 0.8], [8, 1.8]],
                )
                s.set_density("X", 0.7)
                collision.set_grid_cache(s.grid)
                for slope in (-0.5, 0.0, 0.5):
                    with self.subTest(
                        grid=cls.__name__, kind=kind, slope=slope
                    ):
                        f = np.exp(slope * s.cenergy)
                        expected = np.zeros(s.n)
                        exact_rate = 0.0
                        for interval, (a, b) in enumerate(collision.eps):
                            u0, u1 = collision.sigma[interval]
                            value = (
                                solver.GAMMA
                                * quad(
                                    lambda energy: (
                                        (
                                            u0
                                            + (u1 - u0)
                                            * (energy - a)
                                            / (b - a)
                                        )
                                        * energy
                                        * np.exp(slope * energy)
                                    ),
                                    a,
                                    b,
                                    epsabs=1e-11,
                                    epsrel=1e-12,
                                )[0]
                            )
                            exact_rate += value
                            i, j = collision.i[interval], collision.j[interval]
                            expected[i] += 0.7 * collision.in_factor * value
                            expected[j] -= 0.7 * value
                        assert_allclose(
                            s.rate(f, collision), exact_rate, rtol=1e-12
                        )
                        assert_allclose(
                            s.rate(f, collision, weighted=True),
                            0.7 * exact_rate,
                            rtol=1e-12,
                        )
                        actual = s._PQ(f).dot(f)
                        error = np.linalg.norm(actual - expected, 1)
                        self.assertLess(error / exact_rate, 1e-12)


if __name__ == "__main__":
    unittest.main()
