"""Legacy Simpson quadrature for the solver's one-dimensional arrays.

Preserve scipy.integrate.simps(..., even='avg'), including on nonuniform
grids. Modern scipy.integrate.simpson uses a different even-sample rule.
"""
import numpy as np


def _odd_simpson(y, x):
    h = np.diff(x)
    h0, h1 = h[0::2], h[1::2]
    hsum = h0 + h1
    hprod = h0 * h1
    ratio = np.true_divide(h0, h1, out=np.zeros_like(h0), where=h1 != 0)
    inverse = np.true_divide(1.0, ratio, out=np.zeros_like(ratio), where=ratio != 0)
    middle = np.true_divide(hsum, hprod, out=np.zeros_like(hsum), where=hprod != 0)
    return np.sum(hsum / 6.0 * (y[0:-2:2] * (2.0 - inverse)
                               + y[1:-1:2] * hsum * middle
                               + y[2::2] * (2.0 - ratio)))


def simps(y, *, x):
    """Integrate 1-D samples with the original Simpson/averaged-end rule."""
    y, x = np.asarray(y), np.asarray(x, dtype=float)
    if y.ndim != 1 or x.ndim != 1 or y.shape != x.shape or len(y) < 2:
        raise ValueError('simps requires matching 1-D arrays of at least two samples')
    if len(y) % 2:
        return _odd_simpson(y, x)
    first = 0.5 * (x[1] - x[0]) * (y[1] + y[0])
    last = 0.5 * (x[-1] - x[-2]) * (y[-1] + y[-2])
    if len(y) == 2:
        return first
    return 0.5 * (first + last + _odd_simpson(y[:-1], x[:-1])
                  + _odd_simpson(y[1:], x[1:]))
