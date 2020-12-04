# This file is part of ries.

# ries is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ries is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ries.  If not, see <https://www.gnu.org/licenses/>.

r"""
Multidimensional numerical integration with a partitioned first coordinate using `scipy.integrate.nquad`

The goal is to integrate an :math:`n_d`-dimensional function
:math:`f\left( \mathbf{x} \right)` on a hyperrectangle with the side lengths

.. math:: x_{d, \mathrm{f}} - x_{d, \mathrm{i}},

where :math:`x_d` (:math:`0 \leq d < n_d`) is the :math:`d`-th coordinate of
:math:`\mathbf{x}`, and the indices :math:`i` and :math:`f` denote the limits of this coordinate.
For the :math:`x_0` coordinate, a partition :math:`\left\{ x_{0, k} \right\}` (:math:`0 \leq k < n`, :math:`x_{0, 0} = x_{0, i}`, :math:`x_{0, n-1} = x_{0, f}`, and
:math:`x_{0,k} < x_{0, k+1}` for all `k < n - 1`) is assumed to be given.

This function calculates

.. math:: \sum_{k=0}^{n-2} \underbrace{\int_{x_{n_d-1, i}}^{x_{n_d-1, f}} ... \int_{x_{0, k}}^{x_{0, k+1}} f \left( \mathbf{x} \right) \mathrm{d} x_0 ... \mathrm{d} x_{n_d - 1}}_{I_i}.

Mathematically, this is equivalent to

.. math:: \int_{x_{n_d-1, i}}^{x_{n_d-1, f}} ... \int_{x_{0, i}}^{x_{0, f}} f \left( \mathbf{x} \right) \mathrm{d} x_0 ... \mathrm{d} x_{n_d - 1},

but for functions with narrow peaks in the interval :math:`\left[ x_{0,i}, x_{0,f} \right]`,
it may be advantageous to guide a numerical integration algorithm to the 'interesting' range.

The `quad_partition` function wraps the multidimensional integrator `scipy.integrate.nquad` and
integrates each on the :math:`x_0` axis separately.
Each `nquad` integration also reports an estimate of the absolute error
:math:`\left| \Delta I_i \right|` of the integral
:math:`I_i`.
The total error returned by `quad_partition` is:

.. math:: \sqrt{ \sum_{i=0}^{n-2} \left| \Delta I_i \right|^2 }.
"""

import numpy as np
from scipy.integrate import nquad


def quad_partition(f, x0, x1_to_xn_minus_1=None):
    r"""Multidimensional numerical integration with a partitioned first coordinate

    Parameters:

    - `f`, callable, multidimensional function :math:`f`.
    - `x0`, (1,n) array of float, :math:`\left\{ x_{0,k} \right\}`
    - `x1_to_xn_minus_1`, (n_d - 1, 2) array of float in which the l-th line contains the integration
      limits for the l+1-th coordinate (default: empty list, i.e. :math:`f` is a 1D function).

    Returns:

    - (float, float), integral and absolute error estimate.
    """
    integral = 0.0
    uncertainty_estimate = 0.0
    x1_to_xn_minus_1 = x1_to_xn_minus_1 or []

    for i in range(len(x0) - 1):
        nquad_result = nquad(f, [[x0[i], x0[i + 1]], *x1_to_xn_minus_1])
        integral += nquad_result[0]
        uncertainty_estimate += nquad_result[1] * nquad_result[1]

    return (integral, np.sqrt(uncertainty_estimate))
