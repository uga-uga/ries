"""
1D integration using an approximation of the Darboux integral

Given a 1D function :math:`f\\left( x \\right)` and a partition \
:math:`\\left\{ x_i \\right\}` (:math:`0 \\leq i < n` and \
:math:`x_i < x_{i+1}` for all `i < n - 1`), the upper and lower Darboux sums \
:math:`U` and :math:`L` are defined as :cite:`Cortzen2020` [b]_:

.. math:: U = \sum_{i = 0}^{n - 1} \\sup_{x \in \\left[ x_i, x_{i+1} \\right]} f \\left( x \\right) \\left( x_{i+1} - x_{i} \\right)

and

.. math:: L = \sum_{i = 0}^{n - 1} \\inf_{x \in \\left[ x_i, x_{i+1} \\right]} f \\left( x \\right) \\left( x_{i+1} - x_{i} \\right),

where :math:`\\sup` and :math:`\\inf` denote the supremum and the infimum in the given range.
For any partition, the upper and lower sums represent an upper and a lower limit for the integral:

.. math:: \\int_{x_0}^{x_{n-1}} f \\left( x \\right) \\mathrm{d} x,

i.e. it would be straightforward to obtain an estimate of the numerical uncertainty.

However, without knowledge of the entire function :math:`f`, the procedure above cannot be \
implemented as a numerical algorithm.

Here, the following approximations are made:

.. math:: \\sup_{x \in \\left[ x_i, x_{i+1} \\right]} f \\left( x \\right) \\approx \max \\left[ f \\left( x_i \\right), f \\left( x_{i+1} \\right)\\right]

.. math:: \\inf_{x \in \\left[ x_i, x_{i+1} \\right]} f \\left( x \\right) \\approx \min \\left[ f \\left( x_i \\right), f \\left( x_{i+1} \\right)\\right].

If the partition is fine enough so that :math:`f` is monotonous in all intervals, the \
approximations give the same results as the Darboux sums.

.. [b] Using :math:`\\alpha \\left( x \\right) = x` in the notation of Ref. :cite:`Cortzen2020`.
"""

import numpy as np

def darboux(f_or_fx, x):
    """Approximate lower and upper Darboux sums assuming a monotonous function
    
Parameters:

- `f_or_fx`, callable or array of float, either the function :math:`f`, or a set of values \
:math:`\\left\{ f \\left( x_i \\right) \\right\}`.
- `x`, array of float, :math:`\\left\{ x_i \\right\}`

Returns:

- (float, float), the value of the lower sum in this approximation, and the absolute value of the \
difference between the lower and the upper sum.
    """
    dx = (x - np.roll(x, 1))[1:]
    if callable(f_or_fx):
        fx = f_or_fx(x)
    else:
        fx = f_or_fx
    lower_upper = np.zeros((2, len(fx)-1))
    lower_upper[0] = fx[:-1]
    lower_upper[1] = fx[1:]
    lower_sum = np.sum(np.min(lower_upper, axis=0)*dx)
    upper_sum = np.sum(np.max(lower_upper, axis=0)*dx)
    return (
        lower_sum,
        abs(upper_sum-lower_sum),
    )