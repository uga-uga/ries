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

import pytest

import numpy as np

from scipy.stats import cauchy, norm

from ries.resonance.pseudo_voigt import PseudoVoigtDistribution

# The PPF of the Voigt profile is approximated by the inverse of the CDF of a pseudo-Voigt profile.
# Using `scipy.optimize.newton`, a root-finder algorithm, the code first tried to perform the
# inversion numerically.
# `newton` returns different output depending on whether a scalar or an array of start values was
# passed.
# If `newton` does not report convergence, a weighted average of the PPFs of the normal and the
# Cauchy distribution is returned.
# This test explores all different calls and outputs of `newton`.
def test_voigt_ppf():
    resonance_energy = 1e6
    pseudo_voigt = PseudoVoigtDistribution(resonance_energy, 1.0, 1.0, 100.0)

    # Test scalar and array input.
    assert np.isclose(pseudo_voigt.ppf(0.5), resonance_energy, 1e-5)
    assert np.allclose(
        pseudo_voigt.ppf(np.array([0.5, 0.5])),
        np.array([resonance_energy, resonance_energy]),
        1e-5,
    )

    # Test fallback approximation.
    # The result should be somewhere in between the results of the two constituent PPFs.
    extremely_small_quantile = 1e-7
    with pytest.warns(UserWarning):
        assert pseudo_voigt.ppf(extremely_small_quantile) > cauchy.ppf(
            extremely_small_quantile
        )
        assert pseudo_voigt.ppf(extremely_small_quantile) < norm.ppf(
            extremely_small_quantile
        )
