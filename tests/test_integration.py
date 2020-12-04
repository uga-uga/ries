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
from scipy.integrate import quad

from ries.integration.darboux import darboux
from ries.integration.quad_partition import quad_partition
from ries.resonance.gauss import Gauss

from .boron import B11

def stagger(x):
    fx = np.zeros(len(x))
    fx[::2] = 1.
    return fx

class TestIntegration:
    # Integrate a trivial function to see whether the integration algorithm works.
    @pytest.mark.parametrize('n_points', [(2), (10), (100), (1000)])
    def test_quad_partition(self, n_points):
        assert np.isclose(
            quad_partition(lambda x: 1., np.linspace(0., 1., n_points))[0],
            1., rtol = 1e-6
        )

    @pytest.mark.parametrize('n_points', [(2), (10), (100), (1000)])
    def test_nquad_partition(self, n_points):
        assert np.isclose(
            quad_partition(lambda x, y: 1., np.linspace(0., 1., n_points), [[0., 1.]])[0],
            1., rtol = 1e-6
        )

    @pytest.mark.parametrize('n_points', [(2), (10), (100), (1000)])
    def test_darboux(self, n_points):
        # Call darboux by passing the function to be evaluated.
        integral = darboux(lambda x: 1. if np.ndim(x) == 0 else np.ones(len(x)), np.linspace(0., 1., n_points))
        assert np.isclose(integral[0], 1., rtol = 1e-6)
        assert np.isclose(1.+integral[1], 1., rtol = 1e-6)

        # Call darboux by passing the values of the function at the points x.
        assert np.isclose(darboux(np.ones(n_points), np.linspace(0., 1., n_points))[0], 1., rtol=1e-6)

    def test_darboux_lower_upper(self):
        n_points = 100
        x = np.arange(n_points)
        integral = darboux(stagger, x)
        assert integral[0] == 0.
        assert integral[1] == n_points-1.

    def test_resonance_integration(self):
        # Integrate a resonance shape over an extremely large energy range.
        cross_section = Gauss(
            B11.ground_state,
            B11.excited_states['1/2^-_1'],
            B11.amu,
            1.
        )

        limits = (1., 3.)

        # Within this large energy range, scipy.integrate.quad is not able to find the narrow 
        # resonance and returns a value of zero.
        cross_section_integral_analytical = cross_section.energy_integrated_cross_section
        cross_section_integral_numerical = quad(cross_section, *limits)[0]

        assert not np.isclose(cross_section_integral_numerical, cross_section_integral_analytical, rtol=1e-1)

        # The problem can be solved by integrating over subintervals whose lenghts are adapted to
        # the function at hand.
        energies = cross_section.equidistant_probability_grid(limits, 2000)
        cross_section_integral_numerical = quad_partition(cross_section, energies)[0]

        assert np.isclose(cross_section_integral_numerical, cross_section_integral_analytical, rtol=1e-3)

        # The simple Darboux lower sum needs more grid points.
        energies = cross_section.equidistant_probability_grid(limits, 10000)
        cross_section_integral_numerical = darboux(cross_section, energies)[0]

        assert np.isclose(cross_section_integral_numerical, cross_section_integral_analytical, rtol=1e-3)

    def test_resonance_integration_2d(self):
        cross_section = Gauss(
            B11.ground_state,
            B11.excited_states['1/2^-_1'],
            B11.amu,
            1.
        )

        limits = (1., 3.)

        energies = cross_section.equidistant_probability_grid(limits, 250)
        cross_section_integral_analytical = cross_section.energy_integrated_cross_section
        cross_section_integral_numerical = quad_partition(lambda e, z: cross_section(e), energies, [[0., 1.]])[0]

        assert np.isclose(cross_section_integral_numerical, cross_section_integral_analytical, rtol=1e-2)
