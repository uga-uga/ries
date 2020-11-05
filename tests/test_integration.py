import pytest

import numpy as np
from scipy.integrate import quad

from ries.integration.darboux import darboux
from ries.integration.quad_subintervals import quad_subintervals
from ries.resonance.gauss import Gauss

from .boron import B11

def stagger(x):
    fx = np.zeros(len(x))
    fx[::2] = 1.
    return fx

class TestIntegration:
    # Integrate a trivial function to see whether the integration algorithm works.
    @pytest.mark.parametrize('n_points', [(2), (10), (100), (1000)])
    def test_quad_subinterval(self, n_points):
        assert np.isclose(
            quad_subintervals(lambda x: 1., np.linspace(0., 1., n_points)),
            1., rtol = 1e-6
        )

    @pytest.mark.parametrize('n_points', [(2), (10), (100), (1000)])
    def test_darboux(self, n_points):
        integral = darboux(lambda x: 1. if np.ndim(x) == 0 else np.ones(len(x)), np.linspace(0., 1., n_points))
        assert integral[0] == integral[1]
        assert np.isclose(integral[0], 1., rtol = 1e-6)
        assert np.isclose(integral[1], 1., rtol = 1e-6)

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
        cross_section_integral_numerical = quad_subintervals(cross_section, energies)

        assert np.isclose(cross_section_integral_numerical, cross_section_integral_analytical, rtol=1e-3)

        # The simple Darboux lower sum needs more grid points.
        energies = cross_section.equidistant_probability_grid(limits, 10000)
        cross_section_integral_numerical = darboux(cross_section, energies)[0]

        assert np.isclose(cross_section_integral_numerical, cross_section_integral_analytical, rtol=1e-3)