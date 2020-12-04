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
from scipy.integrate import nquad, quad

from ries.nonresonant.klein_nishina import KleinNishina


@pytest.mark.parametrize("E", [(0.5), (1.0), (5.0), (10.0)])
def test_klein_nishina(E):
    # All the tests in this module are self consistent in the sense that the result of one method
    # of the KleinNishina class is tested against the result of another method.
    compton = KleinNishina()

    # Test the Compton-edge calculation by comparing `KleinNishina.compton_edge` and
    # `KleinNishina.Ep_over_E`.
    assert np.isclose(compton.compton_edge(E), E * compton.Ep_over_E(E, np.pi), 1e-5)
    # `KleinNishina.Ep_over_E` (multiplied by E) and KleinNishina.theta` should be the inverse
    # of each other.
    assert np.isclose(
        E * compton.Ep_over_E(E, compton.theta(E, 0.5 * E)), 0.5 * E, 1e-5
    )

    # Compare the analytical expression for the total cross section to numerical integrals of
    # various differential cross sections.
    cs_total_analytical = compton.cs_total(E)

    # Polarized
    cs_total_from_energy_differential = nquad(
        lambda Ep, phi: compton.cs_diff_dEp_dphi(E, Ep, phi),
        [(compton.compton_edge(E), E), (0.0, 2.0 * np.pi)],
    )[0]
    cs_total_from_solid_angle_differential = nquad(
        lambda theta, phi: compton.cs_diff_dOmega(E, theta, phi) * np.sin(theta),
        [(0.0, np.pi), (0.0, 2.0 * np.pi)],
    )[0]

    assert np.isclose(cs_total_analytical, cs_total_from_energy_differential, 1e-5)
    assert np.isclose(cs_total_analytical, cs_total_from_solid_angle_differential, 1e-5)

    # Unpolarized
    cs_total_from_energy_differential = quad(
        lambda Ep: compton.cs_diff_dEp(E, Ep), compton.compton_edge(E), E
    )[0]
    cs_total_from_scattering_angle_differential = quad(
        lambda theta: compton.cs_diff_dtheta(E, theta), 0.0, np.pi
    )[0]

    assert np.isclose(cs_total_analytical, cs_total_from_energy_differential, 1e-5)
    assert np.isclose(
        cs_total_analytical, cs_total_from_scattering_angle_differential, 1e-5
    )
