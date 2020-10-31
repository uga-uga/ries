import pytest

import numpy as np
from scipy.integrate import quad

from ries.nonresonant.klein_nishina import compton_edge, cs_diff, cs_diff_de, cs_diff_dtheta, cs_total

@pytest.mark.parametrize('e0',
[
    (0.5), (1.), (5.), (10.)
])
def test_klein_nishina(e0):
    cs_total_analytical = cs_total(e0)
    cs_total_from_energy_differential = quad(lambda energy: cs_diff_de(e0, energy), compton_edge(e0), e0)[0]
    cs_total_from_polar_angle_differential = quad(lambda theta: cs_diff_dtheta(e0, theta), 0., np.pi)[0]

    assert np.isclose(cs_total_analytical, cs_total_from_energy_differential, 1e-5)
    assert np.isclose(cs_total_analytical, cs_total_from_polar_angle_differential, 1e-5)