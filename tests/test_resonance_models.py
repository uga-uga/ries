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

import numpy as np
import pytest
from scipy.integrate import quad

from ries.resonance.breit_wigner import BreitWigner
from ries.resonance.gauss import Gauss
from ries.resonance.pseudo_voigt import PseudoVoigt
from ries.resonance.voigt import Voigt

from .boron import B11
from ries.resonance.debye_model import (
    effective_temperature_debye_approximation,
    room_temperature_T_D,
)


class TestResonanceModels:
    @pytest.mark.parametrize(
        "Model, parameters, rtol",
        [
            (BreitWigner, [], 1e-4),
            (
                Gauss,
                [
                    B11.amu,
                    effective_temperature_debye_approximation(
                        293.0, room_temperature_T_D["B"]
                    ),
                ],
                1e-4,
            ),
            (
                PseudoVoigt,
                [
                    B11.amu,
                    effective_temperature_debye_approximation(
                        293.0, room_temperature_T_D["B"]
                    ),
                ],
                1e-3,
            ),
            (
                Voigt,
                [
                    B11.amu,
                    effective_temperature_debye_approximation(
                        293.0, room_temperature_T_D["B"]
                    ),
                ],
                2e-2,
            ),
        ],
    )
    def test_coverage(self, Model, parameters, rtol):
        cs = Model(B11.ground_state, B11.excited_states["5/2^-_1"], *parameters)
        if Model in (PseudoVoigt, Voigt):
            with pytest.warns(UserWarning):
                cov_int = cs.coverage_interval(0.5)
        else:
            cov_int = cs.coverage_interval(0.5)

        assert np.isclose(
            quad(cs, cov_int[0], cov_int[1])[0],
            0.5 * cs.energy_integrated_cross_section,
            rtol=rtol,
        )

    def test_limits(self):
        voigt = Voigt(B11.ground_state, B11.excited_states["5/2^-_1"], B11.amu, 1e-3)
        with pytest.warns(UserWarning):
            energy = voigt.equidistant_probability_grid(0.95, 100)
            assert np.allclose(
                voigt.equidistant_probability_grid((energy[0], energy[-1]), 100),
                energy,
                rtol=1e-5,
            )

        breit_wigner = BreitWigner(
            B11.ground_state,
            B11.excited_states["5/2^-_1"],
        )
        assert np.allclose(voigt(energy), breit_wigner(energy), rtol=1e-2)

        voigt = Voigt(B11.ground_state, B11.excited_states["5/2^-_1"], B11.amu, 1e4)
        with pytest.warns(UserWarning):
            energy = voigt.equidistant_probability_grid(0.95, 100)

        gauss = Gauss(B11.ground_state, B11.excited_states["5/2^-_1"], B11.amu, 1e4)
        assert np.allclose(voigt(energy), gauss(energy), rtol=5e-2)

        pseudo_voigt = PseudoVoigt(
            B11.ground_state, B11.excited_states["5/2^-_1"], B11.amu, 1e4
        )

        assert np.allclose(voigt(energy), pseudo_voigt(energy), rtol=1e-2)
        with pytest.warns(UserWarning):
            assert np.allclose(
                pseudo_voigt.equidistant_probability_grid((energy[0], energy[-1]), 100),
                energy,
                rtol=1e-5,
            )
