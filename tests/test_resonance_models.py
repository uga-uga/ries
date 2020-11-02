import numpy as np
import pytest
from scipy.integrate import quad

from ries.resonance.breit_wigner import BreitWigner
from ries.resonance.gauss import Gauss
from ries.resonance.pseudo_voigt import PseudoVoigt
from ries.resonance.voigt import Voigt

from .boron import B11

class TestResonanceModels:
    @pytest.mark.parametrize('Model, parameters, rtol', 
    [
        (BreitWigner, [], 1e-4),
        (Gauss, [B11.amu, 293.], 1e-4),
        (PseudoVoigt, [B11.amu, 293.], 1e-3),
        (Voigt, [B11.amu, 293.], 2e-2),
    ])
    def test_coverage(self, Model, parameters, rtol):
        cs = Model(
            B11.ground_state,
            B11.excited_states['5/2^-_1'],
            *parameters
        )
        cov_int = cs.coverage_interval(0.5)
        assert np.isclose(quad(cs, cov_int[0], cov_int[1])[0], 0.5*cs.energy_integrated_cross_section, rtol=rtol)

    def test_limits(self):
        voigt = Voigt(
            B11.ground_state,
            B11.excited_states['5/2^-_1'],
            B11.amu,
            0.
        )
        energy = voigt.equidistant_probability_grid(0.95, 100)
        assert np.allclose(voigt.equidistant_probability_grid((energy[0], energy[-1]), 100), energy, rtol=1e-5)

        breit_wigner = BreitWigner(
            B11.ground_state,
            B11.excited_states['5/2^-_1'],
        )
        assert np.allclose(voigt(energy), breit_wigner(energy), rtol=1e-2)

        voigt = Voigt(
            B11.ground_state,
            B11.excited_states['5/2^-_1'],
            B11.amu,
            1e4
        )
        energy = voigt.equidistant_probability_grid(0.95, 100)

        gauss = Gauss(
            B11.ground_state,
            B11.excited_states['5/2^-_1'],
            B11.amu,
            1e4
        )
        assert np.allclose(voigt(energy), gauss(energy), rtol=5e-2)

        pseudo_voigt = PseudoVoigt(
            B11.ground_state,
            B11.excited_states['5/2^-_1'],
            B11.amu,
            1e4
        )

        assert np.allclose(voigt(energy), pseudo_voigt(energy), rtol=1e-2)
        assert np.allclose(pseudo_voigt.equidistant_probability_grid((energy[0], energy[-1]), 100), energy, rtol=1e-5)
