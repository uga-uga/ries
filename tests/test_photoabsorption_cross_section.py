import numpy as np

from ries.resonance.photoabsorption_cross_section import PhotoabsorptionCrossSection
from ries.resonance.voigt import Voigt

from .boron import natural_boron, B10, B11
from ries.resonance.debye_model import effective_temperature_debye_approximation, room_temperature_T_D

def test_photoabsorption_cross_section():
    pa_11b = PhotoabsorptionCrossSection(B11, Voigt, (B11.amu, effective_temperature_debye_approximation(293., room_temperature_T_D['B'])))

    assert len(pa_11b.resonances) == len(B11.excited_states)
    first_exc = Voigt(
        B11.ground_state,
        B11.excited_states['5/2^-_1'],
        B11.amu, effective_temperature_debye_approximation(293., room_temperature_T_D['B']))
    first_exc = Voigt(
        B11.ground_state,
        B11.excited_states['5/2^-_1'],
        B11.amu,effective_temperature_debye_approximation(293., room_temperature_T_D['B']))

    for state in B11.excited_states:
        resonance = Voigt(
            B11.ground_state,
            B11.excited_states[state],
            B11.amu, effective_temperature_debye_approximation(293., room_temperature_T_D['B']))
        energy = resonance.equidistant_probability_grid(0.9, 10)

        cs_resonance = resonance(energy)
        cs_isotope = pa_11b(energy)
        assert np.allclose(cs_isotope, cs_resonance, rtol=1e-5)