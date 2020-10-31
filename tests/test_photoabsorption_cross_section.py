import numpy as np

from ries.resonance.photoabsorption_cross_section import PhotoabsorptionCrossSection
from ries.resonance.voigt import Voigt

from .B11 import B11

def test_photoabsorption_cross_section():
    photo_abs = PhotoabsorptionCrossSection(B11, Voigt, (B11.amu, 293.))
    assert len(photo_abs.resonances) == len(B11.excited_states)
    first_exc = Voigt(
        B11.ground_state,
        B11.excited_states['5/2^-_1'],
        B11.amu, 293.)
    first_exc = Voigt(
        B11.ground_state,
        B11.excited_states['5/2^-_1'],
        B11.amu, 293.)

    for state in B11.excited_states:
        resonance = Voigt(
            B11.ground_state,
            B11.excited_states[state],
            B11.amu, 293.)
        energy = resonance.probability_grid(0.9, 3)

        cs_isotope = photo_abs(energy)
        cs_resonance = resonance(energy)
        assert np.allclose(cs_isotope, cs_resonance, rtol=1e-5)
