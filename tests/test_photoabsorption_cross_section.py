import numpy as np

from ries.resonance.photoabsorption_cross_section import PhotoabsorptionCrossSection
from ries.resonance.voigt import Voigt

from .boron import natural_boron, B10, B11

def test_photoabsorption_cross_section():
    pa_11b = PhotoabsorptionCrossSection(B11, Voigt, (B11.amu, 293.))

    assert len(pa_11b.resonances) == len(B11.excited_states)
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
        energy = resonance.probability_grid(0.9, 10)

        cs_resonance = resonance(energy)
        cs_isotope = pa_11b(energy)
        assert np.allclose(cs_isotope, cs_resonance, rtol=1e-5)