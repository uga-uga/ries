import numpy as np

from ries.constituents.element import Element
from ries.resonance.photoabsorption_cross_section import PhotoabsorptionCrossSectionIsotope, PhotoabsorptionCrossSectionElement
from ries.resonance.voigt import Voigt

from .boron import natural_boron, B10, B11

def test_photoabsorption_cross_section():
    pa_10b = PhotoabsorptionCrossSectionIsotope(B10, Voigt, (B11.amu, 293.))
    pa_11b = PhotoabsorptionCrossSectionIsotope(B11, Voigt, (B11.amu, 293.))
    pa_b = PhotoabsorptionCrossSectionElement(
        Element(5, 'B',
            {'10B': B10, '11B': B11},
            {'10B': 0.5, '11B': 0.5}
        ),
        {'10B': Voigt, '11B': Voigt},
        {'10B': (B10.amu, 293.), '11B': (B11.amu, 293.)}
    )

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
        cs_element = pa_b(energy)
        assert np.allclose(cs_element, cs_resonance*pa_b.element.abundances['11B'], rtol=1e-5)