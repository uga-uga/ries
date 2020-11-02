import pytest

import numpy as np
from scipy.constants import physical_constants

from ries.cross_section import SumCrossSection
from ries.constituents.element import natural_elements
from ries.nonresonant.xrmac import cm_to_fm, kg_to_g, xrmac
from ries.resonance.voigt import Voigt

from .boron import B11

class TestCrossSection:
    def test_addition(self):

        ground_state_resonances = [
            Voigt(
                    B11.ground_state,
                    B11.excited_states[excited_state],
                    B11.amu, 293.
            ) for excited_state in B11.excited_states
        ]

        # This does not work in python, because the sum() function is intended for numerical,
        # input.
        # It adds everything to an integer value of zero, and this causes a TypeError because 
        # int + CrossSection is not defined.
        with pytest.raises(TypeError):
            photoabsorption_cross_section = sum(ground_state_resonances)

        photoabsorption_cross_section = ground_state_resonances[0]
        for resonance in ground_state_resonances[1:]:
            photoabsorption_cross_section += resonance
        assert len(photoabsorption_cross_section.reactions) == len(ground_state_resonances)

        for resonance in ground_state_resonances:
            energy = resonance.equidistant_probability_grid(0.9, 10)
            single_cross_section = resonance(energy)
            sum_cross_section = photoabsorption_cross_section(energy)

            assert np.allclose(single_cross_section, sum_cross_section, rtol=1e-3)

        photoabsorption_cross_section = photoabsorption_cross_section + photoabsorption_cross_section
        assert len(photoabsorption_cross_section.reactions) == 2*len(ground_state_resonances)

        for resonance in ground_state_resonances:
            energy = resonance.equidistant_probability_grid(0.9, 10)
            single_cross_section = resonance(energy)
            sum_cross_section = photoabsorption_cross_section(energy)

            assert np.allclose(single_cross_section, 0.5*sum_cross_section, rtol=1e-3)

        photoabsorption_cross_section = photoabsorption_cross_section + xrmac['B']

        photoabsorption_cross_section_value = 5.890e-2*cm_to_fm**2*natural_elements['B'].amu*physical_constants['atomic mass constant'][0]*kg_to_g
        assert np.isclose(photoabsorption_cross_section(1.), photoabsorption_cross_section_value, rtol=1e-6)

    def test_grid(self):
        assert np.allclose(
            xrmac['B'].equidistant_energy_grid((0., 1.), 3),
            np.array([0., 0.5, 1.]), rtol=1e-3
        )