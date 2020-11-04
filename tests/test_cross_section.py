import pytest

import numpy as np
from scipy.constants import physical_constants

from ries.cross_section import ConstantCrossSection, CrossSectionWeightedSum
from ries.constituents.element import natural_elements
from ries.nonresonant.xrmac import cm_to_fm, kg_to_g, xrmac
from ries.resonance.voigt import Voigt
from ries.resonance.debye_model import effective_temperature_debye_approximation, room_temperature_T_D

from .boron import B11

class TestCrossSection:
    def test_algebra(self):

        # Create array of all 11B ground-state transitions.
        ground_state_resonances = [
            Voigt(
                    B11.ground_state,
                    B11.excited_states[excited_state],
                    B11.amu, effective_temperature_debye_approximation(293., room_temperature_T_D['B'])
            ) for excited_state in B11.excited_states
        ]

        # Test CrossSection.__add__()
        photoabsorption_cross_section = ground_state_resonances[0] + ground_state_resonances[1]
        single_cross_section = ground_state_resonances[0](ground_state_resonances[0].resonance_energy)
        sum_cross_section = photoabsorption_cross_section(ground_state_resonances[0].resonance_energy)
        assert np.isclose(single_cross_section, sum_cross_section, rtol=1e-5)

        # Test CrossSection.__mul__() and CrossSection.__rmul__()
        photoabsorption_cross_section = ground_state_resonances[0] + 0.*ground_state_resonances[1]
        single_cross_section = ground_state_resonances[0](ground_state_resonances[0].resonance_energy)
        sum_cross_section = photoabsorption_cross_section(ground_state_resonances[0].resonance_energy)
        assert np.isclose(single_cross_section, sum_cross_section, rtol=1e-5)

        # Test CrossSection.__add__() using the convenient sum() function.
        # CrossSection.radd() are CrossSectionWeightedSum().add() are tested implicitly, 
        # because sum() adds all terms to an integer value of 0, and because any sum of cross 
        # sections yields a CrossSectionWeightedSum object.
        photoabsorption_cross_section = sum(ground_state_resonances)

        for resonance in ground_state_resonances:
            energy = resonance.equidistant_probability_grid(0.9, 10)
            single_cross_section = resonance(energy)
            sum_cross_section = photoabsorption_cross_section(energy)

            assert np.allclose(single_cross_section, sum_cross_section, rtol=1e-3)

        # Test CrossSectionWeightedSum.__add__() and CrossSectionWeightedSum.__mul__().
        photoabsorption_cross_section = photoabsorption_cross_section + 0.5*photoabsorption_cross_section

        for resonance in ground_state_resonances:
            energy = resonance.equidistant_probability_grid(0.9, 10)
            single_cross_section = resonance(energy)
            sum_cross_section = photoabsorption_cross_section(energy)

            assert np.allclose(1.5*single_cross_section, sum_cross_section, rtol=1e-3)

        # Test CrossSectionWeightedSum.__radd__().
        photoabsorption_cross_section = 0. + (ground_state_resonances[0] + ground_state_resonances[1])
        single_cross_section = ground_state_resonances[0](ground_state_resonances[0].resonance_energy)
        sum_cross_section = photoabsorption_cross_section(ground_state_resonances[0].resonance_energy)
        assert np.isclose(single_cross_section, sum_cross_section, rtol=1e-5)

        # Test additivity with a derived class for a nonresonant cross section.
        photoabsorption_cross_section = photoabsorption_cross_section + xrmac['B']

        photoabsorption_cross_section_value = 5.890e-2*cm_to_fm**2*natural_elements['B'].amu*physical_constants['atomic mass constant'][0]*kg_to_g
        assert np.isclose(photoabsorption_cross_section(1.), photoabsorption_cross_section_value, rtol=1e-6)

    def test_grid(self):
        assert np.allclose(
            xrmac['B'].equidistant_energy_grid((0., 1.), 3),
            np.array([0., 0.5, 1.]), rtol=1e-3
        )

        cross_section = ConstantCrossSection(1.)
        assert np.allclose(
            cross_section.equidistant_probability_grid((0., 1.), 3),
            np.array([0., 0.5, 1.]),
            rtol=1e-6
        )

        cross_section = cross_section + cross_section
        assert np.allclose(
            cross_section.equidistant_probability_grid((0., 1.), 3),
            np.array([0., 0.5, 1.]),
            rtol=1e-6
        )

        assert np.allclose(
            cross_section.equidistant_energy_grid((0., 1.), 3),
            np.array([0., 0.5, 1.]),
            rtol=1e-6
        )