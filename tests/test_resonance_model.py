import numpy as np
from scipy.constants import physical_constants

from ries.resonance.resonance_model import ResonanceModel

from .boron import B11

def test_resonance_model():
    cs = ResonanceModel(
        initial_state=B11.ground_state,
        intermediate_state=B11.excited_states['5/2^-_1']
    )

    assert cs(4.44498) == cs.energy_integrated_cross_section
    assert cs(0., input_is_absolute_energy=False) == cs.energy_integrated_cross_section
    assert cs.resonance_energy == cs.intermediate_state.excitation_energy - cs.initial_state.excitation_energy
    assert cs.final_state_branching_ratio == 1.
    assert cs.statistical_factor == 1.5
    assert np.allclose(cs.energy_grid(0.5, 3), np.array([4.44498-0.25, 4.44498, 4.44498+0.25]), rtol=1e-5)
    assert np.allclose(cs.probability_grid(0.5, 3), np.array([4.44498-0.25, 4.44498, 4.44498+0.25]), rtol=1e-5)

    cs = ResonanceModel(
        initial_state=B11.ground_state,
        intermediate_state=B11.excited_states['3/2^-_2'],
        final_state=B11.excited_states['1/2^-_1']
    )

    assert cs.resonance_energy == cs.intermediate_state.excitation_energy - cs.initial_state.excitation_energy
    assert cs.final_state_branching_ratio == 0.144
    assert cs.statistical_factor == 1.

    energy_integrated_cross_section = (
        np.pi*np.pi
        *physical_constants['reduced Planck constant times c in MeV fm'][0]**2
        /(5.02030**2)
        *(3.+1.)/(3.+1.)
        *0.856*1.97e-6
        *0.144
    )
    assert cs.energy_integrated_cross_section == energy_integrated_cross_section