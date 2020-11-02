from scipy.constants import physical_constants

from ries.resonance.resonance import Resonance
from ries.resonance.recoil import FreeNucleusRecoil

from .boron import B11

def test_recoil():
    cs = Resonance(
        initial_state=B11.ground_state,
        intermediate_state=B11.excited_states['5/2^-_1'],
        recoil_correction=FreeNucleusRecoil(11.009305166)
    )

    energy_difference = cs.intermediate_state.excitation_energy - cs.initial_state.excitation_energy
    assert cs.resonance_energy == energy_difference*(1.+energy_difference/(2.*11.009305166*physical_constants['atomic mass constant energy equivalent in MeV'][0]))