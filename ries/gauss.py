import numpy as np
from scipy.constants import physical_constants
from scipy.stats import norm

from ries.cross_section_model import CrossSectionModel

class Gauss(CrossSectionModel):
    def __init__(self, initial_state, intermediate_state, amu, effective_temperature,
        final_state=None):
        CrossSectionModel.__init__(self, initial_state, intermediate_state, final_state)

        self.amu = amu
        self.effective_temperature = effective_temperature
        self.doppler_width = self.get_doppler_width()

        self.probability_distribution = norm
        self.probability_distribution_parameters = (
            self.resonance_energy,
            self.doppler_width/np.sqrt(2.)
        )

    def get_doppler_width(self):
        return self.resonance_energy*np.sqrt(
            2.*physical_constants['Boltzmann constant in eV/K'][0]*1e-6*self.effective_temperature/(
                self.amu*physical_constants['atomic mass constant energy equivalent in MeV'][0]
            )
        )