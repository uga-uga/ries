# This file is part of ries.

# ries is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ries is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ries.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from scipy.constants import physical_constants
from scipy.stats import norm

from ries.resonance.resonance import Resonance

class Gauss(Resonance):
    def __init__(self, initial_state, intermediate_state, amu, effective_temperature,
        final_state=None):
        Resonance.__init__(self, initial_state, intermediate_state, final_state)

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