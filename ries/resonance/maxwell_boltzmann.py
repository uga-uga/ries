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

"""Maxwell-Boltzmann distribution"""

import numpy as np

from scipy.constants import physical_constants


class MaxwellBoltzmann:
    def __init__(self, amu, effective_temperature):
        self.amu = amu
        self.effective_temperature = effective_temperature

    def get_doppler_width(self, resonance_energy):
        return resonance_energy * np.sqrt(
            2.0
            * physical_constants["Boltzmann constant in eV/K"][0]
            * 1e-6
            * self.effective_temperature
            / (
                self.amu
                * physical_constants["atomic mass constant energy equivalent in MeV"][0]
            )
        )
