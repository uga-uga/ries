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

import pytest

import numpy as np
from scipy.constants import physical_constants

from ries.resonance.maxwell_boltzmann import MaxwellBoltzmann

def test_doppler_width():

    maxwell_boltzmann = MaxwellBoltzmann(2.0 * physical_constants["Boltzmann constant in eV/K"][0]*1e-6/physical_constants["atomic mass constant energy equivalent in MeV"][0], 1.)

    assert maxwell_boltzmann.get_doppler_width(2.) == 2.
    assert maxwell_boltzmann.get_effective_temperature(maxwell_boltzmann.get_doppler_width(1.), 2.0 * physical_constants["Boltzmann constant in eV/K"][0]*1e-6/physical_constants["atomic mass constant energy equivalent in MeV"][0], 2.) == 0.25
