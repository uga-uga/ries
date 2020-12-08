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

r"""Doppler-broadened Breit-Wigner resonance approximated using a normal distribution.

In the limit where the Doppler width :math:`\Delta` is much larger than the width :math:`\Gamma_2` 
of the excited state,

.. math:: \Delta \gg \Gamma_2,

the convolution of the Cauchy distribution and the normal distribution asymptotically approaches 
shape of the latter.
The condition above is often a good approximation for nuclei with a strongly fragmented 
low-energy strength function at standard temperatures, and accurate and efficient numerical 
algorithms exist to calculate the PDF, CDF, and PPF of a normal distribution 
(as compared to the Voigt distribution, which would be the exact result of the aforementioned 
convolution).

For more information, see `ries.resonance.voigt`.
"""

import numpy as np
from scipy.constants import physical_constants
from scipy.stats import norm

from ries.resonance.resonance import Resonance
from ries.resonance.maxwell_boltzmann import MaxwellBoltzmann


class Gauss(Resonance):
    r"""Approximation for a Doppler-broadened Breit-Wigner cross section (normal distribution)

    See `ries.resonance.resonance.Resonance`.
    """
    def __init__(
        self,
        initial_state,
        intermediate_state,
        amu,
        effective_temperature,
        final_state=None,
    ):
        r"""Initialization

        Parameters:

        - `amu`, float, mass of the nucleus in atomic mass units.
        - `effective_temperature`, float, effective temperature of the ensemble of nuclei in K.

        See also `ries.resonance.resonance.Resonance`.
        """
        Resonance.__init__(self, initial_state, intermediate_state, final_state)

        self.maxwell_boltzmann = MaxwellBoltzmann(amu, effective_temperature)

        self.probability_distribution = norm
        self.probability_distribution_parameters = (
            self.resonance_energy,
            self.maxwell_boltzmann.get_doppler_width(self.resonance_energy)
            / np.sqrt(2.0),
        )
