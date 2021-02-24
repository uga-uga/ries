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

r"""Maxwell-Boltzmann distribution

The Maxwell-Boltzmann distribution does not only describes the momentum(:math:`p`) distribution 
of particles in an ideal gas, but many other systems like a particle in a harmonic potential 
exhibit normal-distributed [:math:`\exp \left( - c p^2 \right)` with a constant 
:math:`c`] momenta :cite:`LandauLifshitz1980`.

At the moment, this module only provides the distribution of Doppler shifts in an ideal gas 
(see `ries.resonance.voigt` for more information) 
obtained from the 1D ideal-gas expression for the velocity distribution of atoms 
(:math:`\approx` nuclei) with a mass :math:`m \left( ^A\mathrm{X}\right)` at a temperature 
:math:`T_\mathrm{eff}`:

.. math:: P^\mathrm{MB} \left( v \right) = \sqrt{\frac{m \left( ^A\mathrm{X} \right)}{2 \pi k_B T_\mathrm{eff}}} \exp \left[ - \left( \frac{m \left( ^A\mathrm{X}\right) v}{2 k_B T_\mathrm{eff}} \right)^2 \right].

Here, :math:`k_B` denotes the Boltzmann constant.
By using an effective temperature instead of the thermodynamical temperature of the ensemble,
effects of a potential can be taken into account (see, e.g. Refs. :cite:`Lamb1939` :cite:`Metzger1959`).
"""

import numpy as np

from scipy.constants import physical_constants


class MaxwellBoltzmann:
    """Maxwell-Boltzmann distribution of particles in an ideal gas

    Attributes:

    - `amu`, float, mass of the particle in atomic mass units.
    - `effective_temperature`, float, effective temperature of the ensemble of particle in K.
    """

    def __init__(self, amu, effective_temperature):
        """Initialization

        Parameters:

        - `amu`, float, mass of the particle in atomic mass units.
        - `effective_temperature`, float, effective temperature of the ensemble of particle in K.
        """
        self.amu = amu
        self.effective_temperature = effective_temperature

    def get_doppler_width(self, E):
        r"""Doppler width in ensemble with a Maxwell-Boltzmann velocity distribution

        The Doppler width :math:`\Delta` is the square root of two times the standard deviation
        :math:`\sigma` of a normal distribution (:math:`\Delta = \sqrt{2} \sigma`) with the
        reference energy :math:`E`.

        Parameters:

        - `E`, float or array_like, reference energy in MeV.

        Returns:

        - float or array_like, Doppler width in MeV.
        """
        return E * np.sqrt(
            2.0
            * physical_constants["Boltzmann constant in eV/K"][0]
            * 1e-6
            * self.effective_temperature
            / (
                self.amu
                * physical_constants["atomic mass constant energy equivalent in MeV"][0]
            )
        )

    @staticmethod
    def get_effective_temperature(doppler_width, amu, E):
        r"""Return doppler width for a given effective temperature

        This convenience function is the inverse of `MaxwellBoltzmann.get_doppler_width`.
        In the literature, in particular for schematic plots, one often finds that the Doppler
        width is given instead of the effective temperature, since it can directly be compared to
        the total width of a resonance.
        Metzger :cite:`Metzger1959`, for example, makes the relative magnitude of :math:`\Delta`
        :math:`Gamma` a criterion to decide whether the Gaussian 'Doppler form' of the resonance
        shape can be used.

        Parameters
        ----------
        - `doppler_width`, float or array_like, Doppler width :math:`\Delta` in MeV.
        - `amu`, float, mass of the particle in atomic mass units.
        - `E`, float or array_like, reference energy in MeV.

        Returns
        -------

        - float or array_like, Effective temperature in K
        """

        return (
            doppler_width*doppler_width
            * amu
            * physical_constants["atomic mass constant energy equivalent in MeV"][0]
        ) / (E * E * 2.0 * physical_constants["Boltzmann constant in eV/K"][0] * 1e-6)
