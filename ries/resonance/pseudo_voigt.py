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
from scipy.stats import cauchy, norm

from ries.resonance.maxwell_boltzmann import MaxwellBoltzmann
from ries.resonance.resonance import Resonance


class PseudoVoigtDistribution:
    def __init__(self, resonance_energy, width, amu, effective_temperature):
        self.resonance_energy = resonance_energy
        self.maxwell_boltzmann = MaxwellBoltzmann(amu, effective_temperature)

        self.Gamma_L = width
        self.Gamma_L_squared = self.Gamma_L * self.Gamma_L
        self.Gamma_L_cubed = self.Gamma_L_squared * self.Gamma_L
        self.Gamma_L_squared_squared = self.Gamma_L_squared * self.Gamma_L_squared
        self.doppler_width = self.maxwell_boltzmann.get_doppler_width(
            self.resonance_energy
        )
        self.sigma = self.doppler_width / np.sqrt(2.0)
        self.Gamma_G = 2.0 * np.sqrt(np.log(2.0)) * self.doppler_width
        self.Gamma_G_squared = self.Gamma_G * self.Gamma_G
        self.Gamma_G_squared_squared = self.Gamma_G_squared * self.Gamma_G_squared

        self.Gamma = self.get_Gamma()
        self.eta = self.get_eta()
        self.gamma_G = self.Gamma / (2.0 * np.sqrt(np.log(2.0)))
        self.gamma_L = 0.5 * self.Gamma

        self.gamma = 0.5 * width

    def cdf(self, energy):
        return (1.0 - self.eta) * norm.cdf(
            energy, loc=self.resonance_energy, scale=self.gamma_G / np.sqrt(2.0)
        ) + self.eta * cauchy.cdf(
            energy, loc=self.resonance_energy, scale=0.5 * self.gamma_L
        )

    def pdf(self, energy):
        return (1.0 - self.eta) * norm.pdf(
            energy, loc=self.resonance_energy, scale=self.gamma_G / np.sqrt(2.0)
        ) + self.eta * cauchy.pdf(
            energy, loc=self.resonance_energy, scale=0.5 * self.gamma_L
        )

    def ppf(self, quantile):
        return (1.0 - self.eta) * norm.ppf(
            quantile, loc=self.resonance_energy, scale=self.gamma_G / np.sqrt(2.0)
        ) + self.eta * cauchy.ppf(
            quantile, loc=self.resonance_energy, scale=0.5 * self.gamma_L
        )

    def get_Gamma(self):
        return (
            self.Gamma_G_squared_squared * self.Gamma_G
            + 2.69269 * self.Gamma_G_squared_squared * self.Gamma_L
            + 2.42843 * self.Gamma_G_squared * self.Gamma_G * self.Gamma_L_squared
            + 4.47163 * self.Gamma_G_squared * self.Gamma_L_cubed
            + 0.07842 * self.Gamma_G * self.Gamma_L_squared_squared
            + self.Gamma_L_squared_squared * self.Gamma_L
        ) ** 0.2

    def get_eta(self):
        inverse_Gamma = 1.0 / self.Gamma
        inverse_Gamma_squared = inverse_Gamma * inverse_Gamma
        return (
            1.36603 * self.Gamma_L * inverse_Gamma
            - 0.47719 * self.Gamma_L_squared * inverse_Gamma_squared
            + 0.11116 * self.Gamma_L_cubed * inverse_Gamma_squared * inverse_Gamma
        )


class PseudoVoigt(Resonance):
    def __init__(
        self,
        initial_state,
        intermediate_state,
        amu,
        effective_temperature,
        final_state=None,
    ):
        Resonance.__init__(self, initial_state, intermediate_state, final_state)

        self.probability_distribution = PseudoVoigtDistribution(
            self.resonance_energy,
            self.intermediate_state.width,
            amu,
            effective_temperature,
        )
        self.probability_distribution_parameters = ()
