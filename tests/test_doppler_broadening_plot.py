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

import matplotlib.pyplot as plt

from ries.constituents.element import natural_elements
from ries.constituents.state import GroundState, State
from ries.resonance.breit_wigner import BreitWigner
from ries.resonance.gauss import Gauss
from ries.resonance.maxwell_boltzmann import MaxwellBoltzmann
from ries.resonance.voigt import Voigt

from .boron import B11, natural_boron

def test_doppler_broadening_plot():

    amu = 10.
    Delta = 1e-6
    effective_temperature = MaxwellBoltzmann.get_effective_temperature(Delta, amu, 1.)
    Gammas = [0.1*Delta, Delta]
    Gamma_labels = [r'$\Gamma = \Delta / 10$', r'$\Gamma = \Delta$']

    ground_state = GroundState("0^+_1", 0, 1)

    cross_sections_at_rest = []
    cross_sections = []
    cross_section_approximations = []

    for Gamma in Gammas:
        excited_state = State("1^+_1", 2, 1, 1., {"0^+_1": Gamma})

        cross_sections_at_rest.append(
            BreitWigner(
                ground_state,
                excited_state,
            )
        )

        cross_sections.append(
            Voigt(
                ground_state,
                excited_state,
                amu,
                effective_temperature,
            )
        )
        cross_section_approximations.append(
            Gauss(
                ground_state,
                excited_state,
                amu,
                effective_temperature,
            )
        )

    energies = cross_sections_at_rest[0].equidistant_energy_grid(0.99, 1000) - cross_sections_at_rest[0].intermediate_state.excitation_energy
    energies_over_Delta = energies/Delta
    cross_section_at_maximum = cross_sections_at_rest[0](0., input_is_absolute_energy=False)

    fig, ax = plt.subplots(len(Gammas),1, figsize=(5, len(Gammas)*2.2))
    plt.subplots_adjust(hspace=0.)
    for i in range(len(Gammas)):
        if i < len(Gammas)-1:
            pass
        ax[i].tick_params(labelsize=8)
        ax[i].set_ylabel(r"$\sigma (E) / \sigma_0$")
        ax[i].plot(energies_over_Delta, cross_sections_at_rest[i](energies, input_is_absolute_energy=False)/cross_section_at_maximum, color='black', label='$\sigma_a$')
        ax[i].plot(energies_over_Delta, cross_sections[i](energies, input_is_absolute_energy=False)/cross_section_at_maximum, color='royalblue', label=r'$\tilde{\sigma}_a^D$')
        ax[i].plot(energies_over_Delta, cross_section_approximations[i](energies, input_is_absolute_energy=False)/cross_section_at_maximum, '--', color='orange', label=r'$\sigma_a^D$')
        ax[i].text(0.1, 0.7, Gamma_labels[i], transform=ax[i].transAxes, fontsize=12)
        ax[i].legend()
    ax[-1].set_xlabel(r"$(E - E_r) / \Delta$")
    plt.savefig("doppler_broadening_plot.pdf")