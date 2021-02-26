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

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

from ries.constituents.state import GroundState, State
from ries.resonance.voigt import Voigt


class BeamInTarget:
    def __init__(self, sigma, kappa):
        self.sigma = sigma
        self.kappa = kappa

    def photon_flux_density(self, Z, E):
        return np.exp(-(self.kappa(E) + self.sigma(E)) * Z)

    def resonance_absorption_density(self, Z, E):
        return self.sigma(E) * self.photon_flux_density(Z, E)


def test_doppler_broadening_plot():

    sigma = Voigt(
        GroundState("0^+_1", 0, 1),
        State("1^+_1", 2, 1, 1.0, {"0^+_1": 1e-7}),
        10.0,
        100.0,
    )

    Delta = sigma.probability_distribution.doppler_width
    cross_section_at_maximum = sigma(0.0, input_is_absolute_energy=False)

    K = cross_section_at_maximum
    K_over_kappa = 5
    kappa = lambda energy: K / K_over_kappa

    beam_in_target = BeamInTarget(sigma, kappa)

    e = sigma.equidistant_energy_grid(0.995, 100)
    z = np.linspace(0.0, 2.0 / K, 50)

    Z, E = np.meshgrid(z, e)
    Phi = beam_in_target.photon_flux_density(Z, E)
    alpha = beam_in_target.resonance_absorption_density(Z, E)

    _ccount = 10
    _cmap = "rainbow"
    _figsize = (5.5, 5.0)
    _fontsize_text = 14
    _rcount = 10
    _view = (40, -35)
    _wireframe_color = "grey"

    _xlabel = r"$\mathcal{Z} \times K$"
    _xlim = [0.0, z[-1] * K]
    _xticks = [0.0, 0.5, 1.0, 1.5, 2.0]
    _ylabel = r"$(E - E_r) / \Delta$"
    _K_over_kappa_label = r"$\kappa = K/{:d}$".format(K_over_kappa)
    _zlim = [0.0, 1.0]

    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel(_xlabel)
    ax.set_xlim(_xlim)
    ax.set_xticks(_xticks)
    ax.set_ylabel(_ylabel)
    ax.set_zlabel(r"$\Phi_K (\mathcal{Z}, E)$")
    ax.set_zlim(_zlim)
    ax.plot_surface(
        Z * K,
        (E - sigma.probability_distribution.resonance_energy) / Delta,
        Phi,
        cmap=_cmap,
    )
    ax.plot_wireframe(
        Z * K,
        (E - sigma.probability_distribution.resonance_energy) / Delta,
        Phi,
        color=_wireframe_color,
        rcount=_rcount,
        ccount=_ccount,
    )
    ax.view_init(*_view)
    ax.text(1.8, -3.2, 0.0, _K_over_kappa_label, (0, 1, 0), fontsize=_fontsize_text)

    plt.savefig("photon_flux_density.pdf")

    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlabel(_xlabel)
    ax.set_xlim(_xlim)
    ax.set_xticks(_xticks)
    ax.set_ylabel(_ylabel)
    ax.set_zlabel(r"$\alpha_K (\mathcal{Z}, E) / K$")
    ax.set_zlim(_zlim)
    ax.plot_surface(
        Z * K,
        (E - sigma.probability_distribution.resonance_energy) / Delta,
        alpha / K,
        cmap=_cmap,
    )
    ax.plot_wireframe(
        Z * K,
        (E - sigma.probability_distribution.resonance_energy) / Delta,
        alpha / K,
        color=_wireframe_color,
        rcount=_rcount,
        ccount=_ccount,
    )
    ax.view_init(*_view)
    ax.text(0.0, -3.5, 0.73, _K_over_kappa_label, (0, 1, 0), fontsize=_fontsize_text)

    plt.savefig("resonance_absorption_density.pdf")
