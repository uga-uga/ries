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

r"""
At first order in quantum electrodynamics, the cross section
:math:`\sigma_\mathrm{Compton} \left( E \right)` for Compton scattering of a photon with an
initial energy :math:`E` off a free charged particle with mass :math:`m` at rest is given by the
Klein-Nishina formula :cite:`KleinNishina1929`:

.. math:: \sigma_\mathrm{Compton} \left( E \right) = \frac{\pi \alpha^2 \hbar^2}{m^2 c^2} \frac{1}{x^3} \left\{ \frac{2x \left[ 2 + x \left( 1 + x \right) \left( 8 + x \right)\right]}{\left( 1 + 2 x \right)^2} + \left[ \left( x - 2 \right) x - 2 \right] \mathrm{log} \left( 1 + 2x \right) \right\}.

In the equation above, the abbreviation

.. math:: x = \frac{E}{m c^2}

has been used.
The symbols :math:`\alpha`, :math:`\hbar`, :math:`m`, and :math:`c` denote the fine-structure
constant, the reduced Planck constant, the charged particle's rest mass, and the speed of light,
respectively.

The corresponding solid-angle differential cross section is given by :cite:`Depaola2003`:

.. math:: \frac{\mathrm{d} \sigma_\mathrm{Compton}}{\mathrm{d}\Omega} (E, \Omega) = \frac{1}{2}\frac{\alpha^2 \hbar^2}{m^2 c^2} \left( \frac{E^\prime}{E} \right)^2 \left[ \frac{E^\prime}{E} + \frac{E}{E^\prime} - 2 \sin \left( \theta \right)^2 \cos \left( \varphi \right)^2 \right].

Here, the angle :math:`\theta` denotes the polar angle of the scattered photon with respect to
its initial direction of propagation.
The angle :math:`\varphi` is the corresponding azimuthal angle, which is  defined with respect
to the direction of polarization of the initial photon.
The equation above also contains the energy of the scattered photon, :math:`E^\prime`, which is related to :math:`E` via:

.. math:: E^\prime \left( E, \theta \right) = \frac{E}{1 + \frac{E}{m c^2} \left[ 1 - \cos \left(\theta \right) \right]}.

If no information on the polarization of the incident photon is available, the differential
cross section can be averaged over the azimuthal angle:

.. math:: \frac{\mathrm{d} \sigma_\mathrm{Compton, unpolarized}}{\mathrm{d}\Omega} (E, \Omega) = \frac{\int_0^{2\pi} \frac{\mathrm{d} \sigma_\mathrm{Compton}}{\mathrm{d}\Omega} (E, \Omega) \mathrm{d} \varphi}{\int_0^{2\pi} \mathrm{d} \varphi} = \frac{\pi \alpha^2 \hbar^2}{m^2 c^2} \left( \frac{E^\prime}{E} \right)^2 \left[ \frac{E^\prime}{E} + \frac{E}{E^\prime} - \sin \left( \theta \right)^2 \right].

An integration of either :math:`\mathrm{d} \sigma_\mathrm{Compton} / \mathrm{d} \Omega` or
:math:`\mathrm{d} \sigma_\mathrm{Compton, unpolarized} / \mathrm{d} \Omega` over the azimuthal
angle gives the scattering-angle differential cross section:

.. math:: \frac{\mathrm{d} \sigma_\mathrm{Compton}}{\mathrm{d} \theta} (E, \theta) = \int_0^{2\pi} \frac{\mathrm{d} \sigma_\mathrm{Compton}}{\mathrm{d}\Omega} (E, \Omega) \sin \left( \theta \right) \mathrm{d} \varphi = \frac{\pi \alpha^2 \hbar^2}{m^2 c^2} \left( \frac{E^\prime}{E} \right)^2 \left[ \frac{E^\prime}{E} + \frac{E}{E^\prime} - \sin \left( \theta \right)^2 \right] \sin \left( \theta \right).

Note the additional factor :math:`\sin \left( \theta \right)`, which comes from the fact that:

.. math:: \mathrm{d} \Omega = \mathrm{d} \cos \left( \theta \right) \mathrm{d} \varphi =  \sin \left( \theta \right) \mathrm{d} \theta \mathrm{d} \varphi.

The differential cross section with respect to :math:`\cos \left( \theta \right)` can be
converted to the energy-differential cross section using the relation to :math:`E^\prime` above:

.. math:: \mathrm{d} \cos \left( \theta \right) = \frac{m c^2}{E^{\prime 2}} \mathrm{d} E^\prime.

Inserting the relation above in the expression for the scattering-angle differential equation
and abbreviating the inverse relation between :math:`\theta` and :math:`E^\prime` as
:math:`\theta \left( E, E^\prime \right)`, one obtains:

.. math:: \frac{\mathrm{d} \sigma_\mathrm{Compton}}{\mathrm{d} E^\prime} (E, E^\prime) = \frac{\pi \alpha^2 \hbar^2}{E^2} \left\{ \frac{E^\prime}{E} + \frac{E}{E^\prime} - \sin \left[ \theta \left( E, E^\prime \right) \right]^2 \right\}.
"""

import numpy as np
from scipy.constants import physical_constants

from ries.nonresonant.nonresonant import Nonresonant


class KleinNishina(Nonresonant):
    r"""(Differential) Cross section for Compton scattering of a photon off a free charged particle

    Note that the Klein-Nishina formula gives the cross section for the scattering off a single
    charged particle.
    The default particle of this class is the electron, with a charge of :math:`1e` and a mass energy
    equivalent of about :math:`0.511 \mathrm{MeV}`.

    In low-energy nuclear physics, one is often interested in the cross section per atom, to be able
    to compare this electromagnetic process to the cross sections for nuclear reactions.
    To obtain the cross section for scattering off an atom, the cross section needs to be multiplied
    by the number of electrons per atom.

    Attributes:

    - `Z`: int or float, charge in units of elementary charges (default: 1).
    - `mc2`: float, mass energy equivalent (mass times speed of light squared) of the charged particle in :math:`\mathrm{MeV}` (default: rest mass of an electron from `scipy.constants.physical_constants`, i.e.: :math:`m_e c^2 \approx 0.511 \mathrm{MeV}`).
    - `scale_factor`: float, the scale factor of the Klein-Nishina cross section that is independent of the kinematics, i.e. :math:`Z \alpha^2 \hbar^2 / \left( m^2 c^2 \right)`.
    """

    def __init__(
        self, Z=1, mc2=physical_constants["electron mass energy equivalent in MeV"][0]
    ):
        r"""Initialization

        Parameters:

        - `Z`: int or float, charge in units of elementary charges (default: 1).
        - `m`: float, mass of the charged particle in :math:`\mathrm{MeV} c^{-2}` (default: rest mass of an electron from `scipy.constants.physical_constants`, i.e.: :math:`m_e \approx 0.511 \mathrm{MeV}  c^{-2}`).
        """
        self.Z = Z
        self.mc2 = mc2
        self.scale_factor = (
            self.Z
            * physical_constants["fine-structure constant"][0] ** 2
            * physical_constants["Planck constant over 2 pi times c in MeV fm"][0] ** 2
            / (self.mc2 * self.mc2)
        )

    def __call__(self, E):
        r"""Total cross section

        This method has been implemented for convenience and calls `KleinNishina.cs_total()`.

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.

        Returns:

        - array_like or scalar, :math:`\sigma_\mathrm{Compton}` in :math:`\mathrm{fm}^2`.
        """
        return self.cs_total(E)

    def compton_edge(self, E):
        r"""Compton edge for a given initial photon energy

        The Compton edge is the lower limit for the energy :math:`E^\prime` of the photon after the
        scattering process.
        It follows from the conservation of energy and linear momentum.
        The limit corresponds to scattering by :math:`180^\circ`, which means maximum momentum transfer to the electron.
        The expression for the Compton edge is given by:

        .. math:: E^\prime \left( E, 180^\circ \right) = \frac{E}{1+2 \frac{E}{m c^2}}

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.

        Returns:

        - array_like or scalar, energy of the Compton edge in MeV.
        """

        return self.Ep_over_E(E, np.pi) * E

    @staticmethod
    def theta(E, Ep):
        """Scattering angle for a given scattered energy

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - Ep: array_like or scalar, energy of the photon after the scattering process in MeV.

        Returns

        - array_like or scalar, scattering angle of the photon in radians.
        """
        return np.arccos(
            1
            - (E / Ep - 1)
            * physical_constants["electron mass energy equivalent in MeV"][0]
            / E
        )

    @staticmethod
    def Ep_over_E(E, theta):
        """Ratio of final and initial photon energy for a given scattering angle

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - theta: array_like or scalar, scattering angle of the photon in radians.

        Returns:

        - array_like or scalar, ratio of the energy of the scattered photon and its initial energy.
        """
        return 1.0 / (
            1.0
            + E
            / physical_constants["electron mass energy equivalent in MeV"][0]
            * (1.0 - np.cos(theta))
        )

    def cs_diff_dOmega(self, E, theta, phi):
        """Differential cross section w.r.t. the solid-angle

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - theta: array_like or scalar, scattering angle of the photon in radians.
        - phi: array_like or scalar, azimuthal angle of the scattered photon with respect to the polarization vector of the incoming photon.

        Returns:

        - array_like or scalar, solid-angle differential cross section in fm**2.
        """
        relative_energy_change = self.Ep_over_E(E, theta)
        return (
            0.5
            * self.scale_factor
            * relative_energy_change
            * relative_energy_change
            * (
                relative_energy_change
                + 1.0 / relative_energy_change
                - 2.0 * np.sin(theta) ** 2 * np.cos(phi) ** 2
            )
        )

    def cs_diff_dOmega_unpolarized(self, E, theta):
        """Differential cross section w.r.t. the solid-angle for an unpolarized incoming photon

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - theta: array_like or scalar, scattering angle of the photon in radians.

        Returns:

        - array_like or scalar, solid-angle differential cross section in fm**2.
        """

        relative_energy_change = self.Ep_over_E(E, theta)
        return (
            np.pi
            * self.scale_factor
            * relative_energy_change
            * relative_energy_change
            * (
                relative_energy_change
                + 1.0 / relative_energy_change
                - np.sin(theta) ** 2
            )
        )

    def cs_diff_dEp_dphi(self, E, Ep, phi):
        r"""Differential cross section w.r.t. the energy of the scattered photon and the azimuthal angle

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - Ep: array_like or scalar, energy of the photon after the scattering process in MeV.
        - phi: array_like or scalar, azimuthal angle of the scattered photon with respect to the polarization vector of the incoming photon.

        Returns:

        - array_like or scalar, energy-differential cross section in :math:`\mathrm{fm}^2 \mathrm{MeV}^{-1}`.
        """
        return (
            self.cs_diff_dOmega(E, self.theta(E, Ep), phi)
            * physical_constants["electron mass energy equivalent in MeV"][0]
            / (Ep * Ep)
        )

    def cs_diff_dEp(self, E, Ep):
        r"""Differential cross section w.r.t. the energy of the scattered photon for an unpolarized incoming photon

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - Ep: array_like or scalar, energy of the photon after the scattering process in MeV.

        Returns:

        - array_like or scalar, energy-differential cross section in :math:`\mathrm{fm}^2 \mathrm{MeV}^{-1}`.
        """
        return (
            self.cs_diff_dOmega_unpolarized(E, self.theta(E, Ep))
            * physical_constants["electron mass energy equivalent in MeV"][0]
            / (Ep * Ep)
        )

    def cs_diff_dtheta(self, E, theta):
        r"""Differential cross section w.r.t. the polar scattering angle for an unpolarized incoming photon

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.
        - theta: array_like or scalar, scattering angle of the photon in radians.

        Returns:

        - array_like or scalar, scattering-angle differential cross section in :math:`\mathrm{fm}^2`.
        """

        return self.cs_diff_dOmega_unpolarized(E, theta) * np.sin(theta)

    def cs_total(self, E):
        r"""Total cross section

        Parameters:

        - E: array_like or scalar, initial energy of the photon in MeV.

        Returns:

        - array_like or scalar, total cross section in :math:`\mathrm{fm}^2`.
        """
        x = E / self.mc2
        return (
            np.pi
            * self.scale_factor
            * x ** -3
            * (
                (2 * x * (2 + x * (1 + x) * (8 + x))) / ((1 + 2 * x) ** 2)
                + ((x - 2) * x - 2) * np.log(1 + 2 * x)
            )
        )
