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

r"""Doppler-broadened Breit-Wigner cross section assuming Maxwell-Boltzmann velocity distribution (Voigt profile)

The Breit-Wigner cross section is valid only if the target nucleus is at rest.
If it is in motion with a velocity component :math:`v_\parallel` in the direction of 
propagation of the incident photon, the photon's energy :math:`E` will be Doppler shifted to 
:math:`E^\prime` in the rest frame of the nucleus.
In the nonrelativistic approximation, which is appropriate for even the highest 
mechanically attainable velocities of macroscopic samples [p]_:

.. math:: E^\prime \approx (1 + \frac{v_\parallel}{c}) E

Here, :math:`c` denotes the speed of light.
The velocity is defined such that a positive value of :math:`v_\parallel` indicates that the 
nucleus moves in the same direction as the photon (which decreases the photon's energy in the rest 
frame of the nucleus).

In a statistical ensemble, the velocity projections of the atoms follow a probability distribution:

.. math:: P_v \left( v_\parallel \right) \mathrm{d} v_\parallel = \overbrace{P_v \left[ \frac{c}{E} \left( E^\prime - E \right) \right] \underbrace{\frac{c}{E}}_{\frac{\mathrm{d} v_\parallel}{\mathrm{d} E^\prime}}}^{P_E \left( E^\prime - E \right) } \mathrm{d} E^\prime.

In the second part of the equation above, the velocity distribution was converted to a 
distribution of Doppler-shifted energies by performing an integration by substitution 'backwards'.

The Doppler-broadened cross section :math:`\sigma_r^\prime` at an energy :math:`E` is obtained 
from a weighted averge of all possible Doppler-shifts [q]_:

.. math:: \sigma_r^\prime \left( E \right) = \int_{-\infty}^\infty P_E \left( E^\prime - E \right) \sigma_r \left( E^\prime \right) \mathrm{d} E^\prime = I_r \int_{-\infty}^\infty P_E \left( E^\prime - E \right) P_r \left( E^\prime \right) \mathrm{d} E^\prime.

From a computational point of view, the fact that the integral in the equation above is a 
convolution of two function is advantageous, because this mathematical operation can be performed
in an efficient way by using a fast Fourier transform (FFT).

In the present module, it is assumed that :math:`P_v` is a 1D Maxwell-Boltzmann distribution:

.. math:: P_v \left( v_\parallel \right) = \sqrt{\frac{m \left( ^A\mathrm{X}\right)}{2 \pi k_B T}} \exp \left[ - \frac{m \left( ^A\mathrm{X}\right) v_\parallel^2}{2 k_B T} \right].

Here, :math:`m \left( ^A\mathrm{X}\right)` is the mass of 
the nucleus [r]_, 
:math:`k_B` is the Boltzmann constant, and :math:`T` is the temperature of the ensemble.
The assumption of a Maxwell-Boltzmann distribution is obviously valid for gaseous materials that \
are similar to an ideal gas, but calculations in momentum :cite:`Lamb1939` and real 
:cite:`VanHove1954` (see also :cite:`SingwiSjoelander1960`) space show that it is also a good 
approximation for harmonic crystals when the temperature is replaced by an effective temperature 
:math:`T_\mathrm{eff}` that takes into account the impact of the atomic lattice on the motion of 
the atom.

The corresponding distribution of the Doppler shifts is:

.. math:: P_E \left( E^\prime - E \right) = \sqrt{\frac{m \left( ^A\mathrm{X}\right) c^2 }{2 \pi k_\mathrm{B} T E^2}} \exp \left[ - \frac{m \left( ^A\mathrm{X}\right) c^2 \left( E^\prime - E \right)^2}{2 E^2 k_B T} \right] \equiv \frac{1}{\sqrt{\pi} \Delta} \exp \left[ - \left( \frac{E^\prime - E}{\Delta} \right)^2 \right].

In the last step, the 'Doppler width'

.. math:: \Delta = \sqrt{2 k_B T}{ \left( ^A\mathrm{X}\right) c^2 }

was introduced.
The distribution :math:`P_E` is a normal distribution with the mean value :math:`E` and the 
standard deviation :math:`\Delta / \sqrt{2}`.
Consequently, the Doppler-broadened Breit-Wigner cross section is the convolution of a normal
distribution and a cauchy distribution, which results in a so-called Voigt distribution:

.. math:: \sigma^\mathrm{Voigt}_r \left( E \right) = I_r P_E^\mathrm{Voigt} \left( E, \left\{ E_r, \Delta / \sqrt{2}, \Gamma_2 / 2 \right\} \right) = I_r \int_{-\infty}^\infty P_E^{normal} \left( E^\prime - E, \left\{ 0, \Delta / \sqrt{2} \right\} \right) P_E^{Cauchy} \left( E^\prime, \left\{ E_r, \Gamma_2 / 2 \right\} \right) \mathrm{d} E^\prime

The Voigt distribution

.. math:: P_E^\mathrm{Voigt} \left( E, \left\{ \mu, \sigma, \Gamma \right\} \right),

which uses the parameters of the normal- and the Cauchy distribution, interpolates between 
both and approaches them asymptotically as :math:`\sigma / \Gamma` goes to infinity or zero, 
respectively.

There is no closed analytic expression for the Voigt distribution, and, to the knowledge of the 
author, numerical libraries like `scipy` (which is used in the `ries` code), implement at most its 
PDF.
For the CDF and PPF, a pseudo-Voigt distribution is used here, which is a linear combination of a 
normal distribution and a Cauchy distribution instead of a convolution 
(see, e.g., :cite:`Ida2000`).
Since the CDF and the PPF are only used for the determination of an equal-probability grid, the 
relative deviations on the order of few percent :cite:`Ida2000` from the true Voigt distribution
are negligible.

.. [p] For example, modern centrifuges can reach rotational frequencies on the order of :math:`10^4 \mathrm{Hz}` (see, e.g., :cite:`ArabgolSleator2019`, where the outer edge of the sample reached a velocity of about :math:`13.5 \mathrm{kHz} \times 0.004 \mathrm{m} = 54 \mathrm{ms}^{-1}`).
.. [q] Note that the integral requires unphysical negative energies as an argument for :math:`\sigma_r`. This is due to the nonrelativistic approximation for the Doppler shift. In the present case, however, the Breit-Wigner cross section, which will be substituted for :math:`sigma_r`, is also defined on the entire set of real numbers due to an approximation (see `ries.resonance.breit_wigner`). At moderate velocities and for narrow resonances, i.e. when both approximations are applicable, the negative-energy terms should be negligible.
.. [r] At normal conditions, the particles in motion would be atoms instead of bare atomic nuclei, but the contribution of the electrons' masses and their binding energies were neglected here.
"""

from scipy.special import voigt_profile

from ries.resonance.pseudo_voigt import PseudoVoigt, PseudoVoigtDistribution


class SemiPseudoVoigtDistribution(PseudoVoigtDistribution):
    """Class for a Voigt distribution, using pseudo-Voigt expressions for CDF and PPF

    The `pdf()` method calls `scipy.special.voigt_profile`, while the CDF and the PPF are approximated
    by a pseudo-Voigt distribution (see `ries.resonance.pseudo_voigt`).

    See also `ries.resonance.pseudo_voigt.PseudoVoigtDistribution`.
    """

    def pdf(self, E):
        """PDF of Voigt distribution

        Wraps `scipy.special.voigt_profile`.

        Parameter:

        - `E`, float or array_like, energy of the incident beam particle in MeV.

        Returns:

        - float or array_like, PDF
        """
        return voigt_profile(E - self.resonance_energy, self.sigma, self.gamma)


class Voigt(PseudoVoigt):
    r"""Class for a Doppler-broadened Breit-Wigner cross section (Voigt profile)

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
        PseudoVoigt.__init__(
            self,
            initial_state,
            intermediate_state,
            amu,
            effective_temperature,
            final_state,
        )

        self.probability_distribution = SemiPseudoVoigtDistribution(
            self.resonance_energy,
            self.intermediate_state.width,
            amu,
            effective_temperature,
        )
