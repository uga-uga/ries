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

r"""Doppler-broadened Breit-Wigner (Voigt profile) approximated using a pseudo-Voigt distribution.

Whether for simplicity of the implementation or for computational speed, the Voigt distribution 
:math:`P^\mathrm{Voigt}` is often approximated by a pseudo-Voigt distribution 
:math:`P^\mathrm{PV}` :cite:`Ida2000`:

.. math:: P^\mathrm{Voigt} \left( E, \left\{ E_r, \sigma, \Gamma \right\} \right) = \int_{-\infty}^\infty P^\mathrm{normal} \left( E^\prime - E, \left\{ 0, \sigma \right\} \right) P^\mathrm{Cauchy} \left( E^\prime, \left\{ E_r, \Gamma \right\} \right) \mathrm{d} E^\prime

.. math:: \approx P^\mathrm{PV} \left( E, \left\{ E_r, \sigma, \Gamma \right\} \right) = \left[ 1 - \eta \left( \sigma, \Gamma \right) \right] P^\mathrm{normal} \left( E, \left\{ 0, \sigma \right\} \right) + \eta \left( \sigma, \Gamma \right) P^\mathrm{Cauchy} \left( E, \left\{ E_r, \Gamma \right\} \right)

The pseudo-Voigt approximation is a linear combination of a normal distribution and a Cauchy 
distribution instead of a convolution.
In the limiting cases :math:`\sigma \gg \Gamma` and :math:`\Gamma \gg \sigma`, the approximation
is exact, and the largest deviations are expected if their values are on the same order of 
magnitude.
The expression for the pseudo-Voigt distribution above contains a mixing parameter :math:`\eta` 
with an arbitrary dependence on :math:`\sigma` and `\Gamma`.

This implementation uses a comparably simple proposal for 
:math:`\eta \left( \sigma, \Gamma \right)` by Thompson *et al.* :cite:`Thompson1987`, which was 
found to deviate from the exact result by less than 1.2 % even when the width and the scale
parameter are of the same order of magnitude :cite:`Ida2000`.
The authors of the latter publication propose an even better approximation which has 42 parameters 
(Tab. 1 therein) instead of 7 in Thompson *et al.*.

From the PDF of the pseudo-Voigt distribution, its CDF can be obtained as:

.. math:: F^\mathrm{PV} \left( E \right) = \int_{-\infty}^E P^\mathrm{PV} \left( E^\prime \right) \mathrm{d} E^\prime

.. math:: = \int_{-\infty}^E \left[ 1 - \eta \left( \sigma, \Gamma \right) \right] P^\mathrm{normal} \left( E^\prime, \left\{ 0, \sigma \right\} \right) + \eta \left( \sigma, \Gamma \right) P^\mathrm{Cauchy} \left( E^\prime, \left\{ E_r, \Gamma \right\} \right) \mathrm{d} E^\prime

.. math:: = \left[ 1 - \eta \left( \sigma, \Gamma \right) \right] F^\mathrm{normal} \left( E, \left\{ 0, \sigma \right\} \right) + \eta \left( \sigma, \Gamma \right) F^\mathrm{Cauchy} \left( E, \left\{ E_r, \Gamma \right\} \right).

Since the integral with respect to :math:`E` is a linear funtional, the CDF is again a linear
combination of the CDFs of the normal distribution and the Cauchy distribution.

The PPF :math:`\left( F^\mathrm{PV} \right)^{-1}` of the pseudo-Voigt distribution, on the other 
hand, is not a simple linear combination of the PPFs of the normal distribution and the Cauchy 
distribution.
It is defined as:

.. math:: \left( F^\mathrm{PV} \right)^{-1} \left[ F^\mathrm{PV} \left( E \right) \right] = E.

Since the inversion of :math:`F` is nontrivial, apply :math:`F` to both sides of this equation
to obtain a condition that only contains known functions:

.. math:: g \left( E \right) \equiv F^\mathrm{PV} \left( E \right) - Q \overset{!}{=} 0.

Here, the quantile :math:`Q` that corresponds to :math:`E` has been introduced.
Since the CDFs of the normal- and Cauchy distributions are strictly increasing, continuous 
functions there is a unique solution :math:`E` for a given :math:`Q`.
In the present code, this root of the function :math:`g \left( E \right)` is found numerically
using the Newton-Raphson method of `scipy` (`scipy.optimize.newton`).
The algorithm requires knowledge of the first derivative of :math:`g`, which is :math`P`.

When experimenting with the numerical inversion of the CDF above, it was found that the algorithm 
produces a lot of `RuntimeWarnings` due to nan expressions and overflows.
The problem is that for extreme values of the energy far away from the resonance energy, the PDF
and the CDF are 'almost zero' in terms of numerical precision.
This causes the root-finding algorithm to fail in reaching the default precision goal.
At the moment, a warning is issued when this happens.
To avoid continuing with the nonsense values that are produced by overflows, the following 
expression is used as a fallback solution:

.. math:: \left( F^\mathrm{PV} \right)^{-1} \left( Q \right) \approx \left[ 1 - \eta \left( \sigma, \Gamma \right) \right] \left( F^\mathrm{normal} \right)^{-1} \left( Q, \left\{ E_r, \sigma \right\} \right) + \eta \left( \sigma, \Gamma \right) \left( F^\mathrm{Cauchy} \right)^{-1} \left( Q, \left\{ E_r, \Gamma \right\} \right).

It assumes that the PPF is a weighted mean of the PPFs of the normal- and the Cauchy distribution.

This module implements both a pseudo-Voigt distribution as required by 
`ries.resonance.resonance.Resonance` as well as a resonance with a pseudo-Voigt peak.
`scipy.special` provides at least an 'exact' implementation of the PDF, which replaces the
pseudo-Voigt PDF in `ries.resonance.voigt`.
"""

import warnings

import numpy as np
from scipy.optimize import newton
from scipy.stats import cauchy, norm

from ries.resonance.maxwell_boltzmann import MaxwellBoltzmann
from ries.resonance.resonance import Resonance


class PseudoVoigtDistribution:
    """Class for a pseudo-Voigt distribution

    Approximates a Voigt distribution by a linear combination of a normal distribution and a
    Cauchy distribution instead of a convolution.
    For the mixing constant, an expression by Thompson *et al.* :cite:`Thompson1987` is used.

    Attributes:

    - `resonance_energy`, float, resonance energy, i.e. location of the centroid in MeV.
    - `maxwell_boltzmann`, `MaxwellBoltzmann` object, used to calculate the doppler width.
    - `Gamma_L_*`, `Gamma_G_*`, float, parameters used in the approximation of Thompson et al. in units of different powers of MeV.
    - `sigma`, float, standard deviation of the Doppler-broadened normal distribution in MeV.
    - `Gamma`, float, full width at half maximum (FWHM) of the normal distribution in MeV, used in the approximation of Thompson et al..
    - `gamma`, float, scale parameter of the Cauchy distribution.
    - `eta`, float, mixing parameter that controls the relative contributions of the normal- and the Cauchy distribution to the linear combination.
    """

    def __init__(self, resonance_energy, width, amu, effective_temperature):
        """Initialization

        Parameters:

        - `resonance_energy`, float, resonance energy, i.e. location of the centroid in MeV.
        - `width`, float, the width of the excited state in MeV.
        - `amu`, float, mass of the nucleus in atomic mass units.
        - `effective_temperature`, float, effective temperature of the ensemble of nuclei in K.
        """
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

    def pseudo_voigt_expression(self, x, method):
        """Evaluate linear combination of normal- and Cauchy distribution

        To avoid repetitive code, this function replaces 'pdf', 'cdf', and 'ppf' in the general
        linear-combination expression.

        Parameters:

        - `x`, float or array_like, meaning depends on the `method` parameter, but it is most probably an energy in MeV (`method == 'pdf'` or `method == 'cdf'`)  or a quantile (`method == ppf`).
        - `method`, str, method of the distributions that should be called, for example 'pdf', 'cdf', or 'ppf'.

        Returns:

        - float or array_like, meaning depends on the `method` parameter, but it is most probably an energy in MeV (`method == 'ppf'`), a probability (`method == pdf`), or a quantile (`method == cdf`).
        """
        return eval(
            """(1.0 - self.eta) * norm.{}(
            x, loc=self.resonance_energy, scale=self.gamma_G / np.sqrt(2.0)
        ) + self.eta * cauchy.{}(
            x, loc=self.resonance_energy, scale=0.5 * self.gamma_L
        )""".format(
                method, method
            )
        )

    def cdf(self, E):
        """CDF of the pseudo-Voigt distribution

        Parameters:

        - `E`, float or array_like, energy in MeV.

        Returns:

        float or array_like, CDF
        """
        return self.pseudo_voigt_expression(E, "cdf")

    def pdf(self, E):
        """PDF of the pseudo-Voigt distribution

        Parameters:

        - `E`, float or array_like, energy in MeV.

        Returns:

        float or array_like, PDF
        """
        return self.pseudo_voigt_expression(E, "pdf")

    def ppf(self, quantile):
        """PPF of the pseudo-Voigt distribution

        This function tries to invert the pseudo-Voigt CDF numerically first.
        If the root-finding algorithm does not converge for any of the given quantiles, the 
        RuntimeError issued by the algorithm is caught and a fallback approximation is used which 
        is based on a linear combination of PPFs.

        Parameters:

        - `quantile`, float or array_like, quantile.

        Returns:

        float or array_like, PPF
        """

        try:
            # scipy.optimize.newton returns different output depending on whether quantile is a 
            # scalar or an array.
            # For a scalar, the root is `newton_result[0]` and the flag that indicates whether 
            # the algorithm converged is `newton_result[1]`.converged`.
            # In the case of an array, it is `newton_result[0]` and `newton_result[1]`.
            newton_result = newton(
            lambda E: self.pseudo_voigt_expression(E, "cdf") - quantile,
            # Use resonance energy as start value.
            # There are closer guesses for an arbitrary q, but starting at the resonance 
            # energy ensures that the algorithm does not have rounding errors from the start 
            # which might cause it to go in the wrong direction.
            self.resonance_energy if isinstance(quantile, (float, int)) else self.resonance_energy*np.ones(len(quantile)),
            fprime=lambda E: self.pseudo_voigt_expression(E, "pdf"),
            full_output=True) 
            if isinstance(quantile, (int, float)):
                if newton_result[1].converged:
                    return newton_result[0]
            elif False not in newton_result[1]:
                return newton_result[0]

        except RuntimeError:
            warnings.warn(
                "Calculation of pseudo-Voigt PPF by numerical inversion of the CDF failed for at least one value. Using approximation instead.",
                UserWarning,
            )
            return self.pseudo_voigt_expression(quantile, "ppf")

    def get_Gamma(self):
        """Calculate Gamma parameter

        This parameter is used internally in the approximation of Thompson et al..
        """
        return (
            self.Gamma_G_squared_squared * self.Gamma_G
            + 2.69269 * self.Gamma_G_squared_squared * self.Gamma_L
            + 2.42843 * self.Gamma_G_squared * self.Gamma_G * self.Gamma_L_squared
            + 4.47163 * self.Gamma_G_squared * self.Gamma_L_cubed
            + 0.07842 * self.Gamma_G * self.Gamma_L_squared_squared
            + self.Gamma_L_squared_squared * self.Gamma_L
        ) ** 0.2

    def get_eta(self):
        """Calculate mixing parameter"""
        inverse_Gamma = 1.0 / self.Gamma
        inverse_Gamma_squared = inverse_Gamma * inverse_Gamma
        return (
            1.36603 * self.Gamma_L * inverse_Gamma
            - 0.47719 * self.Gamma_L_squared * inverse_Gamma_squared
            + 0.11116 * self.Gamma_L_cubed * inverse_Gamma_squared * inverse_Gamma
        )


class PseudoVoigt(Resonance):
    r"""Approximation for a Doppler-broadened Breit-Wigner cross section (pseudo-Voigt profile)

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

        self.probability_distribution = PseudoVoigtDistribution(
            self.resonance_energy,
            self.intermediate_state.width,
            amu,
            effective_temperature,
        )
        self.probability_distribution_parameters = ()
