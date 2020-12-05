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
Breit-Wigner expression for the cross section of an isolated resonance.

The Breit-Wigner expression :cite:`BreitWigner1936` [m]_ can be derived without making any assumptions 
about the form of the interaction, therefore it is generally valid for resonant reactions 
involving a single isolated state of the nucleus [n]_.
In the case that an incident particle with a spin [o]_ :math:`s` excites a nucleus from a state 
:math:`J_0` to a state :math:`J_2`, the cross section is given by:

.. math:: \sigma_{0 \to 2} \left( E \right) = \pi \underbrace{\left( \frac{\hbar c}{E_r} \right)^2}_{\left( \frac{\Lambda_{2 \to 0}}{2 \pi} \right)^2} \underbrace{\frac{2 J_2 + 1}{\left( 2 s + 1 \right)\left( 2 J_0 + 1\right)}}_{g} \frac{\Gamma_{2 \to 0} \Gamma_2}{\left( E - E_r \right)^2 + \left(\frac{\Gamma_2}{2}\right)^2}.

In this expression, :math:`E_r` denotes the resonance energy, which is equal to the level energy 
difference :math:`E_2 - E_0` up to a recoil correction.
The brace on the left indicates its relation with the de-Broglie wave length 
:math:`\Lambda_{2 \to 0}` of the incident particle.
Note that the approximation of a narrow resonance has already been introduced here 
(see also `ries.resonance.resonance`).
The quantity g is a statistical factor that takes into account all possible magnetic substates of 
the involved angular momenta.
A photon has an intrinsic 'spin' (helicity) of 1, but due to relativistic effects it can only \
have two configurations.
Therefore, its contributions to the statistical factor is:

.. math:: 2s + 1 = 2.

like the one of a spin-:math:`1/2` particle.
At the moment, `ries` uses the photon-spin factor by default.

The symbol :math:`\Gamma_{2 \to 0}` denotes the partial width for the decay from state 2 to 
state 0, and :math:`\Gamma_2` denotes the total width of the state.

The energy dependence of the Breit-Wigner cross section is the one of a Cauchy distribution 

.. math:: P_\mathrm{Cauchy} \left( E, E_r, \Gamma_2 / 2 \right)

with the location :math:`E_r` and the scale `\Gamma_2 / 2`, i.e. it can be written as:

.. math:: \sigma_r \left( E \right) = I_{0 \to 2} P_\mathrm{Cauchy} \left( E, E_r, \Gamma_2 / 2 \right),

with the integrated cross section :math:`I_{0 \to 2}`,

Please note that the domain of the Cauchy distribution is the entire set of real numbers, not only 
the nonnegative real numbers as expected for a physical cross section.
This is because the Breit-Wigner formula is only an approximation for the cross section in the 
vicinity of a resonance.
For this reason, `ries.resonance.resonance.equidistant_probability_grid()` may not be able to
meet a request for a coverage interval close to 1.
If a coverage interval contains a negative number, the function will issue a warning and 
artificially set the lower limit to zero.
The same problem also holds for other cross sections that are based on the Breit-Wigner expression.

.. [m] See also more recent textbooks that deal with scattering theory, for example Refs. :cite:`BlattWeisskopf1979` or :cite:`Iliadis2015`.

.. [n] Here, 'isolated' means that the distance to the neighboring states is much larger than the total width :math:`\Gamma_i` of the state.

.. [o] The notation 'spin' here is just for brevity. The incident particle may also have a nonzero orbital angular momentum w.r.t. the target nucleus.
"""

from scipy.stats import cauchy

from ries.resonance.resonance import Resonance


class BreitWigner(Resonance):
    r"""Class for a Breit-Wigner cross section.

    See `ries.resonance.resonance.Resonance`.
    """

    def __init__(self, initial_state, intermediate_state, final_state=None):
        r"""Initialization

        See `ries.resonance.resonance.Resonance`.
        """
        Resonance.__init__(self, initial_state, intermediate_state, final_state)

        self.probability_distribution = cauchy
        self.probability_distribution_parameters = (
            self.resonance_energy,
            0.5 * self.intermediate_state.width,
        )
