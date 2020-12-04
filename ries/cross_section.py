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
Classes for the description of a general energy-dependent reaction cross section.

The cross section :math:`\sigma_r` for a reaction of type :math:`r` in a beam-on-target
experiment *per target particle* is defined as [d]_:

.. math:: \sigma_r = \frac{\text{number of reactions of type }r}{\text{number of beam particles incident on the target} \times \text{number of target particles} \times \text{target area}}.

If the occurrence of a reaction of type :math:`r` precludes the occurrence of another reaction of type
:math:`r^\prime`, and vice versa, the cross section for any one [e]_ [f]_ of the two types is given by

.. math:: \sigma_{r \veebar r^\prime} = \sigma_r + \sigma_{r^\prime},

i.e. cross sections can be added.
Moreover, cross sections can be 'scaled', i.e. multiplied by a scalar constant :math:`c`.

This module implements classes for cross sections with the aforementioned algebraic properties.
The general cross sections are assumed to be a function of the kinetic energy :math:`E` of the
incident beam particle only.
Derived classes may introduce dependencies on other parameters.
The energy :math:`E` is restricted to the range

.. math::`0 \leq E`,

i.e. negative kinetic energies are unphysical.

*Note: At the moment, many of the classes do not perform type checks before executing the algebraic
operations.
This may lead to runtime `TypeError`s.*

.. [d] More precise definitions can be found in the literature.
    For example, Iliadis :cite:`Iliadis2015` makes the point that the target particles must not
    shadow each other.
    This is the case as long as the total cross section for a reaction of type :math:`r` on the
    target can be expressed as:

.. math:: \text{number of target particles} \times \sigma_r.

.. [e] Logical relation XOR (:math:`\veebar`): either :math:`r` or :math:`r^\prime` happens, but not both.

.. [f] For example: A photon with an energy :math:`E` may either interact via Compton scattering or
    the photoeffect.
    If it is Compton-scattered first, the photon's energy will have changed.
    If the photoeffect occurs, the photon does not exist any more.
    Therefore, the cross sections for Compton scattering and the photoeffect for a photon with an
    energy :math:`E` exclude each other.
"""

from itertools import chain

import numpy as np


class CrossSection:
    """Energy-dependent cross section for a reaction"""

    def __call__(self, E):
        """Evaluate the cross section for a given energy of the incident beam particle

        See `CrossSectionWeightedSum.__call__()`.
        """
        raise NotImplementedError

    def __add__(self, other):
        """Addition of cross sections

        See `CrossSectionWeightedSum.__add__()`.
        """
        if isinstance(other, (int, float)):
            return CrossSectionWeightedSum([self, ConstantCrossSection(other)])
        return CrossSectionWeightedSum([self, other])

    def __radd__(self, other):
        """Addition of cross sections

        See `self.__add__()`.
        """
        return self.__add__(other)

    def __mul__(self, other):
        """Multiplication of cross sections with a scalar

        This will create a `CrossSectionWeightedSum` object.

        See also `CrossSectionWeightedSum.__mul__()`.
        """
        return CrossSectionWeightedSum([self], [other])

    def __rmul__(self, other):
        """Multiplication of cross sections with a scalar

        See `self.__mul__()`.
        """
        return self.__mul__(other)

    def equidistant_energy_grid(self, limits, n_points):
        r"""Create an equidistant grid in a given 1D energy range

        Given a pair of limits :math:`E_0` and :math:`E_{n-1}` and a number of points :math:`n`, this
        function returns a set of :math:`n` energies :math:`\left\{ E_i \right\}` such that

        .. math :: E_{i + 1} - E_i = E_{j + 1} - E_j

        for all :math:`0 \leq i, j < n - 1`, i.e. the range is divided into partitions of equal size.

        This function simply wraps `numpy.linspace`.

        Parameters:

        - `limits`, pair of float, limits of the energy range.
        - `n_points`, int, number of energies that define the grid, i.e. number of partitions plus 1.

        Returns:

        - ndarray, array of grid points
        """
        return np.linspace(limits[0], limits[1], n_points)

    def equidistant_probability_grid(self, limits, n_points):
        r"""Create a grid with equal reaction probabilities per interval

        If the energy-integrated cross section :math:`I_r` of a reaction is finite, i.e.

        .. math:: \int_0^{\infty} \sigma_r \left( E \right) \mathrm{d} E \equiv I_r < \infty,

        the function

        .. math:: \frac{\sigma_r \left( E \right)}{I_r}

        may be interpreted as a probability density function (PDF).
        Given a pair of limits :math:`E_0` and :math:`E_{n-1}` and a number of points :math:`n`, this
        function returns a set of :math:`n` energies :math:`\left\{ E_i \right\}` such that

        .. math :: E_{i + 1} > E_i

        and

        .. math :: \int_{E_i}^{E_{i+1}} \sigma_r \left( E \right) \mathrm{d} E = \int_{E_j}^{E_{j+1}} \sigma_r \left( E \right) \mathrm{d} E

        for all :math:`0 \leq i, j < n - 1`, i.e. the energy range is partitioned into quantiles of the
        (truncated) PDF.

        Since many cross sections of physical processes exhibit sharp resonances and/or strong variations
        with the energy, a numerical evaluation on an equidistant-probability grid may be more efficient
        than an equidistant-energy grid.

        Parameters:

        - `limits`, pair of float, limits of the energy range.
        - `n_points`, int, number of energies that define the grid, i.e. number of partitions plus 1.

        Returns:

        - ndarray, array of grid points
        """
        raise NotImplementedError


class ConstantCrossSection(CrossSection):
    """Class for a constant cross section

    'Constant' means that the cross section is the same for all energies.
    In particular, this means that the energy-integrated cross section is not finite.

    Attributes:

    - `constant`, scalar, constant value of the cross section.
    """

    def __init__(self, constant):
        self.constant = constant

    def __call__(self, energy):
        """Evaluate the cross section for a given energy of the incident beam particle

        See also `CrossSection.__call__()`.
        """
        return self.constant

    def equidistant_probability_grid(self, limits, n_points):
        """Create an equidistant-probability grid in a given 1D energy range

        See also `CrossSection.equidistant_probability_grid()`.
        """
        return np.linspace(*limits, n_points)


class CrossSectionWeightedSum:
    """Weighted sum of multiple energy-dependent cross sections

    Internally, a list of `CrossSection` objects and a list of scale factors with the same length
    are stored.

    Attributes:

    - `reactions`, list of `CrossSection` objects, all cross sections that contribute to the weighted sum.
    - `scale_factors`, list of float. The i-th scale factor in the list is used to scale the i-th cross section in `self.reactions`.
    """

    def __init__(self, reactions=None, scale_factors=None):
        """Initialization

        Parameters:

        - `reactions`, list of `CrossSection` objects, all cross sections that contribute to the weighted sum.
        - `scale_factors`, list of float. The i-th scale factor in the list is used to scale the i-th cross section in `self.reactions`.
        """
        self.reactions = reactions or []
        if scale_factors is None:
            scale_factors = [1.0] * len(reactions)
        self.scale_factors = scale_factors

    def __add__(self, other):
        """Addition of cross sections

        An addition of two `CrossSectionWeightedSum` objects returns a new `CrossSectionWeightedSum`
        object which contains the union of the `reactions` and `scale_factors` of `self` and `other`.
        An addition of a `CrossSectionWeightedSum` object and a `CrossSection` object appends that
        `CrossSection` object to a new `CrossSectionWeightedSum` object's `reactions` and
        a float value of 1 to that one's `scale_factors`.
        An addition of a `CrossSectionWeightedSum` object and a scalar will initialize a
        `ConstantCrossSection` object with the scalar and add that one like a `CrossSection` object.

        It is emphasized again that each of the cases above returns a new `CrossSectionWeightedSum` object.

        Parameters:

        - `other`, scalar, `CrossSection` or `CrossSectionWeightedSum` object.
          If a scalar is given, it is assumed to be an energy-independent cross section.

        Returns:

        - `CrossSectionWeightedSum` object
        """
        reactions = []
        scale_factors = []
        for i, reaction in enumerate(self.reactions):
            reactions.append(reaction)
            scale_factors.append(self.scale_factors[i])

        if isinstance(other, (int, float)):
            reactions.append(ConstantCrossSection(other))
            scale_factors.append(1.0)
        elif isinstance(other, CrossSection):
            reactions.append(other)
            scale_factors.append(1.0)
        else:
            for i, reaction in enumerate(other.reactions):
                reactions.append(reaction)
                scale_factors.append(other.scale_factors[i])

        return CrossSectionWeightedSum(reactions, scale_factors)

    def __radd__(self, other):
        """Addition of cross sections

        See `self.__add__()`.
        """
        return self.__add__(other)

    def __mul__(self, other):
        """Multiplication of cross sections with a scalar

        Returns a new `CrossSectionWeightedSum` object in which all `scale_factors` are multiplied by
        `other`.

        Parameters:

        - `other`, scalar, scale factor.

        Returns:

        - `CrossSectionWeightedSum` object
        """
        reactions = []
        scale_factors = []
        for i, reaction in enumerate(self.reactions):
            reactions.append(reaction)
            scale_factors.append(other * self.scale_factors[i])
        return CrossSectionWeightedSum(reactions, scale_factors)

    def __rmul__(self, other):
        """Multiplication of cross sections with a scalar

        See `self.__mul__()`
        """
        return self.__mul__(other)

    def __call__(self, energy):
        """Evaluate the cross section for a given energy of the incident beam particle

        Invokes the `__call__` method of every `CrossSection` object in `self.reactions` and returns the
        `self.scale_factors`-weighted sum.

        Parameters:

        - `E`, float or array_like, energy of the incident beam particle.

        Returns:

        - float or array_like, value of the cross section at the given energy.
        """
        cross_section = 0.0

        for i, reaction in enumerate(self.reactions):
            cross_section += self.scale_factors[i] * reaction(energy)

        return cross_section

    def equidistant_energy_grid(self, limits, n_points):
        """Create an equidistant grid in a given 1D energy range

        The equidistant energy grid is the same for any cross section in `self.reactions`.

        See also `CrossSection.equidistant_energy_grid()`.

        Parameters:

        - `limits`, pair of float, limits of the energy range.
        - `n_points`, int, number of energies that define the grid, i.e. number of partitions plus 1.

        Returns:

        - ndarray, array of grid points
        """
        return np.linspace(limits[0], limits[1], n_points)

    def equidistant_probability_grid(self, limits, n_points_per_cross_section):
        """Return a union of the equal-reaction-probability grids of all constituent cross sections

        Returns a union of the equidistant-probability grids returned by all individual `CrossSection`
        objects in `self.reactions` in which every energy occurs only once.
        Note that all constituents are treated as if they had the same energy-integrated cross
        section.
        This behavior was chosen on purpose so that narrow resonances are not neglected when compared to
        broad nonresonant processes.

        See also `CrossSection.equidistant_probability_grid()`.
        """
        grid = []
        for reaction in self.reactions:
            grid.append(
                reaction.equidistant_probability_grid(
                    limits, n_points_per_cross_section
                )
            )
        grid = list(chain(*grid))
        return np.unique(grid)
