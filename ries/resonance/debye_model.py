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

r"""Debye model

The Debye model (see, e.g., Ref. :cite:`LandauLifshitz1980`) is an approximation for the specific \
heat of a harmonic solid, which assumes that the frequency distribution 
:math:`P^\mathrm{Debye} \left( \omega \right)` of the oscillators is 
given by:

.. math:: P^\mathrm{Debye} \left( \omega \right) = \frac{9 N k_B \omega^2}{\omega_\mathrm{D}} ~~~ \left( \omega \leq \omega_\mathrm{D}\right)

Here, :math:`N` is the number of particles in the solid, :math:`k_B` is the Boltzmann constant, 
and :math:`\omega_\mathrm{D}` is a cutoff frequency that is specific for the material of interest.
The Debye model interpolates between the :math:`T^3` law at low temperatures to the Dulong-Petit 
law at high temperatures.

The quantity

.. math:: T_\mathrm{D} = \hbar \omega_\mathrm{D}

is called the Debye temperature.
In a solid for which the Debye model is a good approximation, it can be shown :cite:`Lamb1939` 
that the particles move with a Maxwell-Boltzmann velocity distribution at an effective temperature

.. math:: T_\mathrm{eff} = 3 T \left( \frac{T}{T_\mathrm{D}} \right) \int_0^{\frac{T_\mathrm{D}}{T}} t^3 \left( \frac{1}{e^t - 1} + \frac{1}{2} \right) \mathrm{d}t

which is larger than the thermodynamic temperature :math:`T`.

Given values of :math:`T` and :math:`T_\mathrm{D}`, this module provides a function that solves 
the integral above numerically.

In addition, a dictionary of Debye temperatures at room temperature (`room_temperature_T_D`) 
from an online compilation 
:cite:`KnowledgeDoor2020` is provided.
The keys of the dictionary are the respective element symbols.
For example, for the Debye temperature of boron, one would type:

::

    from ries.resonance.debye_model import room_temperature_T_D
    
    print(room_temperature_T_D['B'])
"""

import numpy as np
from scipy.integrate import quad


def effective_temperature_debye_approximation(T, T_D):
    """Debye model for the effective temperature
    
    Given the thermodynamic temperature and the Debye temperature, this function calculates the 
    effective temperature of the material.
    The defining integral is solved numerically using `scipy.integrate.quad`.

    Parameters:

    - `T`, float or array_like, thermodynamic temperature in K.
    - `T_D`, float or array_like, Debye temperature in K.

    Returns:

    - float or array_like, effective temperature in K.
    
    Exceptions:

    - `ZeroDivisionError`, if a value of exactly 0 is entered for the thermodynamic temperature.
    """
    return (
        3.0
        * (T / T_D) ** 3
        * T
        * quad(lambda t: t ** 3 * (1.0 / (np.exp(t) - 1.0) + 0.5), 0.0, T_D / T)[0]
    )


# Source: https://www.knowledgedoor.com/2/elements_handbook/debye_temperature.html
room_temperature_T_D = {
    "Al": 390.0,
    "Am": None,
    "Sb": 200.0,
    "Ar": None,
    "As": 275.0,
    "Ba": 116.0,
    "Be": 1031.0,
    "Bi": 116.0,
    "B": 1362.0,
    "Cd": 221.0,
    "Ca": 230.0,
    "C": 1550.0,
    "Ce": 138.0,
    "Cs": 43.0,
    "Cr": 424.0,
    "Co": 386.0,
    "Cu": 310.0,
    "Cm": None,
    "Dy": 158.0,
    "Er": 163.0,
    "Eu": None,
    "Fr": None,
    "Gd": 155.0,
    "Ga": 240.0,
    "Ge": 403.0,
    "Au": 178.0,
    "Hf": 213.0,
    "Ho": 161.0,
    "H": None,
    "In": 129.0,
    "I": None,
    "Ir": 228.0,
    "Fe": 373.0,
    "Kr": None,
    "La": 135.0,
    "Pb": 87.0,
    "Li": 448.0,
    "Lu": 116.0,
    "Mg": 330.0,
    "Ma": 363.0,
    "Hg": 92.0,
    "Mo": 377.0,
    "Nd": 148.0,
    "Ne": None,
    "Np": 163.0,
    "Ni": 345.0,
    "Nb": 260.0,
    "Os": 400.0,
    "Pd": 275.0,
    "P": 576.0,
    "Pt": 225.0,
    "Pu": 176.0,
    "K": 100.0,
    "Pr": 138.0,
    "Pa": 262.0,
    "Ra": None,
    "Re": 275.0,
    "Rh": 350.0,
    "Rb": 59.0,
    "Ru": 415.0,
    "Sm": 184.0,
    "Sc": 476.0,
    "Se": None,
    "Si": 692.0,
    "Ag": 221.0,
    "Na": 155.0,
    "Sr": 148.0,
    "S": 527.0,
    "Ta": 225.0,
    "Tc": 422.0,
    "Te": None,
    "Tb": 158.0,
    "Tl": 96.0,
    "Th": 100.0,
    "Tm": 167.0,
    "Sn": 254.0,
    "Ti": 380.0,
    "W": 312.0,
    "U": 300.0,
    "Va": 390.0,
    "Xe": None,
    "Yb": None,
    "Y": 214.0,
    "Zn": 237.0,
    "Zr": 250.0,
}
