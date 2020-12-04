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

"""
Module for storing properties of nuclear states.

Any state with an index :math:`i` \
(:math:`i = 0` denotes the ground state, :math:`i > 0` denotes an excited state) \
of a nucleus is characterized by its:

* total angular momentum quantum number ('spin') :math:`J_i`
* parity quantum number :math:`\\pi_i`
* excitation energy :math:`E_i`

Here, it is assumed that if :math:`i > j`, then :math:`E_i > E_j`.
In addition, if the state is not the ground state [c]_, it has the properties:

* total width :math:`\\Gamma_i`
* partial widths for the decay to lower-lying states :math:`\\Gamma_{i \\to j}` (:math:`j < i`)

All state objects (`GroundState`, `State`) have an attribute `J_pi` which is supposed to be their
unique identifier.
For example, the states may simply be labeled by their energies, i.e. the ground state is `0`, \
the first excited state is `1`, and so on.
In the examples of `ries`, the states are labeled by their spin, parity, and the index of that \
combination, i.e. the third state with a spin 2 and a positive parity would be `J_pi = '2^+_3'`.
With these unique identifiers, the partial widths can be implemented as a dictionary with the key \
`J_pi` (see `tests/boron.py` for a detailed example).

By itself, the `State` class makes no assumptions about the units of its properties, but it is \
recommended to give energies and widths in :math:`\\mathrm{MeV}`.

.. [c] For the sake of simplicity, it is assumed that the ground state has no decays, which is \
of course not true in reality. \
However, since photonuclear reactions with real photons are often only feasible with (quasi-) \
stable nuclei, this is a valid approximation.
"""

from scipy.constants import physical_constants

class GroundState:
    """Class representing a ground state

Attributes:

- `J_pi`, str, unique identifier.
- `two_J`, int, two times the total angular momentum quantum number.
- `parity`, int, parity quantum number (`1` means 'positive' and `-1` means 'negative').
- `excitation_energy`, float, excitation energy with respect to the ground state \
(since this class represents the ground state, the excitation energy is zero).
    """
    def __init__(self, J_pi, two_J, parity):
        """Initialization
        
Parameters:

- `J_pi`, str, unique identifier.
- `two_J`, int, two times the total angular momentum quantum number.
- `parity`, int, parity quantum number (`1` means 'positive' and `-1` means 'negative').
        """
        self.J_pi = J_pi
        self.two_J = two_J
        self.parity = parity
        self.excitation_energy = 0.

class State(GroundState):
    """Class representing an excited state

Attributes:

- `J_pi`, str, unique identifier.
- `two_J`, int, two times the total angular momentum quantum number.
- `parity`, int, parity quantum number (`1` means 'positive' and `-1` means 'negative').
- `excitation_energy`, float, excitation energy with respect to the ground state.
- `partial_widths`, dictionary, list of partial widths to lower-lying states. \
The list should be given as a dictionary with `J_pi` of the final states of the decay as keys, \
and the partial widths as values.
- `width`, float, total width. \
This property is inferred from the given partial widths.
    """
    def __init__(self, J_pi, two_J, parity, excitation_energy, partial_widths):
        """Initialization

The total width is inferred from the given partial widths.

Parameters:

- `J_pi`, str, unique identifier.
- `two_J`, int, two times the total angular momentum quantum number.
- `parity`, int, parity quantum number (`1` means 'positive' and `-1` means 'negative').
- `excitation_energy`, float, excitation energy with respect to the ground state.
- `partial_widths`, dictionary, list of partial widths to lower-lying states. \
The list should be given as a dictionary with `J_pi` of the final states of the decay as keys, \
and the partial widths as values.
        """
        GroundState.__init__(self, J_pi, two_J, parity)
        self.excitation_energy = excitation_energy
        self.partial_widths = partial_widths
        self.width = sum([self.partial_widths[state] for state in self.partial_widths])