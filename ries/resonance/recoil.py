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
Recoil corrections to the resonance energy of a photoabsorption process.

At zero-th order, the resonance energy :math:`E_r` for the absorption of a photon is given by the \
energy difference between the initial ('0') and the excited state ('2'):

.. math:: E_r = E_2 - E_0

Due to the conservation of momentum in the absorption process, some of the energy of the photon \
must be transferred to the nucleus or the atomic lattice as kinetic (translational, rotational, \
vibrational, ...) energy.
The recoil correction :math:`R`, introduced as

.. math:: E_r = E_2 - E_0 + R

is usually more than 5 orders of magnitude smaller than the \
resonance energy, but important in many types of experiments.
For example, nuclear resonance fluorescence as an analog of atomic resonance fluorescence could \
not be observed in the laboratory for decades due to the recoil problem :cite:`Metzger1959` :cite:`KneisslPitzZilges1996`.

At the moment, the only nonvanishing recoil correction that this module implements is the one for \
a free nucleus.
However, be aware that the situation for an atom in a condensed-matter state can be drastically \
different :cite:`Lamb1939` :cite:`Moessbauer1958`.
"""
from scipy.constants import physical_constants

class Recoil:
    """Abstract class for a recoil correction"""
    def __init__(self):
        """Initialization"""
        pass

    def __call__(self, energy_difference):
        """Recoil-corrected resonance energy
        
Paramters:

- `energy_difference`, float, energy difference between two nuclear states, which would be the \
resonance energy for an infinitely heavy nucleus.

Returns:

- float, recoil-corrected resonance energy.
        """
        return energy_difference

class NoRecoil(Recoil):
    """Dummy class for negligible recoil correction"""
    def __init__(self):
        """Initialization"""
        Recoil.__init__(self)

class FreeNucleusRecoil(Recoil):
    """Recoil correction for a free nucleus at rest

For a free nucleus at rest, the recoil correction can be derived from 4-momentum conservation:

.. math:: p_\\mu = p_\\mu^\\prime,

where :math:`p_\\mu` and :math:`p_\\mu^\\prime` are the 4 vectors of the system before and after \
the reaction, respectively.
Assuming that the mass of the nucleus is :math:`m \\left( ^A\\mathrm{X}\\right)` and the incident \
photon has the correct energy :math:`E_r` to excite the resonance, the 4-momentum before the \
reaction is (assuming a propagation along the :math:`z` axis w.l.o.g.):

.. math:: p_\\mu = \\left( \\frac{E_r + m \\left( ^A\\mathrm{X}\\right) c^2}{c}, 0, 0, \\frac{E_r}{c} \\right).

After the absorption, the nucleus will be in the excited state and propagating along the \
:math:`z` axis:

.. math:: p_\\mu^\\prime = \\left( \\frac{\\sqrt{ \\left[ m \\left( ^A\\mathrm{X}\\right) c^2 + \\left( E_2 - E_0 \\right) \\right]^2 + E_r^2}}{c}, 0, 0, \\frac{E_r}{c} \\right).

Finding the :math:`z` component of the momentum after the absorption is trivial.
From the conservation of energy, one obtains:

.. math:: E_r = E_2 - E_0 + \\underbrace{\\frac{\\left( E_2 - E_0 \\right)^2}{2 m \\left( ^A\\mathrm{X}\\right) c^2}}_{R}.

Attributes:

- `amu`, float, mass of the nucleus in atomic mass units.
    """
    def __init__(self, amu):
        """Initialization
        
Parameters:

- `amu`, float, mass of the nucleus in atomic mass units.
        """
        self.amu = amu

    def __call__(self, energy_difference):
        """Recoil-corrected resonance energy
        
See `Recoil.__call__()`.
        """
        return (
            energy_difference*(
                1.+energy_difference
                /(2.*self.amu*physical_constants['atomic mass constant energy equivalent in MeV'][0])
            )
        )