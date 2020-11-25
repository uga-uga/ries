"""
At first order in quantum electrodynamics, the cross section \
:math:`\\sigma_\\mathrm{Compton} \\left( E \\right)` for Compton scattering of a photon with an \
initial energy :math:`E` off a free electron at rest is given by the Klein-Nishina formula \
:cite:`KleinNishina1929`:

.. math:: \\sigma_\\mathrm{Compton} \\left( E \\right) = \\frac{\\pi \\alpha^2 \\hbar^2}{m_e^2 c^2} \\frac{1}{x^3} \\left\\{ \\frac{2x \\left[ 2 + x \\left( 1 + x \\right) \\left( 8 + x \\right)\\right]}{\\left( 1 + 2 x \\right)^2} + \\left[ \\left( x - 2 \\right) x - 2 \\right] \\mathrm{log} \\left( 1 + 2x \\right) \\right\\}.

In the equation above, the abbreviation

.. math:: x = \\frac{E}{m_e c^2}

has been used.
The symbols :math:`\\alpha`, :math:`\\hbar`, :math:`m_e`, and :math:`c` denote the fine-structure \
constant, the reduced Planck constant, the electron rest mass, and the speed of light, \
respectively.

The corresponding solid-angle differential cross section is given by :cite:`Depaola2003`:

.. math:: \\frac{\\mathrm{d} \\sigma_\\mathrm{Compton}}{\\mathrm{d}\\Omega} (E, \\Omega) = \\frac{1}{2}\\frac{\\alpha^2 \\hbar^2}{m_e^2 c^2} \\left( \\frac{E^\\prime}{E} \\right)^2 \\left[ \\frac{E^\\prime}{E} + \\frac{E}{E^\\prime} - 2 \\sin \\left( \\theta \\right)^2 \\cos \\left( \\varphi \\right)^2 \\right].

Here, the angle :math:`\\theta` denotes the polar angle of the scattered photon with respect to \
its initial direction of propagation.
The angle :math:`\\varphi` is the corresponding azimuthal angle, which is  defined with respect \
to the direction of polarization of the initial photon.
The equation above also contains the energy of the scattered photon, :math:`E^\\prime`, which is related to :math:`E` via:

.. math:: E^\\prime \\left( E, \\theta \\right) = \\frac{E}{1 + \\frac{E}{m_e c^2} \\left[ 1 - \\cos \\left(\\theta \\right) \\right]}.

If no information on the polarization of the incident photon is available, the differential \
cross section can be averaged over the azimuthal angle:

.. math:: \\frac{\\mathrm{d} \\sigma_\\mathrm{Compton, unpolarized}}{\\mathrm{d}\\Omega} (E, \\Omega) = \\frac{\\int_0^{2\\pi} \\frac{\\mathrm{d} \\sigma_\\mathrm{Compton}}{\\mathrm{d}\\Omega} (E, \\Omega) \\mathrm{d} \\varphi}{\\int_0^{2\\pi} \\mathrm{d} \\varphi} = \\frac{1}{2}\\frac{\\alpha^2 \\hbar^2}{m_e^2 c^2} \\left( \\frac{E^\\prime}{E} \\right)^2 \\left[ \\frac{E^\\prime}{E} + \\frac{E}{E^\\prime} - \\sin \\left( \\theta \\right)^2 \\right].

An integration of either :math:`\\mathrm{d} \\sigma_\\mathrm{Compton} / \\mathrm{d} \\Omega` or \
:math:`\\mathrm{d} \\sigma_\\mathrm{Compton, unpolarized} / \\mathrm{d} \\Omega` over the azimuthal \
angle gives the scattering-angle differential cross section:

.. math:: \\frac{\\mathrm{d} \\sigma_\\mathrm{Compton}}{\\mathrm{d} \\theta} (E, \\theta) = \\int_0^{2\\pi} \\frac{\\mathrm{d} \\sigma_\\mathrm{Compton}}{\\mathrm{d}\\Omega} (E, \\Omega) \\sin \\left( \\theta \\right) \\mathrm{d} \\varphi = \\frac{\\pi \\alpha^2 \\hbar^2}{m_e^2 c^2} \\left( \\frac{E^\\prime}{E} \\right)^2 \\left[ \\frac{E^\\prime}{E} + \\frac{E}{E^\\prime} - \\sin \\left( \\theta \\right)^2 \\right] \\sin \\left( \\theta \\right).

Note the additional factor :math:`\\sin \\left( \\theta \\right)`, which comes from the fact that:

.. math:: \\mathrm{d} \\Omega = \\mathrm{d} \\cos \\left( \\theta \\right) \\mathrm{d} \\varphi =  \\sin \\left( \\theta \\right) \\mathrm{d} \\theta \\mathrm{d} \\varphi.

The differential cross section with respect to :math:`\\cos \\left( \\theta \\right)` can be \
converted to the energy-differential cross section using the relation to :math:`E^\\prime` above:

.. math:: \\mathrm{d} \\cos \\left( \\theta \\right) = \\frac{m c^2}{E^{\\prime 2}} \\mathrm{d} E^\\prime.

Inserting the relation above in the expression for the scattering-angle differential equation \
and abbreviating the inverse relation between :math:`\\theta` and :math:`E^\\prime` as \
:math:`\\theta \\left( E, E^\\prime \\right)`, one obtains:

.. math:: \\frac{\\mathrm{d} \\sigma_\\mathrm{Compton}}{\\mathrm{d} E^\\prime} (E, E^\\prime) = \\frac{\\pi \\alpha^2 \\hbar^2}{E^2} \\left\\{ \\frac{E^\\prime}{E} + \\frac{E}{E^\\prime} - \\sin \\left[ \\theta \\left( E, E^\\prime \\right) \\right]^2 \\right\\}.
"""

import numpy as np
from scipy.constants import physical_constants

from ries.nonresonant.nonresonant import Nonresonant

class KleinNishina(Nonresonant):
    """(Differential) Cross section for Compton scattering of a photon off a free electron

Note that the Klein-Nishina formula gives the cross section for the scattering off a single \
electron.
In low-energy nuclear physics, one is often interested in the cross section per atom, to be able \
to compare this electromagnetic process to the cross sections for nuclear reactions.
To obtain the cross section for scattering off an atom, the cross section needs to be multiplied \
by the number of electrons per atom.

Attributes:

- `Z`: int, number of electrons per atom, which is the proton number for a neutral atom.
    """
    def __init__(self, Z=1):
        """Initialization
        
Parameters:

- `Z`: int, number of electrons per atom, which is the proton number for a neutral atom.
        """
        self.Z = Z

    def __call__(self, E):
        """Total cross section

This method has been implemented for convenience and calls `KleinNishina.cs_total()`.

Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.

Returns:

- array_like or scalar, :math:`\\sigma_\\mathrm{Compton}` in :math:`\\mathrm{fm}^2`.
        """
        return self.cs_total(E)

    def compton_edge(self, E):
        """Compton edge for a given initial photon energy

The Compton edge is the lower limit for the energy :math:`E^\\prime` of the photon after the \
scattering process.
It follows from the conservation of energy and linear momentum. \
The limit corresponds to scattering by :math:`180^\\circ`, which means maximum momentum transfer to the electron. \
The expression for the Compton edge is given by:

.. math:: E^\\prime \\left( E, 180^\\circ \\right) = \\frac{E}{1+2 \\frac{E}{m_e c^2}}
        
Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.

Returns:

- array_like or scalar, energy of the Compton edge in MeV.
        """
        
        return E/(1+2*E/physical_constants['electron mass energy equivalent in MeV'][0])
        
    def theta(self, E, Ep):
        """Scattering angle for a given scattered energy
        
Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.
- Ep: array_like or scalar, energy of the photon after the scattering process in MeV.

Returns

- array_like or scalar, scattering angle of the photon in radians.
        """
        return (np.arccos(
            1-(E/Ep-1)*physical_constants['electron mass energy equivalent in MeV'][0]/E))

    def Ep_over_E(self, E, theta):
        """Ratio of final and initial photon energy for a given scattering angle

Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.
- theta: array_like or scalar, scattering angle of the photon in radians.
            
Returns:

- array_like or scalar, ratio of the energy of the scattered photon and its initial energy.
        """
        return (1./(1.+E/physical_constants['electron mass energy equivalent in MeV'][0]
                    *(1.-np.cos(theta))))

    def cs_diff(self, e0, theta):
        """Differential cross section w.r.t. the solid-angle
        
Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.
- theta: array_like or scalar, scattering angle of the photon in radians.

Returns:

- array_like or scalar, solid-angle differential cross section in fm**2.
        """
        relative_energy_change = self.Ep_over_E(e0, theta)
        return (self.Z*np.pi*physical_constants['fine-structure constant'][0]**2
                *physical_constants['Planck constant over 2 pi times c in MeV fm'][0]**2
                /(physical_constants['electron mass energy equivalent in MeV'][0]**2)
                *relative_energy_change*relative_energy_change*(
                    relative_energy_change + 1./relative_energy_change - np.sin(theta)**2
                )
            )

    def cs_diff_de(self, E, Ep):
        """Differential cross section w.r.t. the energy of the scattered photon
        
Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.
- Ep: array_like or scalar, energy of the photon after the scattering process in MeV.


Returns:

- array_like or scalar, energy-differential cross section in :math:`\\mathrm{fm}^2 \\mathrm{MeV}^{-1}`.
        
        """
        return (self.cs_diff(E, self.theta(E, Ep))
                *physical_constants['electron mass energy equivalent in MeV'][0]/(Ep*Ep))

    def cs_diff_dtheta(self, E, theta):
        """Differential cross section w.r.t. the polar scattering angle

Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.
- theta: array_like or scalar, scattering angle of the photon in radians.
            
Returns:

- array_like or scalar, scattering-angle differential cross section in :math:`\\mathrm{fm}^2`.
        """

        return (self.cs_diff(E, theta)
                *np.sin(theta)
            )

    def cs_total(self, e0):
        """Total cross section

Parameters:

- E: array_like or scalar, initial energy of the photon in MeV.
            
Returns:

array_like or scalar, total cross section in :math:`\\mathrm{fm}^2`.

        """
        x = e0/physical_constants['electron mass energy equivalent in MeV'][0]
        return (self.Z*np.pi*physical_constants['fine-structure constant'][0]**2
                *physical_constants['Planck constant over 2 pi times c in MeV fm'][0]**2
                /(physical_constants['electron mass energy equivalent in MeV'][0]**2*x**3)
                *((2*x*(2+x*(1+x)*(8+x)))/((1+2*x)**2)+((x-2)*x-2)*np.log(1+2*x))
            )