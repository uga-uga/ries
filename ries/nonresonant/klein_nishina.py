import numpy as np
from scipy.constants import physical_constants

from ries.nonresonant.nonresonant import Nonresonant

class KleinNishina(Nonresonant):
    """Cross section for Compton scattering of a photon off a free electron
    
This class implements the total, energy-differential, solid-angle differential, and 
scattering-angle differential cross section for Compton scattering, as given by the Klein-Nishina 
formula [1].

Note that this Klein-Nishina formula gives the cross section for the scattering off a single \
electron.
In low-energy nuclear physics, one is often interested in the cross section per atom to be able \
to compare this electromagnetic process to the cross sections for nuclear reactions.
To obtain the cross section for scattering off an atom, the cross section needs to be multiplied \
by the charge of the atom.

Attributes
----------
Z: int
    Proton number of the atom.

[1] O. Klein and Y. Nishina, 'Über die Streuung von Strahlung durch freie Elektronen nach der \
neuen relativistischen Quantendynamik von Dirac', Z. Phys. 52, 8530868 (1929) https://doi.org/10.1007/BF01366453
    """
    def __init__(self, Z=1):
        """Initialization
        
Parameters
----------
Z: int
    Proton number of the atom. The default value is 1, which corresponds to scattering off a \
single electron.
        """
        self.Z = Z

    def __call__(self, energy):
        """Total cross section

This method has been implemented for convenience and basically calls KleinNishina.cs_total(). \
An additional parameter, the proton number :math:`Z` has been added. \
This is because a user may be interested in the cross section 'per atom' instead of \
the one 'per electron'. \
In this case, the cross section needs to be multiplied by :math:`Z`.

Parameters
----------
e0: array_like or scalar
    Initial energy of the photon in MeV.
Z: int
    Proton number of the scattering atom (default: 1, i.e. scattering on a single electron).
    
Returns
-------
sigma: ndarray or scalar
    Total cross section in fm**2.
        """
        return self.cs_total(energy)

    def compton_edge(self, e0):
        """Compton edge for a given initial photon energy
        
        The Compton edge is the lower limit for the energy :math:`E` of the photon after the \
    scattering process, which follows from the conservation of energy and linear momentum. \
    It corresponds to scattering by 180°, which means maximum momentum transfer to the electron. \
    The expression for the Compton edge is given by:

        .. math:: \frac{E_0}{1+2 \frac{E_0}{m c^2}}
        
        Here, :math:`m` denotes the mass of the charged particle.
        
        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
        
        Returns
        -------
        ec: array_like or scalar
            Energy of the Compton edge in MeV.
        """
        
        return e0/(1+2*e0/physical_constants['electron mass energy equivalent in MeV'][0])
        
    def theta(self, e0, e):
        """Scattering angle for a given energy transfer
        
        This function relates the scattering angle :math:`\theta` of the photon to its \
    initial energy :math:`E_0` and the scattered energy :math:`E`. \
    For the Compton effect, this relation is unique, and it is given by:

        ..math:: \\theta = \\arccos \\left[ 1 - \\left( \\frac{E_0}{E} - 1 \\right) \\frac{m c^2}{E_0}  \\right]
        
        Here, :math:`m` denotes the mass of the charged particle.
        
        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
        e: array_like or scalar
            Energy of the photon after the scattering process in MeV.

        Returns
        -------
        theta: ndarray or scalar
            Polar scattering angle of the photon in radians.
        """
        return (np.arccos(
            1-(e0/e-1)*physical_constants['electron mass energy equivalent in MeV'][0]/e0))

    def e_over_e0(self, e0, theta):
        """Relative energy change of the photon for a given scattering angle
        
        This function relates the scattering angle :math:`\theta` to the ratio of \
    the inital (:math:`E_0`) and the final (:math:`E`) energy of the photon. \
    For the Compton effect, this relation is unique, and it is given by:

        ..math:: \\frac{E}{E_0} = \\frac{1}{1 + \\frac{E_0}{m c^2} \\left[ 1 - \\cos \\left( \\theta \\right) \\right]}
        
        Here, :math:`m` denotes the mass of the charged particle.

        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
        theta: array_like or scalar
            Polar scattering angle of the photon in radians.
            
        Returns
        -------
        e/e0: ndarray or scalar
            Ratio of the energy of the scattered photon and its initial energy.
        """
        return (1./(1.+e0/physical_constants['electron mass energy equivalent in MeV'][0]
                    *(1.-np.cos(theta))))

    def cs_diff(self, e0, theta):
        """Differential cross section w.r.t. the solid-angle
        
        This function implements the so-called Klein-Nishina differential cross section for the \
    scattering of a photon with an initial energy :math:`E_0` off a point-like, charged particle. \
    The scattered photon has an energy :math:`E`, and is scattered at a polar angle :math:`'\theta`. \
    Note that the cross section does not depend on the azimuthal angle :math:`\\varphi`. \
    Explicitly, the differential cross section is given by:
        
        .. math:: \\frac{\\mathrm{d} \\sigma}{\\mathrm{d}\\Omega} = \\frac{1}{2}\\frac{\\alpha^2 \\hbar^2}{m^2 c^2} x^2 \\left[ x + \\frac{1}{x} - \\sin \\left( \\theta^2 \\right) \\right]
        
        Here, :math:`\\sigma` denotes the cross section, :math:`\\Omega` denotes the \
        solid-angle element, :math:`\\alpha` is the fine-structure constant, :math:`m` is the 
        mass of the charged particle, and :math:`x` is an abbreviation for the ratio :math:`E/E_0`.
        
        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
        theta: array_like or scalar
            Polar scattering angle of the photon in radians.
            
        Returns
        -------
        d(sigma)/d(Omega) : ndarray or scalar
            Solid-angle differential cross section for the given initial energy and \
    scattering angle in fm**2.
        """
        relative_energy_change = self.e_over_e0(e0, theta)
        return (self.Z*np.pi*physical_constants['fine-structure constant'][0]**2
                *physical_constants['Planck constant over 2 pi times c in MeV fm'][0]**2
                /(physical_constants['electron mass energy equivalent in MeV'][0]**2)
                *relative_energy_change*relative_energy_change*(
                    relative_energy_change + 1./relative_energy_change - np.sin(theta)**2
                )
            )

    def cs_diff_de(self, e0, e):
        """Differential cross section w.r.t. the energy of the scattered photon
        
        This function implements the energy-differential cross section, which can be obtained \
    from the solid-angle differential cross section as follows: \
    First, the dependence on :math:`\\varphi` can be integrated out:

        .. math: \\int_{0}^{2 \\pi} \\frac{\\mathrm{d} \\sigma}{\\mathrm{d} \\Omega} \\mathrm{d} \\varphi = \\frac{\\mathrm{d} \\sigma}{\\mathrm{d} \\cos \\left( \\theta \\right)}

    The differential :math:`\\mathrm{d} \\cos \\left( \\theta \\right)` can be uniquely related \
    to the energy :math:`E` of the scattered photon [see theta() function], which yields:

        .. math: \\mathrm{d} \\cos \\left( \\theta \\right) = \\frac{m c^2}{E^2} \\mathrm{d} E
        
        See also the function cs_diff_dtheta.
        
        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
        e: array_like or scalar
            Energy of the photon after the scattering process in MeV.

        Returns
        -------
        d(sigma)/d(E) d(E): ndarray or scalar
            Energy-differential cross section for the given initial and final energies \
    in fm**2 MeV**-1.
        
        """
        return (self.cs_diff(e0, self.theta(e0, e))
                *physical_constants['electron mass energy equivalent in MeV'][0]/(e*e))

    def cs_diff_dtheta(self, e0, theta):
        """Differential cross section w.r.t. the polar scattering angle
        
        This function implements the  polar scattering-angle differential cross section, \
    which can be obtained from the solid-angle differential cross section by expanding the \
    differential :math:`\\mathrm{d} \\Omega` in the two angular variables of the sperical \
    coordinate system, the polar angle :math:`\\theta` and the azimuthal angle \
    :math:`\\varphi`:

        ..math:: \\mathrm{d} \\Omega = \\sin \\left( \\theta \\right) \\mathrm{d} \\theta \\mathrm{d} \\varphi.

    Since the cross section is independent of the azimuthal angle

        ..math:: \\frac{\\mathrm{d} \\sigma}{\\mathrm{d} \\varphi} = 0,
        
    the polar-angle differential cross section is given by:

        ..math:: \\int_0^\\pi \\mathrm{d} \\theta \\int_0^{2\\pi} \\mathrm{d} \\varphi \\frac{\\mathrm{d}\\sigma}{\\mathrm{d} \\cos \\left( \\theta \\right) \\mathrm{d} \\varphi} \\sin \\left( \\theta \\right) \\
        = \\int_0^\\pi \\mathrm{d} \\theta \\frac{\\mathrm{d}\\sigma}{\\mathrm{d} \\cos \\left( \\theta \\right)} \\sin \\left( \\theta \\right),

        i.e. the integration over :math:`\\varphi` can simply be dropped, because it results in a :math:`\\theta`-dependent integration constant, which is identified as \
    :math:`\\mathrm{d} \\sigma / \\mathrm{d} \\cos \\left( \\theta \\right)`.

        See also the function cs_diff_de.

        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
        theta: array_like or scalar
            Polar scattering angle of the photon in radians.
            
        Returns
        -------
        d(sigma)/d(theta): ndarray or scalar
            Scattering-angle differential cross section for the given initial energy and scattering angle \
    in fm**2.
        """

        return (self.cs_diff(e0, theta)
                *np.sin(theta) # d(cos(theta)) = sin(theta) dtheta
            )

    def cs_total(self, e0):
        """Total cross section
        
        This function implements the total cross section for the Compton scattering of a photon with a given initial energy :math:`E_0`.
        The total cross section is obtained by integrating the energy- (scattering-angle) differential cross section over the \
    energy (scattering angle) from 0 to infinity (0 to pi). An analytical form for this integral exists, which can be used for \
    cross checks.

        Parameters
        ----------
        e0: array_like or scalar
            Initial energy of the photon in MeV.
            
        Returns
        -------
        sigma: ndarray or scalar
            Total cross section in fm**2.

        """
        x = e0/physical_constants['electron mass energy equivalent in MeV'][0]
        return (self.Z*np.pi*physical_constants['fine-structure constant'][0]**2
                *physical_constants['Planck constant over 2 pi times c in MeV fm'][0]**2
                /(physical_constants['electron mass energy equivalent in MeV'][0]**2*x**3)
                *((2*x*(2+x*(1+x)*(8+x)))/((1+2*x)**2)+((x-2)*x-2)*np.log(1+2*x))
            )