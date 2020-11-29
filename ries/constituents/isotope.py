"""
Module for storing isotope data.

An isotope :math:`^A\\mathrm{X}` with a mass number :math:`\\mathrm{A}` and an element symbol \
:math:`\\mathrm{X}` has the following properties:

* an isotopic mass :math:`m \\left( ^A\\mathrm{X} \\right)`
* a ground state :math:`J_0^{\\pi_0}`
* a set of excited states :math:`\\left\{ J_i^{\\pi_i}\\right\}` (:math:`i > 0`)

For a detailed example of how to create a user-defined element, see `tests/boron.py`.
"""

from scipy.constants import physical_constants

class Isotope:
    """Class representing an isotope
    
Attributes:

- `AX`, str, mass number and element symbol. \
For example, the isotope of lead with a neutron number of 126 would be '208Pb'.
- `amu`, float, isotopic mass in atomic mass units (AMU).
- `ground_state`, `State` object, ground state of the isotope.
- `excited_states`, array of `State` objects, list of excited states of the isotope.
    """
    def __init__(self, AX, amu, 
        ground_state=None, excited_states=[]):
        """Initialization

Parameters:

- `AX`, str, mass number and element symbol. \
For example, the isotope of lead with a neutron number of 126 would be '208Pb'.
- `amu`, float, isotopic mass in atomic mass units (AMU).
- `ground_state`, `State` object, ground state of the isotope (default: None).
- `excited_states`, array of `State` objects, list of excited states of the isotope \
(default: empty list).
        """
        self.AX = AX
        self.amu = amu
        self.ground_state = ground_state
        self.excited_states = excited_states