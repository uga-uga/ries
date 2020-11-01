from pathlib import Path

import numpy as np
from scipy.constants import physical_constants

from ries.constituents.isotope import Isotope
from ries.constituents.natural_element_data import NISTElementDataReader

class Element:
    def __init__(self, Z, X, isotopes, abundances, 
    xrmac=None):
        self.Z = Z
        self.X = X
        self.isotopes = isotopes
        self.abundances = abundances
        self.amu = self.amu_from_isotopic_composition(self.abundances, self.isotopes)

    def amu_from_isotopic_composition(self, abundances, isotopes):
        return np.sum([isotopes[iso].amu*abundances[iso] for iso in isotopes])

nist_element_data_reader = NISTElementDataReader(Path(__file__).parent.absolute() / 'nist_elements/elements.txt')
X = nist_element_data_reader.read_nist_element_symbols()

natural_elements = {}
for Z in range(1, 119):
    abundances, isotopes = nist_element_data_reader.read_nist_element_data(Z, X[Z])
    natural_elements[X[Z]] = Element(
        Z, X[Z],
        isotopes, abundances,
    )