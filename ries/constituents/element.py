from pathlib import Path

import numpy as np
from scipy.constants import physical_constants

from ries.constituents.isotope import Isotope
from ries.nonresonant.xrmac import XRMAC 

class Element:
    def __init__(self, Z, X=None, isotopes=None, abundances=None, 
    xrmac=None):
        self.Z = Z
        self.X = X
        if isinstance(isotopes, str):
            abundances, isotopes = self.read_nist_element_data(isotopes)
        self.isotopes = isotopes
        self.abundances = abundances
        self.amu = self.amu_from_isotopic_composition(self.abundances, self.isotopes)
        self.xrmac = xrmac

    def amu_from_isotopic_composition(self, abundances, isotopes):
        return np.sum([isotopes[iso].amu*abundances[iso] for iso in isotopes])

    def read_nist_element_data(self, element_data_file_name):

        abundances = {}
        isotopes = {}
        with open(element_data_file_name, 'r') as file:
            for line in file:
                if Z_prefix in line:
                    Z = read_nist_element_property(line, Z_prefix, int)
                    if Z < self.Z:
                        continue
                    elif Z > self.Z:
                        break
                elif A_prefix in line and Z == self.Z:
                    A = read_nist_element_property(line, A_prefix, int)
                elif amu_prefix in line and Z == self.Z:
                    amu = read_nist_element_property_with_uncertainty(line, amu_prefix)
                elif abundance_prefix in line and Z == self.Z:
                    abundance = read_nist_element_property_with_uncertainty(line, abundance_prefix, default=1.)
                    AX = '{:d}{}'.format(A, self.X)
                    abundances[AX] = abundance
                    isotopes[AX] = Isotope(AX, amu)
        
        return (abundances, isotopes)

Z_prefix = 'Atomic Number = '
A_prefix = 'Mass Number = '
X_prefix = 'Atomic Symbol = '
amu_prefix = 'Relative Atomic Mass = '
abundance_prefix = 'Isotopic Composition = '

def read_nist_element_property(line, prefix, property_type=str, default=None):
    prop = line[len(prefix):-1] # Take all characters except the last two, which are the newline
    # escape sequence '\n'.
    if prop != '':
        return property_type(prop)
    return default

def read_nist_element_property_with_uncertainty(line, prefix, default=None):
    return read_nist_element_property(
        line[0:line.find('(')+1], prefix, float, default
    )

def read_nist_element_symbols(element_data_file_name):
    X = {}
    with open(element_data_file_name, 'r') as file:
        for line in file:
            if Z_prefix in line:
                Z = read_nist_element_property(line, Z_prefix, int)
            if X_prefix in line and not Z in X:
                X[Z] = read_nist_element_property(line, X_prefix)
    return X

element_data_file_name = Path(__file__).parent.absolute() / 'nist_elements/elements.txt'
xrmac_data_dir = Path(__file__).parent.absolute() / '../nonresonant/nist_xrmac/'

cm_to_fm = 1e13
kg_to_g = 1e3
X = read_nist_element_symbols(element_data_file_name)
natural_elements = {}
for Z in range(1, 119):
    natural_elements[X[Z]] = Element(
        Z, X[Z],
        isotopes=str(element_data_file_name),
    )

    natural_elements[X[Z]].xrmac = XRMAC(
        str(xrmac_data_dir / '{:02d}.txt'.format(Z)),
        xrmac_conversion=lambda xrmac: xrmac*cm_to_fm**2*natural_elements[X[Z]].amu*physical_constants['atomic mass constant'][0]*kg_to_g
    ) if Z < 93 else None # No mass attenuation data for elements with Z > 92 are available from the NIST.