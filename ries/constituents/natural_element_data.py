from ries.constituents.isotope import Isotope

class NISTElementDataReader:
    def __init__(self, element_data_file_name):
        self.element_data_file_name = element_data_file_name

        self.Z_prefix = 'Atomic Number = '
        self.A_prefix = 'Mass Number = '
        self.X_prefix = 'Atomic Symbol = '
        self.amu_prefix = 'Relative Atomic Mass = '
        self.abundance_prefix = 'Isotopic Composition = '

    def read_nist_element_data(self, Z, X):

        abundances = {}
        isotopes = {}
        with open(self.element_data_file_name, 'r') as file:
            for line in file:
                if self.Z_prefix in line:
                    current_Z = self.read_nist_element_property(line, self.Z_prefix, int)
                    if current_Z < Z:
                        continue
                    elif current_Z > Z:
                        break
                elif self.A_prefix in line and current_Z == Z:
                    A = self.read_nist_element_property(line, self.A_prefix, int)
                elif self.amu_prefix in line and current_Z == Z:
                    amu = self.read_nist_element_property_with_uncertainty(line, self.amu_prefix)
                elif self.abundance_prefix in line and current_Z == Z:
                    abundance = self.read_nist_element_property_with_uncertainty(line, self.abundance_prefix, default=1.)
                    AX = '{:d}{}'.format(A, X)
                    abundances[AX] = abundance
                    isotopes[AX] = Isotope(AX, amu)
        
        return (abundances, isotopes)

    def read_nist_element_property(self, line, prefix, property_type=str, default=None):
        prop = line[len(prefix):-1] # Take all characters except the last two, which are the newline
        # escape sequence '\n'.
        if prop != '':
            return property_type(prop)
        return default

    def read_nist_element_property_with_uncertainty(self, line, prefix, default=None):
        return self.read_nist_element_property(
            line[0:line.find('(')+1], prefix, float, default
        )

    def read_nist_element_symbols(self):
        X = {}
        with open(self.element_data_file_name, 'r') as file:
            for line in file:
                if self.Z_prefix in line:
                    Z = self.read_nist_element_property(line, self.Z_prefix, int)
                if self.X_prefix in line and not Z in X:
                    X[Z] = self.read_nist_element_property(line, self.X_prefix)
        return X