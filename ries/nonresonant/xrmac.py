from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.constants import physical_constants

from ries.constituents.element import natural_elements, X
from ries.cross_section import CrossSection
from ries.nonresonant.nonresonant import Nonresonant

class XRMAC(Nonresonant):
    def __init__(self, data,
        energy_conversion=lambda energy: energy,
        xrmac_conversion=lambda xrmac: xrmac
        ):
        self.energy_conversion = energy_conversion
        self.xrmac_conversion = xrmac_conversion
        if isinstance(data, str):
            data = self.read_nist_xrmac(data)
        self.data = data
        self.interpolation_log_log = self.interpolate_log_log(
            self.data)

    def __call__(self, energy):
        return 10**(self.interpolation_log_log(np.log10(energy)))

    def interpolate_log_log(self, data):
        return interp1d(
            np.log10(data[:,0]), np.log10(data[:,1]),
            bounds_error=False,
            fill_value=(np.log10(data[0][1]), np.log10(data[-1][1]))
        )

    def read_nist_xrmac(self, xrmac_file_name):
        data = []
        skip = 0
        # Some data files have a three-character column that indicates the label of the atomic 
        # resonance.
        # Open the file a first time to find out whether this is the case.
        # If yes, skip this column when reading.
        with open(xrmac_file_name, 'r') as file:
            if not file.readline()[0].isdigit():
                skip = 3
        with open(xrmac_file_name, 'r') as file:
            for line in file:
                line = line[skip:-2].split(sep='  ')
                data.append([
                        self.energy_conversion(float(line[0])),
                        self.xrmac_conversion(float(line[1])),
                        self.xrmac_conversion(float(line[2]))
                    ]
                )
        return np.array(data)

xrmac_data_dir = Path(__file__).parent.absolute() / '../nonresonant/nist_xrmac/'

xrmac = {}
cm_to_fm = 1e13
kg_to_g = 1e3

for Z in range(1, 93):
    xrmac[X[Z]] = XRMAC(
        str(xrmac_data_dir / '{:02d}.txt'.format(Z)),
        xrmac_conversion=lambda xrmac: xrmac*cm_to_fm**2*natural_elements[X[Z]].amu*physical_constants['atomic mass constant'][0]*kg_to_g
    )
