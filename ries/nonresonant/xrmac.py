from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

class XRMAC:
    def __init__(self, data):
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
                    float(line[0]), float(line[1]), float(line[2])])
        return np.array(data)