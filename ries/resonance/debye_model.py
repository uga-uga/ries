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

import numpy as np
from scipy.integrate import quad


def effective_temperature_debye_approximation(T, T_D):
    return (
        3.0
        * (T / T_D) ** 3
        * T
        * quad(lambda t: t ** 3 * (1.0 / (np.exp(t) - 1.0) + 0.5), 0.0, T_D / T)[0]
    )


room_temperature_T_D = {
    "Al": 390.0,
    "Am": None,
    "Sb": 200.0,
    "Ar": None,
    "As": 275.0,
    "Ba": 116.0,
    "Be": 1031.0,
    "Bi": 116.0,
    "B": 1362.0,
    "Cd": 221.0,
    "Ca": 230.0,
    "C": 1550.0,
    "Ce": 138.0,
    "Cs": 43.0,
    "Cr": 424.0,
    "Co": 386.0,
    "Cu": 310.0,
    "Cm": None,
    "Dy": 158.0,
    "Er": 163.0,
    "Eu": None,
    "Fr": None,
    "Gd": 155.0,
    "Ga": 240.0,
    "Ge": 403.0,
    "Au": 178.0,
    "Hf": 213.0,
    "Ho": 161.0,
    "H": None,
    "In": 129.0,
    "I": None,
    "Ir": 228.0,
    "Fe": 373.0,
    "Kr": None,
    "La": 135.0,
    "Pb": 87.0,
    "Li": 448.0,
    "Lu": 116.0,
    "Mg": 330.0,
    "Ma": 363.0,
    "Hg": 92.0,
    "Mo": 377.0,
    "Nd": 148.0,
    "Ne": None,
    "Np": 163.0,
    "Ni": 345.0,
    "Nb": 260.0,
    "Os": 400.0,
    "Pa": 275.0,
    "P": 576.0,
    "Pt": 225.0,
    "Pu": 176.0,
    "K": 100.0,
    "Pr": 138.0,
    "Pa": 262.0,
    "Ra": None,
    "Re": 275.0,
    "Rh": 350.0,
    "Rb": 59.0,
    "Ru": 415.0,
    "Sm": 184.0,
    "Sc": 476.0,
    "Se": None,
    "Si": 692.0,
    "Ag": 221.0,
    "Na": 155.0,
    "Sr": 148.0,
    "S": 527.0,
    "Ta": 225.0,
    "Tc": 422.0,
    "Te": None,
    "Tb": 158.0,
    "Tl": 96.0,
    "Th": 100.0,
    "Tm": 167.0,
    "Sn": 254.0,
    "Ti": 380.0,
    "W": 312.0,
    "U": 300.0,
    "Va": 390.0,
    "Xe": None,
    "Yb": None,
    "Y": 214.0,
    "Zn": 237.0,
    "Zr": 250.0,
}
