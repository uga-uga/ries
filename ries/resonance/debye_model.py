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
    return 3.*(T/T_D)**3*T*quad(
        lambda t: t**3*(1./(np.exp(t)-1.)+0.5),
        0., T_D/T
    )[0]

room_temperature_T_D = {
    'Al':  390.,
    'Am':  None,
    'Sb':  200.,
    'Ar':  None,
    'As':  275.,
    'Ba':  116.,
    'Be': 1031.,
    'Bi':  116.,
    'B' : 1362.,
    'Cd':  221.,
    'Ca':  230.,
    'C' : 1550.,
    'Ce':  138.,
    'Cs':   43.,
    'Cr':  424.,
    'Co':  386.,
    'Cu':  310.,
    'Cm':  None,
    'Dy':  158.,
    'Er':  163.,
    'Eu':  None,
    'Fr':  None,
    'Gd':  155.,
    'Ga':  240.,
    'Ge':  403.,
    'Au':  178.,
    'Hf':  213.,
    'Ho':  161.,
    'H' :  None,
    'In':  129.,
    'I' :  None,
    'Ir':  228.,
    'Fe':  373.,
    'Kr':  None,
    'La':  135.,
    'Pb':   87.,
    'Li':  448.,
    'Lu':  116.,
    'Mg':  330.,
    'Ma':  363.,
    'Hg':   92.,
    'Mo':  377.,
    'Nd':  148.,
    'Ne':  None,
    'Np':  163.,
    'Ni':  345.,
    'Nb':  260.,
    'Os':  400.,
    'Pa':  275.,
    'P' :  576.,
    'Pt':  225.,
    'Pu':  176.,
    'K' :  100.,
    'Pr':  138.,
    'Pa':  262.,
    'Ra':  None,
    'Re':  275.,
    'Rh':  350.,
    'Rb':   59.,
    'Ru':  415.,
    'Sm':  184.,
    'Sc':  476.,
    'Se':  None,
    'Si':  692.,
    'Ag':  221.,
    'Na':  155.,
    'Sr':  148.,
    'S' :  527.,
    'Ta':  225.,
    'Tc':  422.,
    'Te':  None,
    'Tb':  158.,
    'Tl':   96.,
    'Th':  100.,
    'Tm':  167.,
    'Sn':  254.,
    'Ti':  380.,
    'W' :  312.,
    'U' :  300.,
    'Va':  390.,
    'Xe':  None,
    'Yb':  None,
    'Y' :  214.,
    'Zn':  237.,
    'Zr':  250.,
}