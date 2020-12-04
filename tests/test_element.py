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
from scipy.constants import physical_constants

from ries.constituents.element import natural_elements, X

def test_element():
    assert natural_elements['Pb'].Z == 82
    assert len(natural_elements['Pb'].isotopes) == 4
    assert natural_elements['Pb'].isotopes['208Pb'].amu == 207.9766525
    assert natural_elements['Pb'].abundances['208Pb'] == 0.524
    assert natural_elements['Pb'].amu == (203.9730440*0.014 + 205.9744657*0.241 + 206.9758973*0.221 + 207.9766525*0.524)