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

from ries.constituents.element import natural_elements
from ries.nonresonant.xrmac import xrmac

def test_xrmac():
    xrmac_test = 7.102e-2*1e26*natural_elements['Pb'].amu*physical_constants['atomic mass constant'][0]*1e3
    assert np.isclose(xrmac['Pb'](1.), xrmac_test, 1e-5)