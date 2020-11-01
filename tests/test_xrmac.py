import numpy as np
from scipy.constants import physical_constants

from ries.constituents.element import natural_elements
from ries.nonresonant.xrmac import xrmac

def test_xrmac():
    xrmac_test = 7.102e-2*1e26*natural_elements['Pb'].amu*physical_constants['atomic mass constant'][0]*1e3
    assert np.isclose(xrmac['Pb'](1.), xrmac_test, 1e-5)