import numpy as np
from scipy.constants import physical_constants

from ries.constituents.element import natural_elements, X

def test_element():
    assert natural_elements['Pb'].Z == 82
    assert len(natural_elements['Pb'].isotopes) == 4
    assert natural_elements['Pb'].isotopes['208Pb'].amu == 207.9766525
    assert natural_elements['Pb'].abundances['208Pb'] == 0.524
    assert natural_elements['Pb'].amu == (203.9730440*0.014 + 205.9744657*0.241 + 206.9758973*0.221 + 207.9766525*0.524)
    xrmac = 7.102e-2*1e26*natural_elements['Pb'].amu*physical_constants['atomic mass constant'][0]*1e3
    assert np.isclose(natural_elements['Pb'].xrmac(1.), xrmac, rtol=1e-5)