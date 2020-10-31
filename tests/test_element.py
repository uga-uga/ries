from ries.element import natural_elements, X

def test_element():
    assert natural_elements['Pb'].Z == 82
    assert len(natural_elements['Pb'].isotopes) == 4
    assert natural_elements['Pb'].isotopes['208Pb'].amu == 207.9766525
    assert natural_elements['Pb'].abundances['208Pb'] == 0.524
    assert natural_elements['Pb'].amu == (203.9730440*0.014 + 205.9744657*0.241 + 206.9758973*0.221 + 207.9766525*0.524)
    assert natural_elements['Pb'].xrmac(1.) == 7.102e-2