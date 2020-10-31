from .B11 import B11

def test_state():
    assert B11.ground_state.excitation_energy == 0.
    assert B11.excited_states['3/2^-_2'].width == 1.97e-6
    assert B11.excited_states['3/2^-_2'].partial_widths['3/2^-_1'] == 0.856*1.97e-6