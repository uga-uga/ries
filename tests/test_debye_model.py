import numpy as np

from ries.resonance.debye_model import effective_temperature_debye_approximation

def test_debye_model():
    T = 1.
    T_D = 100.
    T_eff_low_temperature_limit = 3./8.*T_D
    assert np.isclose(effective_temperature_debye_approximation(T, T_D), T_eff_low_temperature_limit, rtol=1e-3)

    T = 1e3
    T_D = 100.
    T_eff_high_temperature_limit = T
    assert np.isclose(effective_temperature_debye_approximation(T, T_D), T_eff_high_temperature_limit, rtol=1e-3)