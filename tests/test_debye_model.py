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

from ries.resonance.debye_model import effective_temperature_debye_approximation


def test_debye_model():
    T = 1.0
    T_D = 100.0
    T_eff_low_temperature_limit = 3.0 / 8.0 * T_D
    assert np.isclose(
        effective_temperature_debye_approximation(T, T_D),
        T_eff_low_temperature_limit,
        rtol=1e-3,
    )

    T = 1e3
    T_D = 100.0
    T_eff_high_temperature_limit = T
    assert np.isclose(
        effective_temperature_debye_approximation(T, T_D),
        T_eff_high_temperature_limit,
        rtol=1e-3,
    )
