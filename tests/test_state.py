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

from .boron import B11

def test_state():
    assert B11.ground_state.excitation_energy == 0.
    assert B11.excited_states['3/2^-_2'].width == 1.97e-6
    assert B11.excited_states['3/2^-_2'].partial_widths['3/2^-_1'] == 0.856*1.97e-6