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

"""
Abstract class for a nonresonant cross section.

At the moment, this class exists as a counterpiece to the `ries.resonance.Resonance` class, so
that nonresonant processes can inherit from it.
"""

from ries.cross_section import CrossSection

class Nonresonant(CrossSection):
    """Nonresonant cross section"""

    def equidistant_probability_grid(self, limits, n_points):
        """Create a grid with equal reaction probabilities per interval
        
        Since this is still a general class, no assumption is made about the PDF of the cross section,
        i.e. a constant-energy grid is returned.

        See also `CrossSection.equidistant_probability_grid()`.
        """
        return self.equidistant_energy_grid(limits, n_points)
