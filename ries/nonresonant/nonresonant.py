"""
Abstract class for a nonresonant cross section.

At the moment, this class exists as a counterpiece to the `ries.resonance.Resonance` class, so \
that nonresonant processes can inherit from it.
"""

from ries.cross_section import CrossSection

class Nonresonant(CrossSection):
    """Nonresonant cross section"""

    def equidistant_probability_grid(self, limits, n_points):
        """Create a grid with equal reaction probabilities per interval
        
Since this is still a general class, no assumption is made about the PDF of the cross section, \
i.e. a constant-energy grid is returned.

See also `CrossSection.equidistant_probability_grid()`.
        """
        return self.equidistant_energy_grid(limits, n_points)