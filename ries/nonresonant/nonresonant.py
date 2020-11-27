"""
Abstract class for a nonresonant cross section.

At the moment, this class exists as a counterpiece to the `ries.resonance.Resonance` class, so \
that nonresonant processes can inherit from it.
It provides no methods or attributes.
"""

from ries.cross_section import CrossSection

class Nonresonant(CrossSection):
    """Nonresonant cross section"""
    pass