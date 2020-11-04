from warnings import warn

import numpy as np

class CrossSection:
    def __add__(self, other):
        if isinstance(other, (int, float)):
            return CrossSectionWeightedSum([self, ConstantCrossSection(other)])
        return CrossSectionWeightedSum([self, other])

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        return CrossSectionWeightedSum([self], [other])

    def __rmul__(self, other):
        return self.__mul__(other)

    def equidistant_energy_grid(self, limits, n_points):
        return np.linspace(limits[0], limits[1], n_points)

class ConstantCrossSection(CrossSection):
    def __init__(self, constant):
        self.constant = constant

    def __call__(self, energy):
        return self.constant

class CrossSectionWeightedSum:
    def __init__(self, reactions=[], scale_factors=None):
        self.reactions = reactions
        if scale_factors is None:
            scale_factors = [1.]*len(reactions)
        self.scale_factors = scale_factors

    def __add__(self, other):
        reactions = []
        scale_factors = []
        for i, reaction in enumerate(self.reactions):
            reactions.append(reaction)
            scale_factors.append(self.scale_factors[i])

        if isinstance(other, (int, float)):
            reactions.append(ConstantCrossSection(other))
            scale_factors.append(1.)
        elif isinstance(other, CrossSection):
            reactions.append(other)
            scale_factors.append(1.)
        else:
            for i, reaction in enumerate(other.reactions):
                reactions.append(reaction)
                scale_factors.append(other.scale_factors[i])

        return CrossSectionWeightedSum(reactions, scale_factors)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        reactions = []
        scale_factors = []
        for i, reaction in enumerate(self.reactions):
            reactions.append(reaction)
            scale_factors.append(other*self.scale_factors[i])
        return CrossSectionWeightedSum(reactions, scale_factors)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __call__(self, energy):
        cross_section = 0.

        for i, reaction in enumerate(self.reactions):
            cross_section += self.scale_factors[i]*reaction(energy)

        return cross_section