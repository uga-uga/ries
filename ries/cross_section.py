class CrossSection:
    def __add__(self, other):
        return SumCrossSection([self, other])

class SumCrossSection:
    def __init__(self, reactions=[]):
        self.reactions = reactions

    def __add__(self, other):
        sum_cross_section = SumCrossSection(self.reactions)
        if isinstance(other, CrossSection):
            sum_cross_section.reactions.append(other)
        elif isinstance(other, SumCrossSection):
            reactions = []
            for reaction in self.reactions:
                reactions.append(reaction)
            for reaction in other.reactions:
                reactions.append(reaction)
            sum_cross_section = SumCrossSection(reactions)

        return sum_cross_section

    def __call__(self, energy):
        cross_section = 0.

        for reaction in self.reactions:
            cross_section += reaction(energy)

        return cross_section