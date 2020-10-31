from scipy.stats import cauchy

from ries.cross_section_model import CrossSectionModel

class BreitWigner(CrossSectionModel):
    def __init__(self, initial_state, intermediate_state,
        final_state=None):
        CrossSectionModel.__init__(self, initial_state, intermediate_state, final_state)

        self.probability_distribution = cauchy
        self.probability_distribution_parameters = (
            self.resonance_energy,
            0.5*self.intermediate_state.width
        )