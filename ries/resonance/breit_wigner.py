from scipy.stats import cauchy

from ries.resonance.resonance_model import ResonanceModel

class BreitWigner(ResonanceModel):
    def __init__(self, initial_state, intermediate_state,
        final_state=None):
        ResonanceModel.__init__(self, initial_state, intermediate_state, final_state)

        self.probability_distribution = cauchy
        self.probability_distribution_parameters = (
            self.resonance_energy,
            0.5*self.intermediate_state.width
        )