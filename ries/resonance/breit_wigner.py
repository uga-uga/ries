from scipy.stats import cauchy

from ries.resonance.resonance import Resonance

class BreitWigner(Resonance):
    def __init__(self, initial_state, intermediate_state,
        final_state=None):
        Resonance.__init__(self, initial_state, intermediate_state, final_state)

        self.probability_distribution = cauchy
        self.probability_distribution_parameters = (
            self.resonance_energy,
            0.5*self.intermediate_state.width
        )