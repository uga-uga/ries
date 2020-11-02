import numpy as np

class PhotoabsorptionCrossSection:
    def __init__(self, isotope, Resonance,
        resonance_parameters):
        self.isotope = isotope
        self.Resonance = Resonance
        self.resonance_parameters = resonance_parameters
        self.resonances = self.get_resonances()

    def __call__(self, energy):
        cross_section = 0.
        for resonance in self.resonances:
            cross_section += self.resonances[resonance](energy)

        return cross_section

    def get_resonances(self):
        resonances = {}

        for state in self.isotope.excited_states:
            if self.isotope.ground_state.J_pi in self.isotope.excited_states[state].partial_widths:
                resonances[self.isotope.excited_states[state].J_pi] = self.Resonance(
                    self.isotope.ground_state,
                    self.isotope.excited_states[state],
                    *self.resonance_parameters
                )

        return resonances