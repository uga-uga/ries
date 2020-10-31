import numpy as np

class PhotoabsorptionCrossSectionIsotope:
    def __init__(self, isotope, ResonanceModel,
        resonance_model_parameters):
        self.isotope = isotope
        self.ResonanceModel = ResonanceModel
        self.resonance_model_parameters = resonance_model_parameters
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
                resonances[self.isotope.excited_states[state].J_pi] = self.ResonanceModel(
                    self.isotope.ground_state,
                    self.isotope.excited_states[state],
                    *self.resonance_model_parameters
                )

        return resonances

class PhotoabsorptionCrossSectionElement:
    def __init__(self, element, resonance_model,
        resonance_model_parameters):
        self.element = element
        self.resonance_model = resonance_model
        self.resonance_model_parameters = resonance_model_parameters
        self.photoabsorption_cross_sections = self.get_photoabsorption_cross_sections()

    def __call__(self, energy):
        cross_section = 0.
        for isotope in self.element.isotopes:
            cross_section += self.element.abundances[isotope]*self.photoabsorption_cross_sections[isotope](energy)

        return cross_section

    def get_photoabsorption_cross_sections(self):
        photoabsorption_cross_sections = {}

        for isotope in self.element.isotopes:
            photoabsorption_cross_sections[self.element.isotopes[isotope].AX] = PhotoabsorptionCrossSectionIsotope(
                self.element.isotopes[isotope],
                self.resonance_model[isotope],
                self.resonance_model_parameters[isotope]
            )

        return photoabsorption_cross_sections