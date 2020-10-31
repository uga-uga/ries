import numpy as np

from scipy.constants import physical_constants
from scipy.stats import uniform

from ries.recoil import NoRecoil

class CrossSectionModel:
    def __init__(self, initial_state, intermediate_state,
        final_state=None, recoil_correction=NoRecoil()):
        self.initial_state = initial_state
        self.intermediate_state = intermediate_state
        self.final_state = final_state

        self.resonance_energy = recoil_correction(self.intermediate_state.excitation_energy - self.initial_state.excitation_energy)
        self.energy_integrated_cross_section_constant = (np.pi*physical_constants['reduced Planck constant times c in MeV fm'][0])**2
        self.statistical_factor = self.get_statistical_factor()
        self.final_state_branching_ratio = self.get_final_state_branching_ratio()
        self.energy_integrated_cross_section = self.get_energy_integrated_cross_section()

        self.probability_distribution = uniform
        self.probability_distribution_parameters = (self.resonance_energy-0.5, 1.)

    def __call__(self, energy, input_is_absolute_energy=True):
        if not input_is_absolute_energy:
            energy = energy + self.resonance_energy
        return self.energy_integrated_cross_section*self.probability_distribution.pdf(energy, *self.probability_distribution_parameters)

    def coverage_interval(self, coverage):
        return self.probability_distribution.ppf(
            0.5*np.array([1. - coverage, 1.+coverage]),
            *self.probability_distribution_parameters
        )

    def energy_grid(self, coverage, n_points):
        coverage_interval = self.coverage_interval(coverage)
        return np.linspace(coverage_interval[0], coverage_interval[1], n_points)

    def probability_grid(self, coverage, n_points):
        return self.probability_distribution.ppf(
            0.5*np.linspace(1. - coverage, 1.+coverage, n_points),
            *self.probability_distribution_parameters
        )

    def get_energy_integrated_cross_section(self):
        return (
            self.energy_integrated_cross_section_constant
            /((self.resonance_energy)**2)
            *self.statistical_factor
            *self.intermediate_state.partial_widths[self.initial_state.J_pi]
            *self.final_state_branching_ratio
        )

    def get_final_state_branching_ratio(self):
        if self.final_state is None:
            return 1.
        return (
            self.intermediate_state.partial_widths[self.final_state.J_pi]
            /self.intermediate_state.width)

    def get_statistical_factor(self):
        return (
            (self.intermediate_state.two_J+1.)
            /(self.initial_state.two_J+1.)
        )