from scipy.constants import physical_constants

class GroundState:
    def __init__(self, J_pi, two_J, parity):
        self.J_pi = J_pi
        self.two_J = two_J
        self.parity = parity
        self.excitation_energy = 0.

class State:
    def __init__(self, J_pi, two_J, parity, excitation_energy, partial_widths):
        GroundState.__init__(self, J_pi, two_J, parity)
        self.excitation_energy = excitation_energy
        self.partial_widths = partial_widths
        self.width = sum([self.partial_widths[state] for state in self.partial_widths])