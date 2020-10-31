from scipy.constants import physical_constants

class Isotope:
    def __init__(self, AX, amu, 
        ground_state=None, excited_states=[]):
        self.AX = AX
        self.amu = amu
        self.ground_state = ground_state
        self.excited_states = excited_states