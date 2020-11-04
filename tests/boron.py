from ries.constituents.element import Element, natural_elements
from ries.constituents.isotope import Isotope
from ries.constituents.state import GroundState, State

B10 = Isotope(
    AX='10B', 
    amu=natural_elements['B'].isotopes['10B'].amu,
    ground_state=GroundState('3^+_1', 6, -1),
    excited_states={
        '1^+_1': State(
            J_pi='1^+_1', 
            two_J=2,
            parity=1,
            excitation_energy=0.718380,
            partial_widths={'3^+_1': 0.}
        ),
    }
)

B11 = Isotope(
    AX='11B', 
    amu=natural_elements['B'].isotopes['11B'].amu,
    ground_state=GroundState('3/2^-_1', 3, -1),
    excited_states={
        '1/2^-_1': State(
            J_pi='1/2^-_1', 
            two_J=1,
            parity=-1,
            excitation_energy=2.124693, 
            partial_widths={'3/2^-_1': 0.117e-6}
        ),
        '5/2^-_1': State(
            J_pi='5/2^-_1', 
            two_J=5,
            parity=-1,
            excitation_energy=4.44498, 
            partial_widths={'3/2^-_1': 0.55e-6}
        ),
        '3/2^-_2': State(
            J_pi='3/2^-_2', 
            two_J=3,
            parity=-1,
            excitation_energy=5.02030, 
            partial_widths={'3/2^-_1': 0.856*1.97e-6, '1/2^-_1': 0.144*1.97e-6}
        ),
        '7/2^-_1': State(
            J_pi='7/2^-_1', 
            two_J=7,
            parity=-1,
            excitation_energy=6.74185, 
            partial_widths={'3/2^-_1': 0.70*0.030e-6, '5/2^-_1': 0.30*0.030e-6}
        ),
        '1/2^+_1': State(
            J_pi='1/2^+_1', 
            two_J=1,
            parity=1,
            excitation_energy=6.79180, 
            partial_widths={'3/2^-_1': 0.675*0.39e-6, '1/2^-_1': 0.285*0.39e-6, '3/2^-_2': 0.04*0.39e-6}
        ),
        '5/2^+_1': State(
            J_pi='5/2^+_1', 
            two_J=5,
            parity=1,
            excitation_energy=7.28551, 
            partial_widths={'3/2^-_1': 0.87*1.14e-6, '5/2^-_1': 0.055*1.14e-6, '3/2^-_2': 0.075*1.14e-6}
        ),
        '3/2^+_1': State(
            J_pi='3/2^+_1', 
            two_J=3,
            parity=1,
            excitation_energy=7.97784, 
            partial_widths={'3/2^-_1': 0.462*1.15e-6, '1/2^-_1': 0.532*1.15e-6, '5/2^+_1': 0.0085*1.15e-6}
        ),
        '3/2^-_3': State(
            J_pi='3/2^-_3', 
            two_J=3,
            parity=-1,
            excitation_energy=8.5601, 
            partial_widths={'3/2^-_1': 0.56*1.00e-6, '1/2^-_1': 0.30*1.00e-6, '5/2^-_1': 0.05*1.00e-6, '3/2^-_1': 0.09*1.00e-6}
        ),
        '5/2^-_2': State(
            J_pi='5/2^-_2', 
            two_J=5,
            parity=-1,
            excitation_energy=8.92047, 
            partial_widths={'3/2^-_1': 0.95*4.374e-6, '5/2^-_1': 0.045*4.374e-6}
        ),
    }
)

natural_boron = Element(5, 'B', 
    isotopes={'10B': B10, '11B': B11},
    abundances={'10B': natural_elements['B'].abundances['10B'], '11B': natural_elements['B'].abundances['11B']}
)