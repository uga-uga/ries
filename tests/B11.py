from ries.isotope import Isotope
from ries.state import GroundState, State

B11 = Isotope(
    AX='11B', 
    amu=11.009305166,
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
        )
    }
)