from scipy.constants import physical_constants

class Recoil:
    def __init__(self):
        pass

    def __call__(self, energy_difference):
        return energy_difference

class NoRecoil(Recoil):
    def __init__(self):
        Recoil.__init__(self)

class FreeNucleusRecoil(Recoil):
    def __init__(self, amu):
        self.amu = amu

    def __call__(self, energy_difference):
        return (
            energy_difference*(
                1.+energy_difference
                /(2.*self.amu*physical_constants['atomic mass constant energy equivalent in MeV'][0])
            )
        )