import numpy as np
from scipy.special import voigt_profile
from scipy.stats import cauchy, norm

from ries.resonance.gauss import Gauss
from ries.resonance.pseudo_voigt import PseudoVoigt

class VoigtDistribution:
    def __init__(self):
        pass

    def cdf(self, energy, resonance_energy, sigma, gamma, eta, gamma_G, gamma_L):
        return(
            (1.-eta)*norm.cdf(
                energy,
                loc=resonance_energy,
                scale=gamma_G/np.sqrt(2.)
            )
            +eta*cauchy.cdf(
                energy,
                loc=resonance_energy,
                scale=0.5*gamma_L
            )
        )

    def pdf(self, energy, resonance_energy, sigma, gamma, eta, gamma_G, gamma_L):
        return (
            voigt_profile(
                energy-resonance_energy,
                sigma,
                gamma
            )
        )

    def ppf(self, quantile, resonance_energy, sigma, gamma, eta, gamma_G, gamma_L):
        return(
            (1.-eta)*norm.ppf(
                quantile,
                loc=resonance_energy,
                scale=gamma_G/np.sqrt(2.)
            )
            +eta*cauchy.ppf(
                quantile,
                loc=resonance_energy,
                scale=0.5*gamma_L
            )
        )

class Voigt(Gauss):
    def __init__(self, initial_state, intermediate_state, amu, effective_temperature, final_state=None):
        Gauss.__init__(self, initial_state, intermediate_state, amu, effective_temperature, final_state)
        self.pseudo_voigt = PseudoVoigt(initial_state, intermediate_state, amu, effective_temperature, final_state)

        self.probability_distribution = VoigtDistribution()
        self.probability_distribution_parameters = (
            self.resonance_energy,
            self.doppler_width/np.sqrt(2.),
            0.5*self.intermediate_state.width,
            self.pseudo_voigt.eta,
            self.pseudo_voigt.gamma_G,
            self.pseudo_voigt.gamma_L
        )