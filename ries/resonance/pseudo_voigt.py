import numpy as np
from scipy.stats import cauchy, norm

from ries.resonance.gauss import Gauss

class PseudoVoigtDistribution:
    def __init__(self):
        pass

    def pdf(self, energy, resonance_energy, eta, gamma_G, gamma_L):
        return(
            (1.-eta)*norm.pdf(
                energy,
                loc=resonance_energy,
                scale=gamma_G/np.sqrt(2.)
            )
            +eta*cauchy.pdf(
                energy,
                loc=resonance_energy,
                scale=0.5*gamma_L
            )
        )

    def ppf(self, quantile, resonance_energy, eta, gamma_G, gamma_L):
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

class PseudoVoigt(Gauss):
    def __init__(self, initial_state, intermediate_state, amu, effective_temperature, final_state=None):
        Gauss.__init__(self, initial_state, intermediate_state, amu, effective_temperature, final_state)

        self.Gamma_L = self.intermediate_state.width
        self.Gamma_L_squared = self.Gamma_L*self.Gamma_L
        self.Gamma_L_cubed = self.Gamma_L_squared*self.Gamma_L
        self.Gamma_L_squared_squared = self.Gamma_L_squared*self.Gamma_L_squared
        self.Gamma_G = 2.*np.sqrt(np.log(2.))*self.doppler_width
        self.Gamma_G_squared = self.Gamma_G*self.Gamma_G
        self.Gamma_G_squared_squared = self.Gamma_G_squared*self.Gamma_G_squared
        
        self.Gamma = self.get_Gamma()
        self.eta = self.get_eta()
        self.gamma_G = self.Gamma/(2.*np.sqrt(np.log(2.)))
        self.gamma_L = 0.5*self.Gamma

        self.probability_distribution = PseudoVoigtDistribution()
        self.probability_distribution_parameters = (
            self.resonance_energy,
            self.eta,
            self.gamma_G,
            self.gamma_L
        )

    def get_Gamma(self):
        return (
            self.Gamma_G_squared_squared*self.Gamma_G
            +2.69269*self.Gamma_G_squared_squared*self.Gamma_L
            +2.42843*self.Gamma_G_squared*self.Gamma_G*self.Gamma_L_squared
            +4.47163*self.Gamma_G_squared*self.Gamma_L_cubed
            +0.07842*self.Gamma_G*self.Gamma_L_squared_squared
            +self.Gamma_L_squared_squared*self.Gamma_L
        )**0.2
    
    def get_eta(self):
        inverse_Gamma = 1./self.Gamma
        inverse_Gamma_squared = inverse_Gamma*inverse_Gamma
        return (
            1.36603*self.Gamma_L*inverse_Gamma
            -0.47719*self.Gamma_L_squared*inverse_Gamma_squared
            +0.11116*self.Gamma_L_cubed*inverse_Gamma_squared*inverse_Gamma
        )