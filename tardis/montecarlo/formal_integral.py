from astropy import units as u

from tardis.montecarlo.montecarlo import formal_integral
from tardis.montecarlo.spectrum import TARDISSpectrum


class FormalIntegrator(object):

    def __init__(self, model, plasma, runner):
        self.model = model
        self.plasma = plasma
        self.runner = runner

    def calculate_spectrum(self, frequency, N=1000):
        frequency = frequency.to('Hz', u.spectral())
        luminosity = formal_integral(
                self,
                frequency,
                N)
        return TARDISSpectrum(
                frequency,
                luminosity
                )

    def make_source_function(self):
        """
        Calculates the source function using the line absorption rate estimator `Edotlu_estimator`

        Formally it calculates the expression ( 1 - exp(-tau_ul) ) S_ul but this product is what we need later,
        so there is no need to factor out the source function explicitly.

        Parameters
        ----------
        model : tardis.model.Radial1DModel

        Returns
        -------
        Numpy array containing ( 1 - exp(-tau_ul) ) S_ul ordered by wavelength of the transition u -> l
        """
        model = self.model
        plasma = self.plasma
        runner = self.runner
        atomic_data = self.plasma.atomic_data

        Edotlu_norm_factor = (1 / (runner.time_of_simulation * model.volume))
        exptau = 1 - np.exp(- plasma.tau_sobolevs)
        Edotlu = Edotlu_norm_factor * exptau * runner.Edotlu_estimator

        upper_level_index = atomic_data.lines.set_index(['atomic_number', 'ion_number', 'level_number_upper']).index.copy()
        e_dot_lu          = pd.DataFrame(Edotlu, index=upper_level_index)
        e_dot_u           = e_dot_lu.groupby(level=[0, 1, 2]).sum()
        e_dot_u.index.names = ['atomic_number', 'ion_number', 'source_level_number'] # To make the q_ul e_dot_u product work, could be cleaner
        transitions       = atomic_data.macro_atom_data[atomic_data.macro_atom_data.transition_type == -1].copy()
        transitions_index = transitions.set_index(['atomic_number', 'ion_number', 'source_level_number']).index.copy()
        tmp  = plasma.transition_probabilities[(atomic_data.macro_atom_data.transition_type == -1).values]
        q_ul = tmp.set_index(transitions_index)
        t    = model.time_explosion.value
        wave = atomic_data.lines.wavelength_cm[transitions.transition_line_id].values.reshape(-1,1)
        att_S_ul =  ( wave * (q_ul * e_dot_u) * t  / (4*np.pi) )

        result = pd.DataFrame(att_S_ul.as_matrix(), index=transitions.transition_line_id.values)
        return result.ix[atomic_data.lines.index.values].as_matrix(), e_dot_u
