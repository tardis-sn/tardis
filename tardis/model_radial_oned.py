# building of radial_oned_model

import numpy as np
import plasma, atomic, packet_source
import logging
import config_reader
import pandas as pd

logger = logging.getLogger(__name__)


class Radial1DModel(object):
    """
        Class to hold the states of the individual shells (the state of the plasma (as a `~plasma.BasePlasma`-object or one of its subclasses),
        , the plasma parameters (e.g. temperature, dilution factor), the dimensions of the shell).


        Parameters
        ----------

        tardis_configuration : `tardis.config_reader.TardisConfiguration`

        velocities : `np.ndarray`
            an array with n+1 (for n shells) velocities (in cm/s) for each of the boundaries (velocities[0] describing
            the inner boundary and velocities[-1] the outer boundary

        densities : `np.ndarray`
            an array with n densities - being the density mid-shell (assumed for the whole shell)

        abundances : `list` or `dict`
            a dictionary for uniform abundances throughout all shells, e.g. dict(Fe=0.5, Si=0.5)
            For a different abundance for each shell list of abundance dictionaries.


        time_explosion : `float`
            time since explosion in seconds

        atom_data : `~tardis.atom_data.AtomData` class or subclass
            Containing the atom data needed for the plasma calculations

        ws : `None` or `list`-like
            ws can only be specified for plasma_type 'nebular'. If `None` is specified at first initialization the class
            calculates an initial geometric dilution factor. When giving a list positive values will be accepted, whereas
            negative values trigger the usage of the geometric calculation

        plasma_type : `str`
            plasma type currently supports 'lte' (using `tardis.plasma.LTEPlasma`)
            or 'nebular' (using `tardis.plasma.NebularPlasma`)

        initial_t_rad : `float`-like or `list`-like
            initial radiative temperature for each shell, if a scalar is specified it initializes with a uniform
            temperature for all shells




    """

    @classmethod
    def from_config_file(cls, fname, atom_data):
        tardis_config = config_reader.read_config(fname)
        model_obj = cls(tardis_config, atom_data)

        return model_obj


    def __init__(self, configuration_object, atom_data):
        #final preparation for configuration object
        configuration_object.final_preparation()


        #final preparation for atom_data object - currently building data
        atom_data.prepare_atom_data(configuration_object)

        self.atom_data = atom_data

        self.packet_src = packet_source.SimplePacketSource.from_wavelength(configuration_object.spectrum_start,
            configuration_object.spectrum_end)

        self.packet_src.create_packets(configuration_object.spectrum_packets, configuration_object.t_rad_inner)

        self.no_of_shells = configuration_object.no_of_shells

        logger.info('Assuming %d shells' % self.no_of_shells)


        #setting time_explosion
        self.time_explosion = configuration_object.time_explosion

        #initializing velocities and radii
        self.v_inner = configuration_object.velocities[:-1]
        self.v_outer = configuration_object.velocities[1:]

        self.r_inner = self.v_inner * self.time_explosion
        self.r_outer = self.v_outer * self.time_explosion
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)

        self.volumes = (4. / 3) * np.pi * (self.r_outer ** 3 - self.r_inner ** 3)
        #initializing densities
        assert len(configuration_object.densities) == self.no_of_shells

        self.densities_middle = configuration_object.densities




        #Selecting plasma class
        self.plasma_type = configuration_object.plasma_type
        if self.plasma_type == 'lte':
            self.plasma_class = plasma.LTEPlasma
            if configuration_object.ws is not None:
                raise ValueError(
                    "the dilution factor W ('ws') can only be specified when selecting plasma_type='nebular'")

        elif self.plasma_type == 'nebular':
            self.plasma_class = plasma.NebularPlasma
            if not self.atom_data.has_zeta_data:
                raise ValueError("Requiring Recombination coefficients Zeta for 'nebular' plasma_type")
        else:
            raise ValueError("Currently this model only supports 'lte' or 'nebular'")

        if configuration_object.line_interaction_type == 'scatter':
            self.line_interaction_id = 0

        #setting atom data and checking for consistency
        self.atom_data = atom_data

        #initializing abundances
        self.abundances = configuration_object.abundances
        self.number_densities = []
        for abundance, density in zip(self.abundances, self.densities_middle):
            self.number_densities.append(calculate_atom_number_densities(self.atom_data, abundance, density))



        #setting dilution factors
        if configuration_object.ws is None:
            self.ws = 0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle ** 2))
        else:
            self.ws = np.array([(0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle[i] ** 2))) if w < 0\
                                else w for i, w in enumerate(configuration_object.ws)])

        #initializing temperatures

        if np.isscalar(configuration_object.initial_t_rad):
            self.t_rads = [configuration_object.initial_t_rad] * self.no_of_shells
        else:
            assert len(configuration_object.initial_t_rad) == self.no_of_shells
            self.t_rads = np.array(configuration_object.initial_t_rad, dtype=np.float64)

        self.initialize_plasmas(configuration_object.plasma_type)


    @property
    def electron_density(self):
        return np.array([plasma.electron_density for plasma in self.plasmas])


    #    @property
    #    def tau_sobolevs(self):
    #        tau_sobolevs = []
    #        for plasma in self.plasmas:
    #            tau_sobolevs.append(plasma.lines_data[['nu', 'tau_sobolev']].values)
    #        return tau_sobolevs

    def initialize_plasmas(self, plasma_type):
        self.plasmas = []
        self.tau_sobolevs = []
        self.line_list_nu = self.atom_data.lines['nu']

        if plasma_type == 'lte':
            for (current_abundances, current_t_rad) in\
            zip(self.number_densities, self.t_rads):
                current_plasma = self.plasma_class(current_abundances, self.atom_data)
                current_plasma.update_radiationfield(current_t_rad)
                self.tau_sobolevs.append(current_plasma.calculate_tau_sobolev(self.time_explosion))
                self.plasmas.append(current_plasma)


        elif plasma_type == 'nebular':
            for (current_abundances, current_t_rad, w) in\
            zip(self.number_densities, self.t_rads, self.ws):
                current_plasma = self.plasma_class(current_abundances, self.atom_data)
                current_plasma.update_radiationfield(current_t_rad, w)
                self.tau_sobolevs.append(current_plasma.calculate_tau_sobolev(self.time_explosion))
                self.plasmas.append(current_plasma)

            self.tau_sobolevs = np.array(self.tau_sobolevs, dtype=float)

            # update plasmas

    def update_plasmas(self):
        self.tau_sobolevs = []
        for current_plasma in self.plasmas:
            current_plasma.update_radiationfield(20000)
            self.tau_sobolevs.append(current_plasma.calculate_tau_sobolev(self.time_explosion))


def calculate_atom_number_densities(atom_data, abundances, density):
    """
    Calculates the atom number density, using the following formula, where Z is the atomic number
    and X is the abundance fraction

    .. math::
        N_{Z} = \\frac{\\rho_\\textrm{total}\\times \\textrm{X}_\\textrm{Z}}{m_\\textrm{Z}}

    """

    #Converting abundances



    abundance_fractions = pd.Series(abundances.values(),
        index=pd.Index([atom_data.symbol2atomic_number[item] for item in abundances.keys()],
            dtype=np.int, name='atomic_number'), name='abundance_fraction')



    #Normalizing Abundances

    abundance_sum = abundance_fractions.sum()

    if abs(abundance_sum - 1) > 1e-5:
        logger.warn('Abundances do not add up to 1 (Sum = %.4f). Renormalizing', (abundance_sum))

    abundance_fractions /= abundance_sum

    number_densities = (abundance_fractions * density) /\
                       atom_data.atom_data.ix[abundance_fractions.index]['mass']

    return pd.DataFrame({'abundance_fraction': abundance_fractions, 'number_density': number_densities})

