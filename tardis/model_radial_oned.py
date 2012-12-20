# building of radial_oned_model

import numpy as np
import plasma, atomic
import logging
import config_reader

logger = logging.getLogger(__name__)


class Radial1DModel(object):
    """
        Class to hold the states of the individual shells (the state of the plasma (as a `~plasma.BasePlasma`-object or one of its subclasses),
        , the plasma parameters (e.g. temperature, dilution factor), the dimensions of the shell).


        Parameters
        ----------

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
    def from_config(cls, fname, atom_data, plasma_type='lte'):
        tardis_config = config_reader.read_config(fname)

        return cls(tardis_config.velocities, tardis_config.densities, tardis_config.abundances,
            tardis_config.time_explosion, atom_data, plasma_type=tardis_config.plasma_type)


    def __init__(self, velocities, densities, abundances, time_explosion, atom_data, ws=None, plasma_type='lte',
                 initial_t_rad=10000):
        no_of_shells = len(velocities) - 1

        logger.info('Assuming %d shells' % no_of_shells)


        #setting time_explosion
        self.time_explosion = time_explosion

        #initializing velocities and radii
        self.v_inner = velocities[:-1]
        self.v_outer = velocities[1:]

        self.r_inner = self.v_inner * time_explosion
        self.r_outer = self.v_outer * time_explosion
        self.r_middle = 0.5 * (self.r_inner + self.r_outer)

        self.volumes = (4. / 3) * np.pi * (self.r_outer ** 3 - self.r_inner ** 3)
        #initializing densities
        assert len(densities) == no_of_shells

        self.densities_middle = densities


        #initializing abundances
        if isinstance(abundances, dict):
            logger.info('Only one abundance - assuming uniform abundance stratification')
            self.abundances = [abundances] * no_of_shells
        else:
            assert len(abundances) == no_of_shells
            self.abundances = abundances



        #Selecting plasma class
        self.plasma_type = plasma_type
        if plasma_type == 'lte':
            self.plasma_class = plasma.LTEPlasma
            if ws is not None:
                raise ValueError(
                    "the dilution factor W ('ws') can only be specified when selecting plasma_type='nebular'")
        elif plasma_type == 'nebular':
            self.plasma_class = plasma.NebularPlasma
            if not isinstance(atom_data, atomic.NebularAtomData):
                raise ValueError("Requiring Nebular Atom data for 'nebular' plasma_type")
        else:
            raise ValueError("Currently this model only supports 'lte' or 'nebular'")



        #setting atom data and checking for consistency
        self.atom_data = atom_data

        #setting dilution factors
        if ws is None:
            self.ws = 0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle ** 2))
        else:
            self.ws = np.array([(0.5 * (1 - np.sqrt(1 - self.r_inner[0] ** 2 / self.r_middle[i] ** 2))) if w < 0\
                                else w for i, w in enumerate(ws)])

        #initializing temperatures

        if np.isscalar(initial_t_rad):
            self.t_rads = [initial_t_rad] * no_of_shells
        else:
            assert len(initial_t_rad) == no_of_shells
            self.t_rads = np.array(initial_t_rad, dtype=np.float64)

        self.initialize_plasmas(plasma_type)


    @property
    def electron_density(self):
        return np.array([plasma.electron_density for plasma in self.plasmas])


    @property
    def tau_sobolevs(self):
        tau_sobolevs = []
        for plasma in self.plasmas:
            tau_sobolevs.append(plasma.lines_data[['nu', 'tau_sobolev']].values)
        return tau_sobolevs

    def initialize_plasmas(self, plasma_type):
        self.plasmas = []

        if plasma_type == 'lte':
            for (current_abundances, current_t_rad, current_density) in\
            zip(self.abundances, self.t_rads, self.densities_middle):
                current_plasma = self.plasma_class(current_abundances, current_density, self.atom_data)
                current_plasma.update_radiationfield(current_t_rad)
                current_plasma.update_radiationfield(current_t_rad)
                current_plasma.calculate_tau_sobolev(self.time_explosion)
                self.plasmas.append(current_plasma)