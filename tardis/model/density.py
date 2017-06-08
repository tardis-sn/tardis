import numpy as np
import os
from tardis.util import quantity_linspace
from tardis.io.util import HDFReaderWriter

class HomologousDensity(HDFReaderWriter, object):
    """A class that holds an initial density and time

    Parameters
    ----------
    density_0 : astropy.units.Quantity
    time_0 : astropy.units.Quantity

    """
    hdf_properties = ['density_0', 'time_0']
    quantity_attrs = {'density_0':'g/cm^3', 'time_0':'s'}

    def __init__(self, density_0, time_0):
        self.density_0 = density_0
        self.time_0 = time_0

    def calculate_density_at_time_of_simulation(self, time_explosion):
        """
        Scale `density_0` from an `time_0` to the `time_explosion` by ^-3

        Parameters
        ----------
        time_explosion : astropy.units.Quantity

        Returns
        -------
        astropy.units.Quantity
        """
        return calculate_density_after_time(self.density_0, self.time_0,
                                            time_explosion).cgs

    @classmethod
    def from_config(cls, config):
        """
        Create a new HomologousDensity instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration

        Returns
        -------
        HomologousDensity

        """
        d_conf = config.model.structure.density
        velocity = quantity_linspace(config.model.structure.velocity.start,
                                     config.model.structure.velocity.stop,
                                     config.model.structure.velocity.num + 1).cgs

        adjusted_velocity = velocity.insert(0, 0)
        v_middle = (adjusted_velocity[1:] * 0.5 +
                    adjusted_velocity[:-1] * 0.5)
        no_of_shells = len(adjusted_velocity) - 1
        time_explosion = config.supernova.time_explosion.cgs

        if d_conf.type == 'branch85_w7':
            density_0 = calculate_power_law_density(v_middle, d_conf.w7_v_0,
                                                    d_conf.w7_rho_0, -7)
            time_0 = d_conf.w7_time_0
        elif d_conf.type == 'uniform':
            density_0 = (d_conf.value.to('g cm^-3') *
                         np.ones(no_of_shells))
            time_0 = d_conf.get('time_0', time_explosion)
        elif d_conf.type == 'power_law':
            density_0 = calculate_power_law_density(v_middle, d_conf.v_0,
                                                    d_conf.rho_0,
                                                    d_conf.exponent)
            time_0 = d_conf.get('time_0', time_explosion)
        elif d_conf.type == 'exponential':
            density_0 = calculate_exponential_density(v_middle, d_conf.v_0,
                                                      d_conf.rho_0)
            time_0 = d_conf.get('time_0', time_explosion)
        else:
            raise ValueError("Unrecognized density type "
                             "'{}'".format(d_conf.type))
        return cls(density_0, time_0)

    def to_hdf(self, path_or_buf, path=''):
        """
        Store the HomologousDensity to an HDF structure.
        Parameters
        ----------
        path_or_buf
            Path or buffer to the HDF store
        path : str
            Path inside the HDF store to store the HomologousDensity
        Returns
        -------
        None
        """
        homologous_density_path = os.path.join(path, 'homologous_density')
        self.to_hdf_util(path_or_buf, homologous_density_path, {name: getattr(self, name) for name
                                                      in self.hdf_properties})

    @classmethod
    def from_hdf(cls, file_path, path=''):
        """
        This function returns a HomologousDensity object 
        from given HDF5 File.
        Parameters
        ----------
        path : 'str'
            Path to transverse in hdf file
        file_path : 'str'
            Path of Simulation generated HDF file 
        Returns
        -------
        model : `~HomologousDensity`
        """


        buff_path = path + '/homologous_density'
        data = cls.from_hdf_util(file_path, buff_path)
        return cls(data['density_0'], data['time_0'])


def calculate_power_law_density(velocities, velocity_0, rho_0, exponent):
    """

    This function computes a descret exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\left( \\frac{v}{v_0} \\right)^n`

    Parameters
    ----------

    velocities : ~astropy.Quantity
        Array like velocity profile
    velocity_0 : ~astropy.Quantity
        reference velocity
    rho_0 : ~astropy.Quantity
        reference density
    exponent : ~float
        exponent used in the powerlaw

    Returns
    -------

    densities : ~astropy.Quantity

    """
    densities = rho_0 * np.power((velocities / velocity_0), exponent)
    return densities


def calculate_exponential_density(velocities, velocity_0, rho_0):
    """
    This function computes the exponential density profile.
    :math:`\\rho = \\rho_0 \\times \\exp \\left( -\\frac{v}{v_0} \\right)`

    Parameters
    ----------

    velocities : ~astropy.Quantity
        Array like velocity profile
    velocity_0 : ~astropy.Quantity
        reference velocity
    rho_0 : ~astropy.Quantity
        reference density

    Returns
    -------

    densities : ~astropy.Quantity

    """
    densities = rho_0 * np.exp(-(velocities / velocity_0))
    return densities


def calculate_density_after_time(densities, time_0, time_explosion):
    """
    scale the density from an initial time of the model to the
    time of the explosion by ^-3

    Parameters:
    -----------

    densities: ~astropy.units.Quantity
        densities

    time_0: ~astropy.units.Quantity
        time of the model

    time_explosion: ~astropy.units.Quantity
        time to be scaled to

    Returns:
    --------

    scaled_density
    """

    return densities * (time_explosion / time_0) ** -3
