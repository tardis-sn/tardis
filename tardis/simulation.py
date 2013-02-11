import constants
import numpy as np
import montecarlo_multizone
import time
import os
import logging
import synspec
import pdb
import pandas as pd


# Adding logging support
logger = logging.getLogger(__name__)


def run_single_radial1d(radial1d_model):
    """
    Run a single 1D radial simulation

    Parameters
    ----------

    radial1d_model : `~tardis.model_radial_oned.Radial1DModel`


    """

    out_nu, out_energy, j_estimators, nubar_estimators = montecarlo_multizone.montecarlo_radial1d(radial1d_model, 0)


def run_radial1d(radial1d_model, save_history=None):
    for i in range(radial1d_model.iterations):
        out_nu, out_energy, j_estimators, nubar_estimators = montecarlo_multizone.montecarlo_radial1d(radial1d_model)
        updated_t_rads, updated_ws = radial1d_model.calculate_updated_radiationfield(nubar_estimators, j_estimators)

        new_trads = 0.5 * (radial1d_model.t_rads + updated_t_rads)
        new_ws = 0.5 * (radial1d_model.ws + updated_ws)

        print pd.DataFrame({'t_rads': radial1d_model.t_rads, 'updated_t_rads': updated_t_rads, 'new_trads': new_trads,
                            'ws': radial1d_model.ws, 'updated_ws': updated_ws, 'new_ws': new_ws})[::5]

        emitted_energy = radial1d_model.emitted_inner_energy * np.sum(out_energy[out_energy >= 0]) / 1.
        absorbed_energy = radial1d_model.emitted_inner_energy * np.sum(out_energy[out_energy < 0]) / 1.

        print "Luminosity emitted = %.5e Luminosity absorbed = %.5e Luminosity requested = %.5e" % (emitted_energy,
                                                                                                    absorbed_energy,
                                                                                                    radial1d_model.luminosity_outer)
        new_t_inner = radial1d_model.t_inner * (emitted_energy / radial1d_model.luminosity_outer) ** -.25
        print "new t_inner = %.2f" % (new_t_inner,)

        radial1d_model.t_inner = 0.5 * (new_t_inner + radial1d_model.t_inner)
        radial1d_model.create_packets()
        radial1d_model.update_plasmas(new_trads, new_ws)
        spec_nu_flux = np.histogram(out_nu, weights=out_energy, bins=radial1d_model.spec_virt_nu)


        #trying out the virtual packets bit
        out_nu, out_energy, j_estimators, nubar_estimators = montecarlo_multizone.montecarlo_radial1d(radial1d_model,
            10)

        #return out_energy


        if save_history is not None:
            save_history.store_all(radial1d_model)

    return spec_nu_flux



