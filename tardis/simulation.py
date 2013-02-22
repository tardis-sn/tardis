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
    for i in xrange(radial1d_model.iterations - 1):
        logger.info('At run %d of %d', i + 1, radial1d_model.iterations)
        radial1d_model.create_packets()
        out_nu, out_energy, j_estimators, nubar_estimators = montecarlo_multizone.montecarlo_radial1d(radial1d_model)
        radial1d_model.normalize_j_blues()
        if save_history is not None:
            save_history.store_all(radial1d_model, i)

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

        radial1d_model.update_plasmas(new_trads, new_ws)
        #spec_nu_flux = np.histogram(out_nu, weights=out_energy, bins=radial1d_model.spec_virt_nu)

    #Finished second to last loop running one more time
    logger.info('Doing last run')
    if radial1d_model.tardis_config.last_no_of_packets is not None:
        radial1d_model.create_packets(radial1d_model.tardis_config.last_no_of_packets)
    else:
        radial1d_model.create_packets()
    out_nu, out_energy, j_estimators, nubar_estimators = montecarlo_multizone.montecarlo_radial1d(radial1d_model,
                                                                                                  virtual_packet_flag=radial1d_model.tardis_config.no_of_virtual_packets)
    radial1d_model.normalize_j_blues()
    radial1d_model.calculate_spectrum(out_nu, out_energy, distance=radial1d_model.tardis_config.sn_distance)

    if save_history is not None:
        save_history.store_all(radial1d_model, radial1d_model.iterations - 1)
        save_history.finalize()










