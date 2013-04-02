import logging


# Adding logging support
logger = logging.getLogger(__name__)


def run_radial1d(radial1d_model, save_history=None):
    for i in xrange(radial1d_model.iterations - 1):
        logger.info('At run %d of %d', i + 1, radial1d_model.iterations)
        radial1d_model.simulate()        #spec_nu_flux = np.histogram(out_nu, weights=out_energy, bins=radial1d_model.spec_virt_nu)
        if save_history is not None:
            save_history.store_all(radial1d_model, i)

    #Finished second to last loop running one more time
    logger.info('Doing last run')
    if radial1d_model.tardis_config.last_no_of_packets is not None:
        radial1d_model.current_no_of_packets = radial1d_model.tardis_config.last_no_of_packets

    radial1d_model.simulate(enable_virtual=True, update_radiation_field=False)

    if save_history is not None:
        save_history.store_all(radial1d_model, radial1d_model.iterations - 1)
        save_history.finalize()










