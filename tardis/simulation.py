import logging


# Adding logging support
logger = logging.getLogger(__name__)


def run_radial1d(radial1d_model, save_history=None):
    while radial1d_model.iterations_remaining > 0:
        logger.info('Remaining run %d', radial1d_model.iterations_remaining)
        radial1d_model.simulate()
        if save_history is not None:
            save_history.store(radial1d_model)


    #Finished second to last loop running one more time
    logger.info('Doing last run')
    if radial1d_model.tardis_config.last_no_of_packets is not None:
        radial1d_model.current_no_of_packets = radial1d_model.tardis_config.last_no_of_packets

    radial1d_model.simulate(enable_virtual=True, update_radiation_field=False)

    if save_history is not None:
        save_history.store(radial1d_model)
        save_history.finalize()

    logger.info("Finished in %d iterations", radial1d_model.iterations_executed)










