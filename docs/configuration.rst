Configuration Files
===================

.. literalinclude:: example_configs/example_configuration.ini



Simulation Type
---------------

There are two types of simulations.

1. Single run - will run the simulation once with the specified parameters

This requires the setting of ``single_run_packets keyword`` in the config file

2. Temperature convergence - will run the simulation until the temperature structure converged.

This requires the setting of ``calibration_packets``, ``spectrum_packets`` and ``iterations``


