**************
Running TARDIS
**************

To run TARDIS requires two files. The atomic database (for more info refer to :ref:`atomic-data`) and a
configuration file (more info at :ref:`config-file`).

Currently there is no script that can run TARDIS. However it is very easy to set one up:

.. code-block:: python

    from tardis import config_reader, model_radial_oned, simulation

    tardis_config = config_reader.TARDISConfiguration.from_yaml('myconfig.yml')
    radial1d_mdl = model_radial_oned.Radial1DModel(tardis_config)
    simulation.run_radial1d(radial1d_mdl)

