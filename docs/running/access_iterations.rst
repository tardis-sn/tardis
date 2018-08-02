**********************************************
Access information about individual iterations
**********************************************

Currently we store information about the plasma state and the inner boundary
for each iteration. This is saved in the simulation object and can be accessed
via

.. code-block:: python

  import tardis
  mdl = tardis.run_tardis("tardis_config.yml")
  mdl.iterations_w

in case of the dilution factor. 

Currently, the following properties are available:

.. code-block:: python

   mdl.iterations_w  # dilution factor in each cell
   mdl.iterations_t_rads  # radiation temperature in each cell
   mdl.iterations_electron_densities  # electron density in each cell
   mdl.iterations_t_inner  # inner boundary temperature
