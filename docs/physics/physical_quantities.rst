.. _physical_quantities:

*****************************
Accessing Physical Quantities
*****************************

In order to compute the synthetic spectrum, Tardis must either be told
or must calculate many physical properties of the model. To understand and
test the code it can be important to look at these values. One
easy way to do this is to run Tardis in an interactive mode and then
inspect the model properties.


Runing in interactive Python session
------------------------------

With iPython installed launch a session using

.. code-block:: none

    ipython --pylab

and then run your calculation, which is based on your "my_config.yml" file

.. code-block:: python

    from tardis import run_tardis
    import yaml

    configuration_dict = yaml.load(open('myconfig.yml')
    model = run_tardis(configuration_dict)


If all goes well, the simulation should run as usual. Afterwards, the
information from the simulation will all exist in "radial1d" and
can be examined. Some examples for useful/interesting quantities are
given below (but much more information is available: contact us via 
`tardis-sn-users <http://groups.google.com/forum/#!forum/tardis-sn-users>`_ if you need
further help).

Examples of finding physical quantities
---------------------------------------

For example, two of our important quantities are the parameters of the
radiation field model, :math:`T_{\rm rad}` and :math:`W`. These exist
as Astropy `Quantities <http://astropy.readthedocs.org/en/v0.2.1/_generated/astropy.units.quantity.Quantity.html>`_. Thus

.. code-block:: python

    model.t_rads.cgs

will give you a list of the :math:`T_{\rm rad}`-values for the model zones
in cgs units. To obtain an array of the values (without units) use

.. code-block:: python

    model.t_rads.cgs.value

Similarly, the :math:`W`-values can be accessed using

.. code-block:: python

    model.ws.cgs


Several important quantities that were setup when the model was defined
by the configuration file are located in the "tardis_config"
section. For example, the inner and outer velocity boundaries of the
zones in the model

.. code-block:: python

    model.tardis_config.structure.v_inner.cgs
    model.tardis_config.structure.v_outer.cgs

and the average density in the zones

.. code-block:: python

    model.tardis_config.structure.mean_densities.cgs

Many other interesting quantities are stored in the
"plasma_array". For example the calculated ion populations or level
populations:

.. code-block:: python

    model.plasma_array.ion_populations
    model.plasma_array.level_populations

These are stored as Pandas `DataFrames
<http://pandas.pydata.org/pandas-docs/version/0.13.1/generated/pandas.DataFrame.html>`_.
An index can be supplied to obtain the population in a particular
zone. E.g., for the ion populations of the innermost zone (index = 0)

.. code-block:: python

    model.plasma_array.ion_populations[0]

Ion populations for a particular ionization stage of a particular
element can be accessed by specifying an appropriate tuple :math:`(Z,C)`, which
identifies the element (via atomic number :math:`Z` ) and the charge
(via the ion charge :math:`C` ). Thus, 

.. code-block:: python

    model.plasma_array.ion_populations.ix[(14,1)]

will identify the ion popuations for  Si II (:math:`Z=14, C=1`) in all
the zones. The above examples can be combined to obtain e.g. the Si II
population in the innermost zone

.. code-block:: python

    model.plasma_array.ion_populations[0].ix[(14,1)]

The level populations are stored (and can be accessed) in a similar
way - a third label can be used to pick out a particular atomic
level. E.g., to pull out the population of the ground state (index 0)
of Si II

.. code-block:: python

    model.plasma_array.level_populations.ix[(14,1,0)]

.. note::

    If you prefer to work in SI units, all the astropy Quantities may
    instead by accessed with "xxx.si".

.. note::

    Information that is not stored as astropy Quantities (e.g. the ion
    an level populations used in the example above) are usually stored
    in cgs units (i.e. :math:`{\rm cm}^{-3}` for the populations).
