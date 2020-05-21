Running TARDIS in the commandline
=================================

.. warning::

    This option will be removed in the next versions of TARDIS


After installing TARDIS, just download the configuration file from the
`tardis-setups <https://github.com/tardis-sn/tardis-setups>`_ and the standard
atomic data set from the `tardis-refdata
<https://github.com/tardis-sn/tardis-refdata>`_ repository and run TARDIS.
Assuming you have ``wget``, you could follow the procedure:


.. code-block:: none

    mkdir tardis_example
    cd tardis_example
    wget https://raw.githubusercontent.com/tardis-sn/tardis-setups/master/2014/2014_kerzendorf_sim/appendix_A1/tardis_example.yml
    wget https://github.com/tardis-sn/tardis-refdata/raw/master/atom_data/kurucz_cd23_chianti_H_He.h5
    tardis tardis_example.yml output_spectrum.dat


Then plot the output_spectrum.dat with your favourite plotting program. Here's an example of how to do this with python
(the only thing you need to install is ipython and matplotlib --- in addition to TARDIS's requirements).

.. code-block:: python

    ipython --pylab
    >>> tardis_spec = loadtxt('output_spectrum.dat')
    >>> plot(tardis_spec[:,0], tardis_spec[:,1])

More atomic datasets can be downloaded from :ref:`atomic-data-download`.