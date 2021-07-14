.. _tardis-example:

****************************
Tardis Example Configuration
****************************

The simple and fast TARDIS setup is provided by ``tardis_example.yml`` which
may be obtained `here
<https://raw.githubusercontent.com/tardis-sn/tardis/master/docs/tardis_example.yml>`_. We suggest every new user of TARDIS to run this
setup first, which can be done using the :doc:`quickstart guide <../../quickstart/quickstart>`.

It calculates a spectrum for a Type Ia supernova model 13 days after explosion,
requesting a total output luminosity of

.. math::
    
    L = 10^{9.44}\, \mathrm{L}_{\odot}

A simple power-law density profile (seventh order polynomial fit to the Nomoto
et al. 1984 W7 model) is used together with a uniform composition, involving
only six elements. To avoid long run times only a moderate number of real and
virtual Monte Carlo packets are used. Also, very simple ionization and
excitation assumptions are adopted.

The following YAML file summarizes the tardis_example setup:

.. literalinclude:: tardis_example.yml
    :language: yaml

.. note::
    Due to the low number of packets, the simplistic ionization and excitation
    treatments and the reduced abundance set, this TARDIS setup serves for
    illustrative purposes and not for detailed SNe Ia spectral synthesis
    calculations.
