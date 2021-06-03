tardis_example
--------------

The simple and fast TARDIS setup is provided by ``tardis_example.yml`` which
may be obtained from `tardis-setups
<https://github.com/tardis-sn/tardis-setups>`_ repository. It is located in
the ``tardis-setups/2014/2014_kerzendorf_sim/appendix_A1`` subfolder. We suggest every new user of TARDIS to run this
setup first.

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
    
See the following link for an example of running TARDIS with this setup in a Jupyter notebook:

.. toctree::
    :maxdepth: 1
    
    run_simple_example

