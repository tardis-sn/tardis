**************************
Debugging numba_montecarlo
**************************
To facilitate more in-depth debugging when interfacing with the `montecarlo_numba`
module, we provide a set of debugging configurations. PyCharm debugging
configurations, in addition to related scripts and .yml files, are contained in
`tardis.scripts.debug`. Currently, these include the ability to run TARDIS
in asingle-packet mode, with the packet seed identified at debug time.
`tardis_example_single.yml` is the configuration filethat is used to set up the
single-packet TARDIS run; `run_numba_single.py` is thePython script that runs
this .yml file; `run_numba_single.xml` is the PyCharmdebug configuration file
that can be used in conjunction with the above files.

Note that this method is EXPERIMENTAL.
