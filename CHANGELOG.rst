0.9.2 (unreleased)
------------------

Bugfixes
^^^^^^^^

- "print" replaced with logger for classical nebular
- logger statement added for coronal approximation
- warning added to documentation since plasma is out of date (temp
  solution only) #108
- fix to binary search to deal with packets at end of line list
- warning added for density file readin outside tabulated range
- fix to NLTE solver (treats stimulated emission directly); tests included


New Features
^^^^^^^^^^^^
- added data_sources and version to `AtomData` set
- added atomic database table with downloads to documentation
- added more conversion routines for Species strings (e.g. Si IX) to util
- added new density models (power law and exponential) [mklauser]
- reimplementation of binary-search in C (towards faster/profileable code) [V. Jancauskas]


0.9.1 (2014-02-03)
------------------

New Features
^^^^^^^^^^^^

- bugfix release of TARDIS
- updated example section in the documentation
- now with working test coverage setup (https://coveralls.io/r/tardis-sn/tardis)


Bugfixes
^^^^^^^^

- missing ez_setup.py and setuptools_bootstrap.py included now
- reading of ascii density and abundance profiles
- several fixes in the documentation


