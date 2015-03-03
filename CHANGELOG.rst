1.1 (unreleased)
----------------

- Nothing changed yet.


1.0 (2015-03-03)
----------------

New Features
^^^^^^^^^^^^
- the main feature is the re-write of the montecarlo code to C, which comes
with a 4x speed increase.

Pull Requests
^^^^^^^^^^^^^

- [#210](https://github.com/tardis-sn/tardis/pull/210) a restructure of the C-file file structure (@wkerzendorf, @mklauser)
- [#7](https://github.com/tardis-sn/tardis/pull/7) include .c in montecarlo/src/ for the randomkit files (@mklauser)
- [#209](https://github.com/tardis-sn/tardis/pull/209) Testing/pandas version info (@wkerzendorf)
- [#205](https://github.com/tardis-sn/tardis/pull/205) As simpler `run_tardis` (@wkerzendorf)
- [#206](https://github.com/tardis-sn/tardis/pull/206) Minor but critical fix to binary search (reverse) (@ssim)
- [#199](https://github.com/tardis-sn/tardis/pull/199) Gnuify cmontecarlo.c and cmontecarlo.h files (@orbitfold)
- [#203](https://github.com/tardis-sn/tardis/pull/203) Docs/reporting bugs (@wkerzendorf)
- [#201](https://github.com/tardis-sn/tardis/pull/201) A critical warning was added to the model indicating that no packet has left the simulation through the outer boundary. (@mklauser)
- [#200](https://github.com/tardis-sn/tardis/pull/200) added a new location for the atomic data. (@wkerzendorf)
- [#196](https://ggithub.com/tardis-sn/tardis/pull/196) New binary search added (@mklauser)
- [#169](https://github.com/tardis-sn/tardis/pull/169) fix for #168 (@mklauser)
- [#194](https://github.com/tardis-sn/tardis/pull/194) Setup/fix requirements (@wkerzendorf)
- [#188](https://github.com/tardis-sn/tardis/pull/188) A macro to automagically fix inline issues (@orbitfold)
- [#190](https://github.com/tardis-sn/tardis/pull/190) General/fixing unit problems (@wkerzendorf)
- [#179](https://github.com/tardis-sn/tardis/pull/179) making sure that if last_no_of_packets is not specified that it is set t... (@wkerzendorf)
- [#178](https://github.com/tardis-sn/tardis/pull/178) Atomic/fix reprepare (@wkerzendorf)
- [#177](https://github.com/tardis-sn/tardis/pull/177) added from_yaml and from_config_dict to ConfigurationNameSpace (@wkerzendorf)
- [#175](https://github.com/tardis-sn/tardis/pull/175) Fix a few problems in cmontecarlo.c (@orbitfold, @wkerzendorf)
- [#174](https://github.com/tardis-sn/tardis/pull/174) fixes added for clang compile on mac (@wkerzendorf)
- [#171](https://github.com/tardis-sn/tardis/pull/171) Fix issues when using clang compiler (issue #170) (@orbitfold)
- [#167](https://github.com/tardis-sn/tardis/pull/167) fix for #165 (@mklauser, @wkerzendorf)
- [#3](https://github.com/tardis-sn/tardis/pull/3) testing the new t_inner fix (@wkerzendorf)
- [#164](https://github.com/tardis-sn/tardis/pull/164) Config/toggle validation (@wkerzendorf)
- [#151](https://github.com/tardis-sn/tardis/pull/151) WIP Montecarlo C Rewrite (@orbitfold)
- [#160](https://github.com/tardis-sn/tardis/pull/160) how to get constant density? broken? (@wkerzendorf)


0.9.2 (2014-06-12)
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
- fix for zeta factor outside of temperature + test included [#134, wkerzendorf]


New Features
^^^^^^^^^^^^
- added data_sources and version to `AtomData` set
- added atomic database table with downloads to documentation
- added more conversion routines for Species strings (e.g. Si IX) to util
- added new density models (power law and exponential) [mklauser]
- reimplementation of binary-search in C (towards faster/profileable code) [V. Jancauskas]
- added new astropy setup using astropy-helpers (thanks to @embray and @astrofrog for debugging help) [wkerzendorf #144]
- added a test that runs the full calculation and compares the spectrum output [wkerzendorf #144]
- added the new documentation validator [mklauser & wkerzendorf #134, #136]
- added a new configuration object for TARDIS [wkerzendorf #143]
- changed TARDISConfiguration -> Configuration [wkerzendorf #154]
- changed most small functions to C [V. Jancauskas #142]

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


