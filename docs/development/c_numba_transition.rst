******************
C-Numba transition
******************

Packet-by-packet comparison for simple tardis model (tardis_configv1_verysimple.yml) 
[PR929](https://github.com/tardis-sn/tardis/pull/929)

**Things to do:**
* line scattering (complete)
* electron scattering (complete)
* virtual packets (complete)
* macroatom (complete)
* packet logging (last_packet_interaction)

There are 4 code versions. There is upstream/master there is numba_montecarlo both of which have different seeding strategies and will never be numerically equal (those we will summarize as main). There are two versions that have additional code that will allow logging - those we will call c_compare and numba_compare (and will call the compare branch).

* virtual - only compare - sica (complete)
* virtual - only compare - ddc10
* integrated - main only - sica
* integrated - main only - ddc10
* estimators - compare only -sica (complete)
* estimators - compare only -ddc10 (complete)
* iterations - (compare/main) -sica (complete)
* iterations - (compare/main) -ddc10
* Need to fix negative boundary distance issue for vpackets: specific case occurs when a packet reaches a boundary with a close line so distance_trace_line = 0 and distance_boundary can be negative (shell ID does not advance). The models that shall be run is low velocity Si/Ca and DDC10.

**Post-merge**

* Start writing unit tests for the numba version (in progress)
* Improve speed.
* Close-lines need fixing in the numba code (C version not correct)
* Port formal integrator to numba 
