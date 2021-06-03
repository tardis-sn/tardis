*************************
Update the Reference Data
*************************

A special kind of tests are executed only when `pytest` is called alongside the `--refdata` flag.
These tests compares the output of the TARDIS code (mostly arrays) against the information stored
in the reference data files.

TARDIS stores reference data in the `tardis-sn/tardis-refdata <https://github.com/tardis-sn/tardis>`_
repository. This repository also has a mirror hosted in Azure Pipelines (synced automatically by a 
GitHub workflow) since the Microsoft service does not have limitations in bandwith nor storage.

Sometimes, this data needs to be updated. The workflow to update these files manually is not trivial
and has been automatized thanks to the `NumFOCUS <https://numfocus.org/>`_ support.


=================
Default Procedure
=================

Imagine you are working on a new feature (or fix) for TARDIS and you have opened a pull request. Then,
reference data tests are failing in the testing pipeline. This could happen for many reasons:

A. There's a problem in your code.
B. Your code is good, but the reference data is outdated.
C. The pipeline is broken.

If you suspect you are dealing with scenario B, then:

#. Write ``/azp run compare-refdata`` in a comment on your PR. If you don't have enough privileges, tag
    (@) a TARDIS developer.
#. Analyze the results.
#. If your suspicions are correct, update the reference data by writing a commenting
    ``/azp run update-refdata`` on a new comment (privileges are required too).
    
If everything went well, the reference data will have been updated by the TARDIS Bot, and the commit
message should include the TARDIS pull request number. If one of these two pipelines fail, tag the 
CI/CD responsible.

================
Manual Procedure
================

`--generate-refdata`.