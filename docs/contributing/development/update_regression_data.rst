.. _update regression-data:

*************************
Update the Regression Data
*************************

A special kind of tests are executed only when ``pytest`` is called alongside the ``--tardis-regression-data`` flag. These tests compare the output of the TARDIS code (mostly arrays) against the information stored in the regression data files.

TARDIS stores regression data in the `tardis-regression-data <https://github.com/tardis-sn/tardis-regression-data>`_ repository. Sometimes, this data needs to be updated. The procedure to update these files has been simplified, allowing for a more straightforward process.

Imagine you are working on a new feature (or fix) for TARDIS, and you have opened a pull request. If the regression data tests are failing, this could happen for various reasons:

A. There's a problem in your code.
B. Your code is OK, but the regression data is outdated.
C. The pipeline is broken.

If you suspect scenario B, please follow these instructions:

#. Activate the ``tardis`` environment.
#. Fork and clone the ``tardis-regression-data`` repository.
#. Follow any necessary instructions within your local copy.
#. Go to your local ``tardis`` repository and ensure you are working on the branch from which you want to generate new regression data.
#. Generate new regression data with ``pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference``.
#. Check your results and ensure everything is correct.
#. Make a new branch in ``tardis-regression-data``, push your new regression data, and open a pull request.

If any issues arise during this process, please tag a `TARDIS team member <https://tardis-sn.github.io/people/collaboration/>`_ responsible for CI/CD.