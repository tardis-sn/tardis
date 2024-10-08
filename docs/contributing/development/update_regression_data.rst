.. _update regression-data:

*************************
Update the Regression Data
*************************

A special kind of tests are executed only when ``pytest`` is called alongside the ``--regression-data`` flag. These tests compare the output of the TARDIS code (mostly arrays) against the information stored in the regression data files.

TARDIS stores regression data in the `tardis-regression-data <https://github.com/tardis-sn/tardis-regression-data>`_ repository. Sometimes, this data needs to be updated. The procedure to update these files has been simplified, allowing for a more straightforward process.

=================
Default Procedure
=================

Imagine you are working on a new feature (or fix) for TARDIS, and you have opened a pull request. If the regression data tests are failing in the testing pipeline, this could happen for various reasons:

A. There's a problem in your code.  
B. Your code is OK, but the regression data is outdated.  
C. The pipeline is broken.

If you suspect scenario B, then:

#. Analyze the results to determine if the regression data requires an update.
#. Update your fork of the ``tardis-regression-data`` repository by pulling the latest changes and merging them into your local branch.
#. Push your updated fork to GitHub.

.. note::

    - If you do not have enough privileges to update the repository, tag a TARDIS developer capable of doing so.
    - If any issues arise during this process, please tag a `TARDIS team member <https://tardis-sn.github.io/team/community_roles/>`_ responsible for CI/CD.

If everything goes smoothly, your regression data will be updated in your fork, and you can proceed with your development.

================
Manual Procedure
================

The manual procedure is documented for debugging purposes and should not be used generally.

#. Activate the ``tardis`` environment.
#. Fork and clone the ``tardis-regression-data`` repository.
#. Follow any necessary instructions within your local copy.
#. Go to your local ``tardis`` repository and ensure you are working on the branch from which you want to generate new regression data.
#. Generate new regression data with ``pytest tardis --regression-data=/path/to/tardis-regression-data --generate-reference``.
#. Check your results and ensure everything is correct.
#. Make a new branch in ``tardis-regression-data``, push your new regression data, and open a pull request.

By following these updated procedures, you can efficiently manage and update regression data within your TARDIS project setup.
