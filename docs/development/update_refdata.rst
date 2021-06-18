.. _update refdata:

*************************
Update the Reference Data
*************************

A special kind of tests are executed only when ``pytest`` is called alongside the ``--refdata`` flag.
These tests compares the output of the TARDIS code (mostly arrays) against the information stored
in the reference data files.

TARDIS stores reference data in the `tardis-refdata <https://github.com/tardis-sn/tardis-refdata>`_
repository. This repository also has a mirror hosted in Azure Pipelines (synchronized automatically by a 
GitHub workflow) since this Microsoft service does not have limitations in bandwith nor storage.

Sometimes, this data needs to be updated. The procedure to update these files manually is not trivial
and has been automatized recently thanks to the `NumFOCUS <https://numfocus.org/>`_ support.


=================
Default Procedure
=================

Imagine you are working on a new feature (or fix) for TARDIS, you have opened a pull request and the
reference data tests are failing in the testing pipeline. This could happen for many reasons:

A. There's a problem in your code.
B. Your code is OK, but the reference data is outdated.
C. The pipeline is broken.

If you think your could be dealing with scenario B, then:

#. Write ``/azp run compare-refdata`` in a comment on your PR.
#. Analyze the results and discuss if the reference data effectively requires an update.
#. Update the reference data by writing ``/azp run update-refdata`` on a new comment.

.. note::

    - If you don't have enough privileges to run the pipelines, tag a TARDIS developer capable of doing so.
    - If any of these two pipelines fail, please tag the :ref:`CI/CD responsible <team>`.

If everything went well, the reference data will have been updated by the TARDIS bot and the commit
message should include the pull request number that triggered the update.

================
Manual Procedure
================

The manual procedure is documented for debugging purposes and should not be used in general.

#. Activate the ``tardis`` environment.
#. Fork and clone the ``tardis-refdata`` repository.
#. Follow the instructions at the top of the notebook ``tardis-refdata/notebooks/ref_data_compare.ipynb``.
#. Go to your local ``tardis`` repository and make sure you are working on the branch you want to generate new reference data from.
#. Generate new reference data with ``pytest tardis --refdata=/path/to/tardis-refdata --generate-reference``.
#. Run the ``ref_data_compare.ipynb`` notebook and check the results.
#. Make a new branch in ``tardis-refdata``, push your new reference data and open a pull request.
