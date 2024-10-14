.. include:: links.inc

**********************
Continuous Integration
**********************

We use a so-called `continuous integration`_ workflow with TARDIS.
This means that each time a change is proposed (via pull request)
or a change is merged into the *master* branch, a service will
clone the repository, checkout to the current commit and execute
all the TARDIS tests. This helps us to detect bugs immediately.

GitHub Actions
--------------

A pipeline (or a workflow) is essentially a :term:`YAML` configuration file 
with different sections such as variables, jobs and steps. These files
run commands or tasks when they are triggered by some event, like a 
commit being pushed to a certain branch.

Currently, we use GitHub Actions to run all of our pipelines. Making changes to an existing 
pipeline is as easy as making a pull request. To create a new GitHub Action workflow, 
just create a new YAML file in ``.github/workflows``.

TARDIS Pipelines
----------------

Brief description of pipelines already implemented on TARDIS

Documentation build pipeline
============================

A GitHub Actions workflow that builds and deploys the TARDIS documentation website.


Documentation preview pipeline
==============================

This workflow does not run on the main repository, just on forks.
See the :ref:`Documentation Preview <doc-preview>` section for more information.


Testing pipeline
================

The `testing pipeline`_ (CI) comprises of six concurrent jobs. Each of these jobs runs three types of tests across two distinct categories—continuum and non-continuum—and supports three different operating system platforms. Additionally, there are extra steps involved in executing the tests and uploading the coverage results


Authors pipeline
================

This pipeline runs a notebook located in ``tardis-zenodo`` repository and
pushes a new version of ``.zenodo.json`` to the root of ``tardis``
repository if new committers are found (or author order changes). The
rendered notebook is uploaded to the pipeline results as an artifact.

.. warning :: Fails if some author name is incomplete (due to an incomplete
          GitHub profile) or duplicated (committed with more than one 
          email address). In both cases update ``.mailmap`` to fix it.

In the near future we want to auto-update the citation guidelines in the
``README.rst`` and the documentation.


Release pipeline
================

Publishes a new release of TARDIS every sunday at 00:00 UTC.


TARDIS Carsus Compatibility Check
=================================
The TARDIS Carsus Compatibility Check or the "Bridge" compares reference data 
generated with different versions of Carsus. It consists of two jobs- a "carsus-build" job to 
generate an atomic file with the latest version of Carsus and a "tardis-build" job 
to generate a new reference data with it. These two reference data files are compared using the 
`this notebook <https://github.com/tardis-sn/tardis-refdata/blob/master/notebooks/ref_data_compare_from_paths.ipynb>`_.
The workflow has a ``workflow_dispatch`` event so that it can be triggered manually, but is also 
triggered every week due to the "save-atomic-files" workflow. 


The Save Atomic Files Workflow
==============================
The Save Atomic Files workflow runs every week but can also be triggered manually.
It runs the "Bridge" and sends an artifact containing the generated atomic data file
and the comparison notebook to Moria. This workflow has a separate job to indicate if the 
bridge has failed.