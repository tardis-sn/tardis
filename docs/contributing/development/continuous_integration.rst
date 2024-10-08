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
pipeline is as easy as making a pull request. To create a new workflow on GitHub, 
just create a new YAML file in ``.github/workflows``.

TARDIS Pipelines
----------------

Brief description of pipelines already implemented on Tardis

Documentation build pipeline
============================

A GitHub Actions workflow that builds and deploys the TARDIS documentation website.


Documentation preview pipeline
==============================

This workflow does not run on the main repository, just on forks.
See the :ref:`Documentation Preview <doc-preview>` section for more information.


Testing pipeline
================

The `testing pipeline`_ (CI) consists basically in the same job running twice
in parallel (one for each OS), plus extra steps to run the tests and upload the coverage results.


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