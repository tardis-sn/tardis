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
================

Brief description of pipelines already implemented on TARDIS

Cache Keys in TARDIS CI
-----------------------

TARDIS uses specific cache key formats to efficiently store and retrieve data during CI runs:

1. **Regression Data Cache Keys**
   - Format: ``tardis-regression-<data-type>-<hash>-v1``
   - Examples:
     - ``tardis-regression-atom-data-sparse-<hash>-v1`` - For atomic data cache
     - ``tardis-regression-full-data-<hash>-v1`` - For full TARDIS regression data cache
   - Used in: ``setup_lfs`` action

2. **Environment Cache Keys**
   - Format: ``tardis-conda-env-<os-label>-<hash>-v1``
   - Examples:
     - ``tardis-conda-env-linux-<hash>-v1`` - For Linux conda environment
     - ``tardis-conda-env-macos-<hash>-v1`` - For macOS conda environment
   - Used in: ``setup_env`` action

.. warning::
   - The version suffix (-v1) allows for future cache invalidation if needed.
   - The `lfs-cache` workflow will fail if the cache is not available and will not pull LFS data by default. 
   - However, if the `allow_lfs_pull` label is added to the PR, the workflow will pull LFS data. Please note that this is to be used sparingly and only with caution.

Streamlined Steps for TARDIS Pipelines
========================================

We have a common set of steps which are utilized in TARDIS pipelines to streamline the process:

Common Steps
------------

1. **Use `setup_lfs` Action and `lfs-cache` workflow**
   - If you need access to regression or atomic data, incorporate the `setup_lfs` action to ensure proper handling of large file storage.
   - The `lfs-cache` workflow is used to cache the regression data and atomic data and to check if the cache is available.

2. **Use `setup_env` Action**
   - To configure your environment effectively, utilize the `setup_env` action. This will help establish the necessary variables and settings for your pipeline.

3. **Run Configuration**
   - Ensure that your pipeline runs with the appropriate shell settings. You can define this in your YAML configuration as follows:

    .. code-block:: yaml
        
        defaults:
          run:
            shell: bash -l {0}


Documentation build pipeline
============================

A GitHub Actions workflow that builds and deploys the TARDIS documentation website.


Documentation preview pipeline
==============================

This workflow does not run on the main repository, just on forks.
See the :ref:`Documentation Preview <doc-preview>` section for more information.


Testing pipeline
================

The `testing pipeline`_ (CI) comprises multiple concurrent jobs. Each of these jobs runs tests across two distinct categories—continuum and rpacket tracking—and supports two different operating systems. Additionally, there are extra steps involved in executing the tests and uploading the coverage results


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

Publishes a new release of TARDIS every Sunday at 00:00 UTC.


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

The Regression Data Comparison Workflow
======================================

The Regression Data Comparison workflow compares the regression data between the current branch and the base branch on pull requests. It only runs on pull requests and not on the master branch. The workflow generates regression data for the latest commit on the pull request and compares it with the master branch using the comparison notebook. The notebook is then uploaded as an artifact and pushed to reg-data-comp repository for previews in the bot comment.

.. note :: The workflow exports images from the comparison notebook and embeds them in the bot comment. Unless there are any key changes to any of the HDF files in the regression data the bot will only show two images, one containing the spectrum change and another containing relative changes in the keys. If there are any key changes, the bot will show three images, the additional one visualizing the key changes.

The `LFS-cache` workflow
=======================

The `LFS-cache` workflow caches the regression data and atomic data and can be triggered either manually or when there is a push to the main branch of the regression data repository. This is mainly responsible for doing LFS pulls when necessary and caching objects while the `setup-lfs` action is used to restore the cached objects. Both fail if the cache is not available.
