******************
Developer Workflow
******************

If you want to contribute to the development of TARDIS, please read through
this page and follow the guidelines outlined here.  Many of the Development
workflow is taken from `Astropy <http://docs.astropy.org
/en/stable/development/workflow/development_workflow.html>`_ and credit belongs
to the Astropy team for designing it.

Preparation and Working with Git
================================

The first step is to setup up a python environment. :doc:`python_environment`
explains how to do that. Some of this information is also contained in the
:doc:`Installation guide <../installation>`. Once the python environment is
running consult :doc:`git_workflow` to learn about interacting with git and
github for the purpose of developing TARDIS.


General Workflow to add a new feature
=====================================

In TARDIS we aim to stick to a test driven development. This uses the testing
framework extensively starting with a test that shows this feature lacking via
the implementation of the feature until the merging of the code to the main
repository.

In most cases we try to break down big features into small, quantifiable goals
which are then acted upon.

* Document feature to be added in an issue and maybe ask the mailing
  list if this feature exists
* Write a test that demonstrates what feature will be added.
* Run the test to verify that it fails in the way you think it should.
* If it fails in an unexpected way, your test may be wrong. This is a
  great time to ask the group for guidance
* If it passes, you are done! You just added test coverage to an
  already existing feature, and that is great! (unlikely)
* Add the feature (also known as "a simple matter of programming").
* Run the test to verify that it passes.
* Write documentation about your feature.
* close issue/partial PR and add to changelog.


.. include:: release_workflow.inc
