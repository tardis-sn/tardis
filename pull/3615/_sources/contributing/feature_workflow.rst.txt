Feature Workflow
------------------

.. _new-tutorial-add-a-new-feature:

Tutorial: Add a New Feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS development aims to follow test-driven development. A new feature should
move from a documented need, to a failing test, to implementation, to passing
tests, to documentation and changelog updates.

1. Document the feature to be added in an issue.
2. If you are unsure whether the feature already exists, ask the mailing list.
3. Write a test that demonstrates the feature that will be added.
4. Run the test to verify that it fails in the way you expect.
5. If the test fails unexpectedly, the test may be wrong. Ask the group for
   guidance.
6. If the test passes, the feature may already exist. In that case, you have
   added test coverage for existing behavior.
7. Add the feature.
8. Run the test to verify that it passes.
9. Write documentation about the feature.
10. Close the issue or partial pull request and add the change to the changelog.

Large features should usually be broken into small, quantifiable goals that can
be acted on separately.

For example, if a change adds support for a new density configuration, a small
development path would be: open an issue describing the new density option, add
a failing test near the existing model or configuration tests, update
``tardis/io/model/parse_density_configuration.py``, run the targeted tests, and
then document the new configuration behavior.

.. _new-how-to-guide-introduce-a-new-feature:

How-To Guide: Introduce a New Feature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before implementing a new feature, make the feature reviewable:

1. Open or find an issue that states the science requirement, the expected
   behavior, and any scientific or data assumptions.
2. Identify the smallest code path that owns the behavior. For example, a new
   model input format likely belongs near ``tardis/io/model/``, while a new
   visualization behavior likely belongs near ``tardis/visualization/``.
3. Add a failing unit test or regression test that captures the new behavior.
4. Implement the smallest code change that makes the test pass.
5. Add or update documentation for the feature.
6. Run targeted tests, Ruff, and any relevant documentation build.
7. Open a pull request with the issue link, test commands, and reviewer notes.

Example feature checklist:

.. code-block:: markdown

   Issue
   - Adds support for a new density option in CSVY model input.

   Code
   - Update tardis/io/model/parse_density_configuration.py.
   - Add or update tests under tardis/io/model/readers/tests.

   Validation
   - pytest tardis/io/model/readers/tests
   - ruff check tardis/io/model

   Documentation
   - Update the relevant configuration or model reader page.


.. _new-explanation-developer-workflow:

Explanation: Developer Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The developer workflow is intended for anyone contributing to TARDIS
development. Much of the workflow is taken from Astropy's development workflow,
with credit to the Astropy team.

The workflow emphasizes:

- Filing issues for bugs, inconsistencies, missing functionality, and planned
  features.
- Working from a fork.
- Keeping an ``upstream`` remote connected to the main TARDIS repository.
- Avoiding development on local ``master``.
- Starting each separable task from current upstream trunk.
- Keeping work organized through feature branches.
- Running tests before review.
- Updating documentation.
- Opening pull requests for review and merge.

This helps maintain readable history and helps maintainers understand what
changed and why.

.. _new-explanation-test-driven-development:

Explanation: Test-Driven Development
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS aims to use test-driven development. The preferred path starts by
describing a feature, then writing a test for the missing behavior, then
implementing the feature until the test passes.

Unexpected test failures can indicate that the test is wrong. A passing test
before implementation can mean the feature already exists and only lacked
coverage.

Learn More About How To Work As A Test-Driven Developer Here
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use these resources to connect the TARDIS workflow to common test-driven
development practices:

- Pytest getting started:
  https://docs.pytest.org/en/stable/getting-started.html
- Pytest fixtures:
  https://docs.pytest.org/en/stable/explanation/fixtures.html
- Pytest parametrized tests:
  https://docs.pytest.org/en/stable/how-to/parametrize.html
- Python ``unittest`` documentation:
  https://docs.python.org/3/library/unittest.html
- Astropy testing guidelines:
  https://docs.astropy.org/en/stable/development/testguide.html
