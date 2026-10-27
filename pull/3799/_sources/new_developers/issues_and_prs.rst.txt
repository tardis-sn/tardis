Reporting Problems with TARDIS
------------------------------

.. _how-to-guide-report-issues:

How-To Guide: Report Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS is under constant development, so bugs, inconsistencies, and missing
functionality can occur. File issues on the official GitHub repository:

https://github.com/tardis-sn/tardis

Interested users are encouraged to contribute to TARDIS development by following
the development workflow.

.. _matterbridge:


Pull Requests
-------------

.. _how-to-guide-open-a-pull-request:

How-To Guide: Open a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you are ready for review or merge consideration:

1. Go to your fork, for example ``https://github.com/your-user-name/tardis``.
2. Use the branch dropdown to select the branch with your changes.
3. Click the pull request button.
4. Enter a title and an explanation of what you have done.
5. Include anything that needs particular attention, such as a complicated
   change or code you are unsure about.
6. If the request is not ready to merge, say so in the pull request message 
   and set the pull request to Draft using the drop-down menu.
   Opening an unfinished pull request can still be a good way to start a
   preliminary code review.
7. Build documentation for the pull request using the documentation preview
   process so reviewers can see how notebooks, docstrings, and API docs render.

Example pull request description:

.. code-block:: markdown
   ### :pencil: Description

   **Type:** :beetle: `bugfix` 

   The NLTE mask calculation for stimulated emission factor uses a lambda function acting on rows of a dataframe. This is extremely slow. Instead here is a vectorized version.

   Also includes various other performance improvements across the IIP plasma and the macroatom. Needs regression data.

   ### :pushpin: Resources

   [nlte_line_mask_benchmark.pdf](https://github.com/user-attachments/files/29320046/nlte_line_mask_benchmark.pdf)

   https://github.com/tardis-sn/tardis-regression-data/pull/103

   ### :vertical_traffic_light: Testing

   How did you test these changes?

   - [x] Testing pipeline
   - [ ] Other method (describe)
   - [ ] My changes can't be tested (explain why)


   ### :ballot_box_with_check: Checklist

   - [x] I requested two reviewers for this pull request
   - [ ] I updated the documentation according to my changes
   - [ ] I built the documentation by applying the `build_docs` label

   > **Note:** If you are not allowed to perform any of these actions, ping (@) a contributor.


.. _how-to-guide-review-pull-requests:

How-To Guide: Review Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GitHub supports line comments and discussion threads on pull requests. Sometimes
it is easier to demonstrate a suggested change by making a pull request on a
pull request.

Add pull request refs to the ``upstream`` remote in ``.git/config``:

.. code-block:: ini

   [remote "upstream"]
       url = git@github.com:tardis-sn/tardis.git
       fetch = +refs/heads/*:refs/remotes/upstream/*
       fetch = +refs/pull/*/head:refs/remotes/upstream/pr/*


Fetch upstream:

.. code-block:: shell

   git fetch upstream


You can then check out a pull request branch, for example pull request 116:

.. code-block:: shell

   git checkout upstream/pr/116


This leaves you in a detached ``HEAD`` state. Create a branch for your changes:

.. code-block:: shell

   git checkout -b helping-with-PR116


After making and committing your changes, push them to your fork:

.. code-block:: shell

   git push origin helping-with-PR116


.. _how-to-guide-update-a-pull-request:

How-To Guide: Update a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pull requests often remain open while feedback is added and other pull requests
are merged. If changes on trunk affect your work, merge the trunk into your feature branch.

Fetch trunk and switch to your branch:

.. code-block:: shell

   git fetch upstream
   git checkout cool-feature


Merge upstream trunk into your feature branch:

.. code-block:: shell

   git merge upstream/master

Resolve any merge conflicts and commit the changes.


If the feature branch is already on GitHub, push after merging:

.. code-block:: shell

   git push origin cool-feature



