Issues
------

.. _new-how-to-guide-report-issues:

How-To Guide: Report Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS is under constant development, so bugs, inconsistencies, and missing
functionality can occur. File issues on the official GitHub repository:

https://github.com/tardis-sn/tardis

Interested users are encouraged to contribute to TARDIS development by following
the development workflow.

.. _new-matterbridge:


Pull Requests
-------------

.. _new-how-to-guide-open-a-pull-request:

How-To Guide: Open a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you are ready for review or merge consideration:

1. Go to your fork, for example ``https://github.com/your-user-name/tardis``.
2. Use the branch dropdown to select the branch with your changes.
3. Click the pull request button.
4. Enter a title and an explanation of what you have done.
5. Include anything that needs particular attention, such as a complicated
   change or code you are unsure about.
6. If the request is not ready to merge, say so in the pull request message.
   Opening an unfinished pull request can still be a good way to start a
   preliminary code review.
7. Build documentation for the pull request using the documentation preview
   process so reviewers can see how notebooks, docstrings, and API docs render.

Example pull request description:

.. code-block:: markdown

   Summary
   - Update density parsing validation.
   - Add a regression or unit test for the new behavior.

   Checks
   - pytest tardis/io/model/readers/tests
   - ruff check tardis/io/model

   Notes for reviewers
   - Please check whether the error message is clear for invalid density types.


.. _new-how-to-guide-review-pull-requests:

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


.. _new-how-to-guide-update-a-pull-request:

How-To Guide: Update a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pull requests often remain open while feedback is added and other pull requests
are merged. If changes on trunk affect your work, rebase your feature branch on
top of the current trunk.

Fetch trunk and switch to your branch:

.. code-block:: shell

   git fetch upstream
   git checkout cool-feature


Create a temporary backup branch:

.. code-block:: shell

   git branch tmp cool-feature


Rebase onto upstream trunk:

.. code-block:: shell

   git rebase --onto upstream/master upstream/master cool-feature


If you are already on ``cool-feature``, the command can be shortened:

.. code-block:: shell

   git rebase upstream/master


When everything looks good, delete the backup branch:

.. code-block:: shell

   git branch -D tmp


If the feature branch is already on GitHub, force push after rebasing:

.. code-block:: shell

   git push -f origin cool-feature


Force pushing overwrites the branch on GitHub and can lose commits if used
incorrectly. Never force push to the main TARDIS repository, typically called
``upstream``, because that rewrites shared history.

