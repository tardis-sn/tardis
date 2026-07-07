Git Workflow
------------

.. _explanation-git-branching-strategy:

Explanation: Git Branching Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the developer workflow, ``master`` is referred to as the trunk. Feature work
should begin from ``upstream/master``.

Feature branches make reviews easier because each branch contains a related set
of edits. Contributors should avoid merging trunk or other branches into a
feature branch when possible. If trunk changes need to be incorporated, rebasing
is preferred because it replays feature commits on top of the latest trunk and
keeps history easier to read.

A rebase transforms a history like:

.. code-block:: text

         A---B---C cool-feature
        /
   D---E---F---G trunk


into:

.. code-block:: text

                 A'--B'--C' cool-feature
                /
   D---E---F---G trunk


Merge conflicts can occur if the feature branch and trunk changed the same
files. Resolve conflicts using the Git rebase documentation and related merge
resolution guidance.

.. _reference-git-command-reference:

Reference: Git Command Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clone a fork:

.. code-block:: shell

   git clone git@github.com:your-user-name/tardis.git


Inspect branches:

.. code-block:: shell

   git branch -a


Inspect remotes:

.. code-block:: shell

   git remote -v
   git remote -v show


Add upstream:

.. code-block:: shell

   git remote add upstream https://github.com/tardis-sn/tardis.git


Fetch upstream:

.. code-block:: shell

   git fetch upstream


Start a branch from trunk:

.. code-block:: shell

   git checkout upstream/master
   git checkout -b my-new-feature


Push a branch:

.. code-block:: shell

   git push origin my-new-feature
   git push --set-upstream origin my-new-feature


Check status and diff:

.. code-block:: shell

   git status
   git diff


Stage files:

.. code-block:: shell

   git add new_file_name
   git add modified_file_name


Commit:

.. code-block:: shell

   git commit -m "A commit message"


Rebase:

.. code-block:: shell

   git fetch upstream
   git checkout cool-feature
   git branch tmp cool-feature
   git rebase upstream/master


Abort a rebase:

.. code-block:: shell

   git rebase --abort


Force push a rebased branch to your fork:

.. code-block:: shell

   git push -f origin cool-feature


Reset to a backup branch:

.. code-block:: shell

   git reset --hard tmp


Inspect reflog:

.. code-block:: shell

   git reflog show cool-feature

.. _how-to-guide-recover-from-git-mistakes:

How-To Guide: Recover from Git Mistakes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a rebase goes wrong while it is in progress:

.. code-block:: shell

   git rebase --abort


If you notice the problem after the rebase and made a backup branch:

.. code-block:: shell

   git reset --hard tmp


If you forgot to make a backup branch, inspect the reflog:

.. code-block:: shell

   git reflog show cool-feature


Example reflog:

.. code-block:: text

   8630830 cool-feature@{0}: commit: BUG: io: close file handles immediately
   278dd2a cool-feature@{1}: rebase finished: refs/heads/my-feature-branch onto 11ee694744f2552d
   26aa21a cool-feature@{2}: commit: BUG: lib: make seek_gzip_factory not leak gzip obj


Reset to the point before the bad rebase:

.. code-block:: shell

   git reset --hard cool-feature@{2}

.. _how-to-guide-start-a-feature-branch:

How-To Guide: Start a Feature Branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the developer workflow, the TARDIS ``master`` branch is the trunk.

1. Do not use your local ``master`` branch for development. Consider deleting it
   to reduce confusion.
2. Fetch the latest upstream changes:

   .. code-block:: shell

      git fetch upstream


3. Start from the current upstream trunk:

   .. code-block:: shell

      git checkout upstream/master


4. Create a new feature branch:

   .. code-block:: shell

      git checkout -b my-new-feature


Use a new branch for each separable set of changes: one task, one branch. Choose
an informative name, such as ``bugfix-for-issue-14``,
``refactor-density-parser``, or ``update-regression-data-docs``.

Push the branch to your fork:

.. code-block:: shell

   git push origin my-new-feature


With Git 1.7 or newer, you can set the upstream branch:

.. code-block:: shell

   git push --set-upstream origin my-new-feature


.. _how-to-guide-use-the-editing-workflow:

How-To Guide: Use the Editing Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Make your changes.
2. Run tests to check for regressions:

   .. code-block:: shell

      pytest tardis


3. If Sphinx is installed, check the documentation build:

   .. code-block:: shell

      cd docs
      make html


   The build should succeed and should not report warnings.

4. Check changed files:

   .. code-block:: shell

      git status


5. Inspect the actual changes:

   .. code-block:: shell

      git diff


6. Add new files:

   .. code-block:: shell

      git add new_file_name


7. Add modified files you want to commit:

   .. code-block:: shell

      git add modified_file_name


8. Check what will be committed:

   .. code-block:: shell

      git status


9. Commit:

   .. code-block:: shell

      git commit -m "A commit message"


10. Push to your fork:

    .. code-block:: shell

       git push


For a small parser change, a typical local loop might be:

.. code-block:: shell

   ruff check tardis/io/model/parse_density_configuration.py
   pytest tardis/io/model/readers/tests
   git diff
   git add tardis/io/model/parse_density_configuration.py
   git commit -m "Improve density parser validation"
   git push
