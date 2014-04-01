.. _development-workflow:

=======================
Workflow for Developers
=======================

In this document, we refer to the Tardis ``master`` branch as the *trunk*.

.. _forking:

Creating a fork
===============

You need to do this only once for each package you want to contribute to. The
instructions here are very similar to the instructions at
http://help.github.com/fork-a-repo/ |emdash| please see that page for more
details. We're repeating some of it here just to give the specifics for the
TARDIS_ project, and to suggest some default names.

Set up and configure a GitHub account
-------------------------------------

If you don't have a GitHub account, go to the GitHub_ page and make one.

You then need to configure your account to allow write access |emdash| see
the `Generating SSH keys
<http://help.github.com/articles/generating-ssh-keys>`_ help on `GitHub Help`_.

Create your own fork of a repository
------------------------------------

The following example shows how to fork the core ``astropy`` repository, but
the same applies to other packages:

#. Log into your GitHub_ account.

#. Go to the `TARDIS GitHub`_ home page.

#. Click on the *fork* button:

   .. image:: forking_button.png

   After a short pause and an animation of Octocat scanning a book on a flatbed
   scanner, you should find yourself at the home page for your own forked copy
   of TARDIS_.

Setting up the fork to work on
------------------------------

.. _linking-to-upstream:

Overview
^^^^^^^^

This is done using::

    git clone git@github.com:your-user-name/tardis.git
    cd astropy
    git remote add upstream git://github.com/tardis-sn/tardis.git

In detail
^^^^^^^^^

#. Clone your fork to the local computer::

    git clone git@github.com:your-user-name/tardis.git

#. Change directory to your new repo::

    cd astropy

   Then type::

    git branch -a

   to show you all branches.  You'll get something like::

    * master
    remotes/origin/master

   This tells you that you are currently on the ``master`` branch, and
   that you also have a ``remote`` connection to ``origin/master``.
   What remote repository is ``remote/origin``? Try ``git remote -v`` to
   see the URLs for the remote connections.  They will point to your GitHub
   fork.

   Now connect to the TARDIS repository, so you can merge in changes from the
   trunk::

    cd
    git remote add upstream git://github.com/tardis-sn/tardis.git

   ``upstream`` is just the arbitrary name we're using to refer to the main
   TARDIS_ repository.

   Note that we've used ``git://`` for the URL rather than ``git@``. The
   ``git://`` URL is read only. This means that we can't accidentally (or
   deliberately) write to the upstream repo, and we are only going to use it
   to merge into our own code.

   Just for your own satisfaction, show yourself that you now have a new
   remote connection with ``git remote -v show``, which should give you
   something like::

    upstream   git://github.com/tardis-sn/tardis.git (fetch)
    upstream   git://github.com/tardis-sn/tardis.git (push)
    origin     git@github.com:your-user-name/tardis.git (fetch)
    origin     git@github.com:your-user-name/tardis.git (push)

   Your fork is now set up correctly, and you are ready to hack away.

Installing TARDIS in develop mode
==================================

TARDIS is designed so that it can generally be used directly out of the source
tree by using ``import `` when running Python in the source of an
TARDIS repository clone.

#. Install TARDIS_ in develop mode::

       $ python setup.py develop

   This semi-permanently installs TARDIS on your path in such a way that
   ``tardis`` is always imported from your repository clone regardless of your
   working directory.  This way any edits you make to the code in your
   repository will always be immediately available next time you start a Python
   interpreter and ``import tardis``.

Workflow summary
================

This section gives a summary of the workflow to follow once you have
successfully forked the repository. The details for each of these steps are
given in the following sections.

* Don't use your ``master`` branch for anything.  Consider deleting it.

* When you are starting a new set of changes, fetch any changes from the
  trunk, then start a new *feature branch* from that.

* Make a new branch for each separable set of changes |emdash| "one task, one
  branch" (`ipython git workflow`_).

* Name your branch for the purpose of the changes, for example
  ``bugfix-for-issue-14`` or ``refactor-database-code``.

* If you can possibly avoid it, don't merge the trunk or any other branches into
  your feature branch while you are working.

* If you do find yourself merging from the trunk, consider
  :ref:`rebase-on-trunk`

* Ask on the `tardis-sn-dev mailing list`_ if you get stuck.

* Once your code is nearing completion, run the test suite to ensure
  you have not accidentally caused regressions, and add new tests to ensure
  your contribution behaves correctly (see :ref:`testing-guidelines`).

* Issue a pull request on github!

* As the code is converging to a final state, ensure your
  documentation follows the guidelines (see :ref:`documentation-guidelines`).

* Once your code is ready to be accepted, please add an entry to the changelog
  (see :ref:`changelog-format`).  If you're sure where to put this, please at
  least suggest a brief (one or two sentence) description of your change so
  that another Astropy developer can add it to the changelog.

This way of working helps to keep work well organized, with readable history.
This in turn makes it easier for project maintainers (that might be you) to
see what you've done, and why you did it.

See `linux git workflow`_ and `ipython git workflow`_ for some explanation.

Deleting your master branch
===========================

It may sound strange, but deleting your own ``master`` branch can help reduce
confusion about which branch you are on.

.. _update-mirror-trunk:

Updating the mirror of trunk
============================

From time to time you should fetch the upstream (trunk) changes from GitHub::

   git fetch upstream

This will pull down any commits you don't have, and set the remote branches to
point to the right commit. For example, 'trunk' is the branch referred to by
(remote/branchname) ``upstream/master``, and if there have been commits since
you last checked, ``upstream/master`` will change after you do the fetch.

.. _make-feature-branch:

Making a new feature branch
===========================

When you are ready to make some changes to the code, you should start a new
branch. Branches that are for a collection of related edits are often called
'feature branches'.

Making a new branch for each set of related changes will make it easier for
someone reviewing your branch to see what you are doing.

Choose an informative name for the branch to remind yourself and the rest of
us what the changes in the branch are for. For example ``add-ability-to-fly``,
or ``buxfix-for-issue-42``.

::

    # Update the mirror of trunk
    git fetch upstream

    # Make new feature branch starting at current trunk
    git checkout upstream/master # checking out the newest
    git checkout -b my-new-feature

Generally, you will want to keep your feature branches on your public GitHub_
fork. To do this, you `git push`_ this new branch up to your
github repo. Generally (if you followed the instructions in these pages, and
by default), git will have a link to your GitHub repo, called ``origin``. You
push up to your own repo on GitHub with::

   git push origin my-new-feature

In git >= 1.7 you can ensure that the link is correctly set by using the
``--set-upstream`` option::

   git push --set-upstream origin my-new-feature

From now on git will know that ``my-new-feature`` is related to the
``my-new-feature`` branch in the GitHub repo.

.. _edit-flow:

The editing workflow
====================

Overview
--------

Make changes, test, and::

   git add my_new_file
   git commit -m 'NF - some message'
   git push

In more detail
--------------

#. Make some changes

#. Once you are a bit further along, test your changes do not lead to
   regressions, and add new tests (see :ref:`testing-guidelines`).::

     python setup.py test

   If you have sphinx installed, you can also check that the documentation
   builds and looks correct::

     python setup.py build_sphinx

   The last line should just state ``build succeeded``, and should not mention
   any warnings.  (For more details, see :ref:`documentation-guidelines`.)

#. See which files have changed with ``git status`` (see `git status`_).
   You'll see a listing like this one::

     # On branch ny-new-feature
     # Changed but not updated:
     #   (use "git add <file>..." to update what will be committed)
     #   (use "git checkout -- <file>..." to discard changes in working directory)
     #
     #    modified:   README
     #
     # Untracked files:
     #   (use "git add <file>..." to include in what will be committed)
     #
     #    INSTALL
     no changes added to commit (use "git add" and/or "git commit -a")

#. Check what the actual changes are with ``git diff`` (see `git diff`_).

#. Add any new files to version control with ``git add new_file_name`` (see
   `git add`_).

#. Add any modified files that you want to commit using
   ``git add modified_file_name``  (see `git add`_).

#. Once you are ready to commit, check with ``git status`` which files are
   about to be committed::

    # Changes to be committed:
    #   (use "git reset HEAD <file>..." to unstage)
    #
    #    modified:   README

   Then use ``git commit -m 'A commit message'``. The ``m`` flag just
   signals that you're going to type a message on the command line. The `git
   commit`_ manual page might also be useful.

#. Push the changes up to your forked repo on GitHub with ``git push`` (see
   `git push`_).
