Commit Message Guidelines
=========================

When to commit? Ask yourself the following questions:

* Do the changes you have implemented represent a single cohesive feature? Break a large task down in to several subtasks.
* Is the system in a build-able state after the commit? Don't do commits that you know will leave the system in a broken state.
* Can your commit be safely rolled-back? Will the system be left in a broken state if you roll-back one commit? If so then maybe several commits should be merged.

Preferably a commit should affect one file and do one thing only. For example resist the temptation to fix indentation in an unrelated file or in an unrelated class/method.

When committing changes do not use -m option. This is only suitable for really short commit messages. Instead omit the message flag and you will be taken to your default text editor where you can enter the message. If you have already pushed the changes to a remote repository or if you have committed a change locally with a message that does not conform to the guidelines given below please use the following instructions to change commit messages:

https://help.github.com/articles/changing-a-commit-message/

The following are the requirements for the commit messages:

* The commit message consists of the following sections:

  * One-line description
  * Empty line
  * Body of the commit message

* One-line description consists of the following:

  * The commit message has a one-line description. Keep it short, preferablypreferrably under 80 characters. If you're squashing a bug, the bug number you're squashing should be included in this one-line description, like so at the end: (bug #12345). Omit final period (full stop).
  * Sometimes the one-line description is sufficient, for small changes. In that case, we use the following prefixes to give a description of what the small change is:

    * maint: for reorganisation of the sources that do not change the source. Regular merge commits are a prominent example.
    * doc: for changes to the documentation.

  * If your change only touches one file, then the filename of that file should be the prefix of the one-line description. The function or class that is being modified should be included in the parenthetical remark, as in the full body of the commit message.

* If there is more than one file touched in different ways and the one-line description isn't enough to describe all changes, the commit message needs a full-body description. The full-body description would consist of a description of the changes and a list of changes to individual files. The description of changes should answer the following three questions:

  * Why is this change necessary?
  * How does it address the issue?
  * What side effects does this change have?

If the changes are simple and/or the answers to the above questions are obvious from the context then you can ommit it.

  Individual files you touched get their changes separately enumerated. If there a particular C or C++ function you changed, the general format is:

  * file.py (class1.function1): Make change foo.
     (class2.function2): Make change bar.
  * Note the newline for describing the change to a different function. General GNU style guidelines can be followed here, so that similar changes to different files can be grouped.
  * For style points, please also write "New function" instead of "Added function" or "Return retval" instead of "Changed to return retval".
  * And please, never "Fixed bug" or similar. That doesn't add any specific information about what was changed.

* The changelog info in the commit message should say what changed, not why. If you need to explain why, that info should probably go in comments in the code.

Commit Message Examples
-----------------------

Standard commit messages:

Example 1:

::

   ensure PEP8 is followed in our code (issue #427)

   * datahandler.py, docs/conf.py, interface.py, widgets.py: changes to make sure that PEP8 is followed

Example 2:

::

   allow user to specify the number of threads in the configuration file

   * data/tardis_config_definition.yml: added nthreads attribute under montecarlo
   * model.py (Radial1DModel.simulate): pass the required number of threads to the C part of montecarlo routines

Example 3:

::

   added --with-vpacket-logging build option (bug #428)

   * setup.py: added --with-vpacket-logging flag for install, build and develop commands
   * montecarlo/setup_package.py: define macro WITH_VPACKET_LOGGING if --with-vpacket-logging flag is used


For short and very simple commit messages under 80 characters:

Example 1:

::

   tests/test_atomic.py: added a test for atomic.read_levels_data

Example 2:

::

   montecarlo/src/cmontecarlo.c: fix indentation

Maintenance or docs commits, no actual code changed:

Example 1:

::

   maint: updated .travis.yml file to do a rehash before executing conda
