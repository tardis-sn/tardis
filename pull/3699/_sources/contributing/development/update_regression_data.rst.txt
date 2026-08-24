.. _update regression-data:

**************************
Update the Regression Data
**************************

A special kind of tests are executed only when ``pytest`` is called alongside the ``--tardis-regression-data`` flag. These tests compare the output of the TARDIS code (mostly arrays) against the information stored in the regression data files.

TARDIS stores regression data in the `tardis-regression-data <https://github.com/tardis-sn/tardis-regression-data>`_ repository. Sometimes, this data needs to be updated. The procedure to update these files has been simplified, allowing for a more straightforward process.

Imagine you are working on a new feature (or fix) for TARDIS, and you have opened a pull request. If the regression data tests are failing, this could happen for various reasons:

A. There's a problem in your code.
B. Your code is OK, but the regression data is outdated.
C. The pipeline is broken.

If you suspect scenario B, please follow these instructions:

#. Activate the ``tardis`` environment.
#. Fork and clone the ``tardis-regression-data`` repository.
#. Follow any necessary instructions within your local copy.
#. Go to your local ``tardis`` repository and ensure you are working on the branch from which you want to generate new regression data.
#. Generate new regression data with ``pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference``.
#. Check your results and ensure everything is correct.
#. Make a new branch in ``tardis-regression-data``, push your new regression data, and open a pull request.

If any issues arise during this process, please tag a `TARDIS team member <https://tardis-sn.github.io/people/collaboration/>`_ responsible for CI/CD.


FAQ
===

In case you encounter an error similar to this:

.. code-block:: shell

    $ git push origin <branch-name>
    Uploading LFS objects: 100% (1/1), 13 MB | 0 B/s, done.                                                                                      
    Enumerating objects: 17, done.
    Counting objects: 100% (17/17), done.
    Delta compression using up to 192 threads
    Compressing objects: 100% (7/7), done.
    Writing objects: 100% (9/9), 732 bytes | 366.00 KiB/s, done.
    Total 9 (delta 4), reused 0 (delta 0), pack-reused 0
    remote: Resolving deltas: 100% (4/4), completed with 4 local objects.
    remote: error: GH008: Your push referenced at least 1 unknown Git LFS object:
    remote:     <file-hash>
    remote: Try to push them with 'git lfs push --all'.
    To https://github.com/<username>/tardis-regression-data.git
    ! [remote rejected] <branch-name> -> <branch-name> (pre-receive hook declined)
    error: failed to push some refs to 'https://github.com/<username>/tardis-regression-data.git'

Please check your LFS endpoint like so:

.. code-block:: shell

    $ git lfs env

And you should get:

.. code-block:: shell

    $ git lfs env
    git-lfs/3.4.1 (GitHub; linux amd64; go 1.22.2)
    git version 2.43.0
    Endpoint=https://tardis:tardis-2025-lfs@registry.moria.egr.msu.edu/repository/tardis-lfs/info/lfs (auth=basic)

If the endpoint is not this, please check your ``.git/config`` file and the ``.lfsconfig`` file. Git always prioritizes ``.git/config`` over ``.lfsconfig``, so please make sure there are no duplicates and your ``.lfsconfig`` matches with that on the TARDIS Regression Data repository. See the `git-lfs-config documentation <https://github.com/git-lfs/git-lfs/blob/main/docs/man/git-lfs-config.adoc>`_ for more details.

