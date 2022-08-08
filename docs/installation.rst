.. _installation:

************
Installation
************


.. warning::
    
    - TARDIS is only supported on macOS and GNU/Linux. Windows users can run TARDIS 
      from our official Docker image (*coming soon*), `WSL <https://docs.microsoft.com/en-us/windows/wsl/>`_ 
      or a Virtual Machine.

    - TARDIS packages and dependencies are distributed only through the `conda <https://docs.conda.io/en/latest/>`_ 
      package management system, therefore installation requires `Anaconda <https://docs.anaconda.com/anaconda/install/index.html>`_ 
      or `Miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
      to be installed on your system.


Install with lockfiles
======================

Conda lockfiles are platform-specific dependency files that produce repeatable environments.
These files are generated on every new release. We strongly recommend installing TARDIS using
this method by following the steps described below.

1. Download the latest lockfile for your operating system from our 
   `releases section <https://github.com/tardis-sn/tardis/releases>`_, or run
   the following command while replacing ``{platform}`` with ``linux`` or ``osx`` as appropriate.

  ::

    $ wget -q https://github.com/tardis-sn/tardis/releases/latest/download/conda-{platform}-64.lock

2. Create and activate the ``tardis`` environment.

  ::

    $ conda create --name tardis --file conda-{platform}-64.lock
    $ conda activate tardis

3. a. Non-developers can install the latest release from ``conda-forge`` with the ``--no-deps`` flag,

      ::

        $ conda install tardis-sn --channel conda-forge --no-deps

      or trying the most recent, unreleased changes from upstream.

      ::

        $ pip install git+https://github.com/tardis-sn/tardis.git@master

   b. Instead, developers should `fork the repository <https://github.com/tardis-sn/tardis/fork>`_, configure
      GitHub to `work with SSH keys <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_,
      set up the `upstream remote <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/configuring-a-remote-for-a-fork>`_,
      and install the package in development mode.

      ::

        $ git clone git@github.com:<username>/tardis.git
        $ cd tardis
        $ git remote add upstream git@github.com:tardis-sn/tardis.git
        $ git fetch upstream
        $ git checkout upstream/master
        $ python setup.py develop

      .. note::

        The complete developer guidelines can be found :ref:`here <developer_guidelines>`.


4. Once finished working, you can deactivate your environment.

  ::

    $ conda deactivate

From now on, just activate the ``tardis`` environment before working with the TARDIS package.

You have successfully installed TARDIS! 🎉 Please refer to `Quickstart for TARDIS <quickstart.ipynb>`_ 
to start running simulations.


Install from package
====================

It's also possible to install TARDIS by pulling the `conda-forge package <https://anaconda.org/conda-forge/tardis-sn>`_
into a clean environment. However, we still encourage using lockfiles to ensure
reproducibility of scientific results.

::

    $ conda create --name tardis-forge tardis-sn --channel conda-forge


Environment update
==================

To update the environment after a new release, download the latest lockfile and run ``conda update``.

::

    $ conda update --name tardis --file conda-{platform}-64.lock

.. note::

  If you have installed the tardis package in development mode, you should *ideally* update your environment whenever you pull latest tardis code because the new code added might be using updated (or new) dependencies. If you don't do that and your tardis installation seems broken, you can check if your environment requires update by comparing it against the latest environment file:

  ::

      $ conda compare --name tardis tardis_env3.yml
