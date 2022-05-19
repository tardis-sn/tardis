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


Installing with lockfiles
=========================

We encourage all users to install TARDIS by following these steps.

1. Download the latest lockfile file for your operating system from our 
   `releases section <https://github.com/tardis-sn/tardis/releases>`_, or run
   the following command while replacing ``{platform}`` with ``linux`` or ``osx`` as appropriate.

  ::

    $ wget -q https://github.com/tardis-sn/tardis/releases/latest/download/conda-{platform}-64.lock

2. Create and activate the ``tardis`` environment.

  ::

    $ conda create --name tardis --file conda-{platform}-64.lock
    $ conda activate tardis

3. a. Non-developers can install the latest release from `conda-forge <https://anaconda.org/conda-forge/tardis-sn>`_
      with the ``--no-deps`` flag,

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
        $ pip install -e .

      .. note::

        The complete developer guidelines can be found :ref:`here <developer_guidelines>`.


4. Once finished working, you can deactivate your environment.

  ::

    $ conda deactivate

You are ready! From now on, just activate the ``tardis`` environment before working with the 
TARDIS package.


Update an existing environment
------------------------------

Just download a new environment file and run ``conda update``.

::

    $ conda update --name tardis --file conda-{platform}-64.lock


conda-forge package
===================

It's also possible to install the TARDIS environment just by pulling the `conda-forge <https://conda-forge.org/>`_
package to a dedicated environment.

.. warning::
    
    We do not recommend using this method.

::

    $ conda create --name tardis-forge tardis-sn --channel conda-forge
