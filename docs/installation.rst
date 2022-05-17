.. _installation:

************
Installation
************


.. warning::
    
    - TARDIS is only supported on macOS and GNU/Linux. Windows users can run TARDIS 
      from our official Docker image (*coming soon*), `WSL <https://docs.microsoft.com/en-us/windows/wsl/>`_ 
      or a Virtual Machine.

    - TARDIS packages and dependencies are distributed only through the `conda <https://docs.conda.io/en/latest/>`_ 
      package management system, therefore requires having `Anaconda <https://docs.anaconda.com/anaconda/install/index.html>`_ 
      or `Miniconda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
      installed on your system.


Development version
===================

We encourage users to use the latest TARDIS development version by reproducing the following steps:

1. Download the latest dependency file for your platform (`linux` or `osx`) from our 
   `releases page <https://github.com/tardis-sn/tardis/releases>`_.

  ::

    $ wget -q https://github.com/tardis-sn/tardis/releases/latest/download/conda-{platform}-64.lock

2. Create and activate the `tardis` environment.

  ::

    $ conda create -n tardis --file conda-{platform}-64.lock
    $ conda activate tardis

3. a. Non-developers can install the package with latest changes from upstream.

      ::

        $ pip install git+https://github.com/tardis-sn/tardis.git

   b. Instead, developers should :ref:`fork the repository <forking>`, configure
      GitHub to `work with SSH <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_,
      set up the upstream remote, and install the package in development mode.

      ::

        $ git clone git@github.com:<username>/tardis.git
        $ cd tardis
        $ git remote add upstream git@github.com:tardis-sn/tardis.git
        $ git fetch upstream
        $ git checkout upstream/master
        $ pip install -e .

      .. note::

        The complete developer guidelines can be found :ref:`here <developer_guidelines>`.


4. Once finished working, deactivate your environment.

  ::

    $ conda deactivate

You are ready! From now on, just activate the `tardis` environment before working with the 
TARDIS package.


Update an existing environment
------------------------------

To update an existing environment, download the latest environment file and run ``conda update``.

::

    $ wget -q https://github.com/tardis-sn/tardis/releases/latest/download/conda-{platform}-64.lock
    $ conda update -n tardis --file conda-{platform}-64.lock


conda-forge package
===================

TARDIS is also distributed as a `conda-forge <https://conda-forge.org/>`_ package, though we still
recommend using the latest development version for better reproducibility of scientific results.

To use this package, just create a dedicated environment with the ``tardis-sn`` package from the
``conda-forge`` channel.  

::

    $ conda create -n tardis-forge tardis-sn -c conda-forge
