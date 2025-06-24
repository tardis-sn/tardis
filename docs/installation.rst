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

Conda lockfiles are platform-specific dependency files that produce reproducible environments.
We strongly recommend installing TARDIS using this method by following the steps below.

.. note::

  Please note that you don't need to install separate environments for each TARDIS package (STARDIS, Carsus, or TARDISBase). However, for scientific reproducibility, we recommend creating a new environment whenever you start a new project with TARDIS.

1. Download the lockfile for your platform:

   .. tabs:: 

      .. group-tab:: OSX (Arm CPU)

        .. code-block:: bash

            wget -q https://github.com/tardis-sn/tardisbase/master/conda-osx-arm64.lock
        
      
      .. group-tab:: Linux

        .. code-block:: bash

            wget -q https://github.com/tardis-sn/tardisbase/master/conda-linux-64.lock
      
      .. group-tab:: OSX (Intel CPU)

        .. code-block:: bash

            wget -q https://github.com/tardis-sn/tardisbase/master/conda-osx-64.lock
        
        .. warning::
            Use at your own risk. This lockfile is not tested, so we recommend :ref:`running the test <running-tests>` before using any of the TARDIS ecosystem packages with this environment.

2. Create the environment:

Replace ``{project_name}`` with a name for your TARDIS project.

   .. tabs:: 

      .. group-tab:: OSX (Arm CPU)

        .. code-block:: bash

            conda create --name tardis-{project_name} --file conda-osx-arm64.lock
        
      
      .. group-tab:: Linux

        .. code-block:: bash

            conda create --name tardis-{project_name} --file conda-linux-64.lock
      
      .. group-tab:: OSX (Intel CPU)

        .. code-block:: bash

            conda create --name tardis-{project_name} --file conda-osx-64.lock


3. Activate the environment:

   .. code-block:: bash

       conda activate tardis-{project_name}

4. The installation process differs for developers and non-developers:

   a. Developers should `fork the repository <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ , configure
      GitHub to `work with SSH keys <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_,
      set up the `upstream remote <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/configuring-a-remote-for-a-fork>`_ and `origin` (pointing to your fork),
      and install TARDIS in development mode.

      .. code-block:: bash

        $ git clone git@github.com:tardis-sn/tardis.git
        $ cd tardis
        $ git remote add upstream git@github.com:tardis-sn/tardis.git
        $ git fetch upstream
        $ git checkout upstream/master
        $ pip install -e ".[tardisbase,viz]" # or pip install -e ".[viz]" if tardisbase is already installed in editable mode

      .. note::

        The complete developer guidelines can be found :ref:`here <developer_guidelines>`.
        
   b. Non-developers can install from specific releases using pip:

      .. code-block:: bash

        $ pip install git+https://github.com/tardis-sn/tardis.git@{tag}

      For example, to install the latest release:

      .. code-block:: bash
      
        $ pip install git+https://github.com/tardis-sn/tardis.git@release-latest

      or to install the most recent, unreleased changes from upstream:

      .. code-block:: bash

        $ pip install git+https://github.com/tardis-sn/tardis.git@master
        
    .. note::
      Running specific modules or tests will require additional optional dependencies. 
      The `tardisbase` package is required for running TARDIS Regression Tests.
      The `viz` package is required for running the TARDIS visualization tools.
      These optional dependencies can be installed by running:

      .. code-block:: bash
      
        $ pip install -e ".[tardisbase,viz]" 

      To update optional dependencies, use:

      .. code-block:: bash
      
          $ pip install -e ".[tardisbase,viz]" --upgrade --force-reinstall


5. Once finished working, you can deactivate your environment.

  ::

    $ conda deactivate

From now on, just activate the ``tardis-{project_name}`` environment before working with the TARDIS package.

You have successfully installed TARDIS! 🎉 Please refer to `Quickstart for TARDIS <quickstart.ipynb>`_ 
to start running simulations.


Environment update
==================


**Recommended approach:**
We highly recommend deleting your existing environment and creating a new one using the latest lockfile whenever you need to update your environment.

Use the following ``conda`` command to remove your current ``tardis`` environment:

.. code-block:: bash

    $ conda remove --name tardis-{project_name} --all

Now, you can create a new environment by following the steps given `here <https://tardis-sn.github.io/tardis/installation.html#install-with-lockfiles>`_.

To update the environment, download the latest lockfile and run ``conda update``.

.. code-block:: bash

    $ wget -q https://github.com/tardis-sn/tardisbase/master/conda-{platform}-64.lock
    $ conda update --name tardis --file conda-{platform}.lock

.. note::

  If you have installed TARDIS in development mode, you should *ideally* update your environment whenever you pull the latest code because the new code added might be using updated (or new) dependencies. If you don't do that and your installation seems broken, you can check if your environment requires update by comparing it against the latest environment file:

  .. code-block:: bash

      $ conda compare --name tardis-{project_name} env.yml
   
  We also recommend updating optional dependencies whenever you pull the latest code.

