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
We strongly recommend installing TARDIS ecosystem packages using this method by following the steps below.

1. Download the lockfile for your platform:

   .. code-block:: bash

       wget -q https://github.com/tardis-sn/tardisbase/master/conda-{platform}-64.lock

   Replace ``{platform}`` with ``linux-64`` or ``osx-arm64`` based on your operating system.

2. Create the environment:

   .. code-block:: bash

       conda create --name tardis --file conda-{platform}.lock

3. Activate the environment:

   .. code-block:: bash

       conda activate tardis

4. To install TARDIS, first execute these commands:

   .. code-block:: bash

      $ git clone git@github.com:tardis-sn/tardis.git
      $ cd tardis
      $ git remote add upstream git@github.com:tardis-sn/tardis.git
      $ git fetch upstream
      $ git checkout upstream/master
    
   The installation process differs for developers and non-developers:

   a. Developers should `fork the repository <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ of the package to be installed, configure
      GitHub to `work with SSH keys <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_,
      set up the `upstream remote <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/configuring-a-remote-for-a-fork>`_,
      and install the package in development mode.

      .. code-block:: bash

        $ pip install -e .

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
      Running specific modules or tests for some packages might require additional optional dependencies. 
      The tardisbase package can also be installed as an optional dependency.
      These optional dependencies can be installed by running:
      
      .. code-block:: bash
      
        $ pip install -e ".[optional_dependencies]"
        # for example:
        # pip install -e ".[tardisbase]" # installs the package with tardisbase optional dependency group
        # for multiple optional dependencies
        # $ pip install -e ".[tardisbase,viz]" # installs the package with tardisbase and viz optional dependency groups
      
      .. note::
        The tardisbase package is required for running TARDIS Regression Tests. The viz package is required for running the TARDIS visualization tools.

      To update optional dependencies, use:

      .. code-block:: bash
      
          $ pip install -e ".[optional_dependency]" --upgrade --force-reinstall
          # for example:
          $ pip install -e ".[tardisbase]" --upgrade --force-reinstall # forces reinstall of tardisbase dependencies group
      
      See the package documentation for a complete list of optional dependencies for that package.

5. Once finished working, you can deactivate your environment.

  ::

    $ conda deactivate

From now on, just activate the ``tardis`` environment before working with the TARDIS package.

You have successfully installed TARDIS! 🎉 Please refer to `Quickstart for TARDIS <quickstart.ipynb>`_ 
to start running simulations.


.. Install from package
.. ====================

.. It's also possible to install TARDIS by pulling the `conda-forge package <https://anaconda.org/conda-forge/tardis-sn>`_
.. into a clean environment. However, we still encourage using lockfiles to ensure
.. reproducibility of scientific results.

.. ::

..     $ conda create --name tardis-forge tardis-sn --channel conda-forge


Environment update
==================

To update the environment, download the latest lockfile and run ``conda update``.

.. code-block:: bash

    $ wget -q https://github.com/tardis-sn/tardisbase/master/conda-{platform}-64.lock
    $ conda update --name tardis --file conda-{platform}.lock

.. note::

  If you have installed the package in development mode, you should *ideally* update your environment whenever you pull latest package code because the new code added might be using updated (or new) dependencies. If you don't do that and your installation seems broken, you can check if your environment requires update by comparing it against the latest environment file:

  .. code-block:: bash

      $ conda compare --name tardis env.yml
   
  We also recommend updating optional dependencies whenever you pull latest code.


**Recommended approach:**

We highly recommend deleting your existing environment and creating a new one using the latest lockfile whenever you need to update your environment.

Use the following ``conda`` command to remove your current ``tardis`` environment:

.. code-block:: bash

    $ conda remove --name tardis --all

Now, you can create a new environment by following the steps given `here <https://tardis-sn.github.io/tardis/installation.html#install-with-lockfiles>`_.

