.. _git_lfs_instructions:

******************
Installing git-lfs
******************

**Another useful installation guide can be found at** `<https://git-lfs.github.com>`_

1. Download and install git-lfs using a package manager
    For example, using conda

    .. code-block:: bash

        conda install git-lfs
    Or using MacPorts

    .. code-block:: bash

        sudo port install git-lfs

2. Set up git lfs using:

    .. code-block:: bash

        git lfs install


3. Download files from lfs repository (If needed)

     .. code-block:: bash

        git lfs fetch upstream

4. Pull files from git lfs (if needed)

     .. code-block:: bash

        git lfs pull


