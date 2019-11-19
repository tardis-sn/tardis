.. _git_lfs_instructions:

******************
Installing git-lfs
******************

**Another useful installation guide can be found at** `<https://git-lfs.github.com>`_

1. Download and install git-lfs using a package manager
    For example, using conda

    .. code-block:: None

        conda install git-lfs
    Or using MacPorts

    .. code-block:: None

        sudo port install git-lfs

2. Set up git lfs using:

    .. code-block:: None

        git lfs install


3. Download files from lfs repository (If needed)

     .. code-block:: None

        git lfs fetch upstream

4. Pull files from git lfs (if needed)

     .. code-block:: None

        git lfs pull


