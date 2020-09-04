**************************************
Procedure to update the reference data
**************************************

We assume that you have added the necessary changes to TARDIS and have an open pull request.

#. Fork the ``tardis-refdata`` repository as described in :ref:`development-workflow`.

#. Clone ``tardis-refdata`` to your computer. Then make a new branch named the same as your new TARDIS feature branch. Make sure to have correctly set up the ``upstream`` remote as explained in :ref:`reviewing-and-helping-with-pr`. 

#. Activate your ``tardis`` conda environment.

    .. code-block:: None
    
        conda activate tardis

#. Navigate to your ``tardis-refdata`` directory and type:

    .. code-block:: None
    
        git lfs install
        git lfs fetch upstream
        git lfs pull

#. Go to your ``tardis`` directory. Make sure you are working on the correct branch and generate new reference data using:

    .. code-block:: None

        python setup.py test --args="--tardis-refdata=<path to refdata repo on the right branch> --generate-reference"

#. Re-run the tests to make sure it does not fail using:

    .. code-block:: None

        python setup.py test --args="--tardis-refdata=<path to refdata repo on the right branch>"

#. Go back to your ``tardis-refdata`` folder, and run a Jupyter Notebook session inside the ``notebook`` folder. 

#. Open the ``ref_data_compare.ipynb`` file notebook and look for the cell with the following code:

    .. code-block:: None

        comparer = ReferenceComparer(ref2_hash='upstream/pr/24')
        
   Replace '24' for the number of the last merged pull request showed `here <https://github.com/tardis-sn/tardis-refdata/pulls?utf8=%E2%9C%93&q=is%3Apr+is%3Aclosed>`_.

#. Run all cells and check if there are any unexpected changes in the updated reference data respect the previous one.

#. Commit the changes made to the reference data and open a PR on ``tardis-refdata``.

#. Switch back to the ``tardis`` directory. Open the ``.travis.yml`` file and change the following lines:

    .. code-block:: None

        - if [[ $TEST_MODE == 'spectrum' ]]; then git fetch origin pull/<your tardis-refdata PR number; not the TARDIS PR number>/head:<some descriptive name>; fi
        - if [[ $TEST_MODE == 'spectrum' ]]; then git checkout <some descriptive name>; fi```

#. Commit the ``.travis.yml`` to your pull request.

#. Make sure that your TARDIS pull request now passes all tests on Travis-CI.

#. Ask for review for your PR on ``tardis-refdata`` and wait for merge.

#. Then change the ``.travis.yml`` in ``tardis`` directory to:

    .. code-block:: None

        - if [[ $TEST_MODE == 'spectrum' ]]; then git fetch origin; fi
        - if [[ $TEST_MODE == 'spectrum' ]]; then git checkout origin/master; fi```

#. Ensure TARDIS pull request passes Travis-CI again and ping someone to merge your PR to the TARDIS master branch.


Congratulations! You have updated TARDIS to be better. Have a beer and steak (or Tofu if you are vegetarian/vegan).


Troubleshooting
###############

* Unable to generate reference data
    * If generating fails due to an inability to open ``chianti_He.h5``, make sure that you have activated your `tardis` conda environment and that ``git-lfs`` is installed. Fetch and pull the files from GitHub as explained in step 4.

* Error when running ``comparer = ReferenceComparer(ref2_hash='upstream/pr/XX')`` on the comparer notebook: ``No such file or directory: '.../unit_test_data.h5'``
    * If notebook file is unable to find the file ``unit_test_data.h5``, make sure you have correctly set your upstream as explained in :ref:`reviewing-and-helping-with-pr`.
