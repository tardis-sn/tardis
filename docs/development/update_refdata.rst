**************************************
Procedure to update the reference data
**************************************

We assume that you have added the necessary changes to TARDIS and have a
PR open.

#. Fork the tardis ref-data repository like described in :ref:`development-workflow`.

#. Clone tardis-refdata. Then make a new branch named the same as your new tardis feature branch. Make sure to set the upstream of the repository to `https://github.com/tardis-sn/tardis-refdata.git`
    * You can check that you have the right upstream set using:
    .. code-block:: None

        git checkout upstream/pr/11

#. Make sure to have sure to have git-lfs installed. If not already done, the process can be seen here: :ref:`git_lfs_instructions`

#. Open tardis-refdata directory. Fetch from upstream and then use git-lfs to pull ref-data files from github using:
    .. code-block:: None

        git lfs fetch upstream
        git lfs pull

#. Open into tardis directory. Generate new reference data (in the correct branch) offline using:

    .. code-block:: None

        python setup.py test --args="--tardis-refdata=<path to refdata repo/with the right branch> --generate-reference"

#. Rerun the tests to make sure it does not fail using:

    .. code-block:: None

        python setup.py test --args="--tardis-refdata=<path to refdata repo/with the right branch>"

#. Open and run the refdata comparer notebook provided in TARDIS-refdata to check if there are any unexpected changes in the updated reference data and the previous reference data.

#. Switch to tardis-refdata. Commit the changed ref-data and open a PR on tardis-refdata

#. Switch back to the TARDIS directory. Open `.travis.yml`

#. Change the following lines

    .. code-block:: None

        - if [[ $TEST_MODE == 'spectrum' ]]; then git fetch origin pull/<your TARDIS-refdata PR number; not the TARDIS PR number>/head:<some descriptive name>; fi
        - if [[ $TEST_MODE == 'spectrum' ]]; then git checkout <some descriptive name>; fi```

#. Commit the `.travis.yml` to your Pull request
#. Make sure that your TARDIS PR now passes on TRAVIS.
#. Then merge the PR on tardis-refdata.
#. Then change `.travis.yml` to

    .. code-block:: None

        - if [[ $TEST_MODE == 'spectrum' ]]; then git fetch origin; fi
        - if [[ $TEST_MODE == 'spectrum' ]]; then git checkout origin/master; fi```
#. Then make sure that your TARDIS PR passes again.
#. Then merge your PR to TARDIS master
#. Congratulations - you updated TARDIS to be better. Have a beer and steak (or Tofu if you are vegetarian/vegan)

Troubleshooting
###############

* Unable to generate reference data
    * If generating fails due to an inability to open chianti_He.h5, make sure that you have installed git-lfs and have pulled the files from github (See steps 3 and 4).

* Error in running `comparer = ReferenceComparer(ref2_hash='upstream/pr/11')` on the comparer notebook: `No such file or directory: '.../unit_test_data.h5'`
    * If notebook file is unable to find the file /unit_test_data.h5, make sure you have correctly set your upstream. To check this, use:


    * If this fails, make sure that your upstream is set correctly to `https://github.com/tardis-sn/tardis-refdata.git` (See step 2).