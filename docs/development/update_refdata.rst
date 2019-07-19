**************************************
Procedure to update the reference data
**************************************

We asume that you have added the necessary changes to TARDIS and have a
PR open.

#. clone tardis-refdata (`git clone tardis-refdata`; you need to have git lfs installed) then make a new branch named
    the same as your new tardis feature branch.
#. Generate new reference data (in your tardis directory and right branch) offline using

    .. code-block:: None

        python setup.py test --args="--tardis-refdata=<path to refdata repo/with the right branch> --generate-reference"

#. Rerun the tests and see if it does not fail using

    .. code-block:: None
    
        python setup.py test --args="--tardis-refdata=<path to refdata repo/with the right branch>"

#. Switch to tardis-refdata. Commit the changed ref-data and open a PR on tardis-refdata
#. Make a copy of the refdata comparer notebook provided in TARDIS-refdata to check if there are
    any unexpected changes in the updated reference data and the previous reference data
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
