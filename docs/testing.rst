Testing TARDIS
--------------

Testing a code like TARDIS is an important step to making sure that the output is trustworthy. The TARDIS team has opted
to use the py.test testing framework for use with TARDIS (and some astropy extensions of it). To run the tests download
the desired version of TARDIS and run it with::

    python setup.py test

This will run the currently implemented tests (which are not as many as there should be)

To quickly test TARDIS with any atomic dataset (it will only see if it in general runs)::

    python setup.py test --args="--atomic-dataset=<path_to_dataset> -s"

