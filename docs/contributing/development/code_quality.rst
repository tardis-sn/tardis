***********************
Code Quality Guidelines
***********************

Code quality ensures that new developers will have an easier time understanding what previous developers have written. Hence, in an open-source software like TARDIS, writing quality code is essential. Quoting from `this RealPython article <https://realpython.com/python-code-quality>`_, a high-quality code is identified by:

- **It does what it is supposed to do** - code should perform the functions that it is written for.

- **It does not contain defects or problems** - things shouldn't break on edge cases and defects should throw exceptions instead of causing unwanted behavior.

- **It is easy to read, maintain, and extend** - code should be easy to comprehend and it should be easy to add a new feature in it without disrupting previous features.


Code Style Conventions
======================

TARDIS follows the `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ style guide written by the author of the Python programming language. It defines a consistent way to write your code making, it easier to read and maintain.

Ruff
----
`Ruff <https://docs.astral.sh/ruff/>`_ is a code linter and formatter that checks for common mistakes and automatically fixes them. It is currently not installed in the TARDIS conda environment, so you will have to install it manually: ::

    conda install -c conda-forge ruff

To run Ruff, use the following command: ::

    ruff check <source_file_or_directory> # Lints the code
    ruff check <source_file_or_directory> --fix # Lints and fixes any fixable errors

.. note :: We adopt the linting rules utilized by astropy. Permanent rules are defined in the ``pyproject.toml``, non-permanent rules are defined in the ``.ruff.toml`` file. If you want to add a new rule, please add it to the ``.ruff.toml`` file. If you want to add a permanent rule, please open a PR to the ``pyproject.toml``.

Pre-commit (Optional)
---------------------
`Pre-commit <https://pre-commit.com/>`_ hooks are tools that help enforce quality standards by running checks on your code before you commit. If you choose to use pre-commit on your local machine, please follow these steps:

Install pre-commit by running: ::

    pip install pre-commit

Set up the pre-commit hooks with: ::

    pre-commit install

This needs to be done only once per repository. The pre-commit hooks will now automatically run on each commit to ensure your changes meet our code quality standards.

Naming Conventions
------------------

While Ruff automatically conforms your code to a majority of the PEP 8 style guide (like whitespace usage, string quotes, code layout, etc.), your code still needs to conform to the `PEP 8 naming conventions <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_. The main things to keep in mind are:

- Function names should be lowercase, with words separated by underscores as necessary to improve readability (i.e. snake_case).

- Variable names follow the same convention as function names.

- Class names should use the CapWords convention.


.. _docstrings:

Docstrings
==========

A docstring (short for documentation string) is a string that describes a module's, function's, class's, or method's definition. The docstring is a special attribute of the object (``object.__doc__``) and, for consistency, is surrounded by triple double quotes. Besides helping others understand what your code does, docstrings are used by Sphinx to auto-generate `API documentation <https://tardis-sn.github.io/tardis/api/modules.html>`_.

At TARDIS, we follow the `Numpy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`_ - please go through this format guide before writing docstrings. The following is an example of a properly formatted docstring from the `model_reader <https://github.com/tardis-sn/tardis/blob/master/tardis/io/model_reader.py>`_ module of TARDIS:

.. code-block:: python

    def read_density_file(filename, filetype):
        """
        Read different density file formats.

        Parameters
        ----------
        filename : str
            file name or path of the density file
        filetype : str
            type of the density file

        Returns
        -------
        time_of_model : astropy.units.Quantity
            time at which the model is valid
        velocity : np.ndarray
            the array containing the velocities
        unscaled_mean_densities : np.ndarray
            the array containing the densities
        """

        # Code goes here

Some of the important formatting conventions to note here are:

- The docstring should have no leading or trailing carriage returns, and there should be a carriage return between each segment.

- At the start of the docstring there is a summary explaining the purpose of the function/class/module/method. This summary should follow standard English syntax, starting with a capitalized letter and ending with appropriate punctuation.

- The docstring summary should not explain individual lines or the returns, it should summarize the purpose of the function/class/module. Comments on how individual lines work should be written using inline comments (``# comment``).

- Variable, module, function, and class names should be written between single back-ticks \` \`.

- In the above example the return variable and type is specified. For the "Returns" section, the type must always be stated, even if the variable is not. The "Returns" section should follow the format of:

.. code-block:: python

    """
    Returns
    -------
    (`optional variable name` : )type
        (optional descriptor)
    """

- The "Returns" section should not be included if the function/module/class does not have a return value(s).

- Always list the full path for a variable type if it is not a built-in type, like in above example it is shown for ``time_of_model``.


Edge Cases and Exception Handling
=================================

Code should be written with a bit of foresight to handle errors that can occur during its execution. If you know that an `exception <https://docs.python.org/3/tutorial/errors.html>`_ is likely to occur in a certain case and can be dealt with accordingly, then your code should `handle <https://docs.python.org/3/tutorial/errors.html#handling-exceptions>`_ that exception. In another scenario, you may know that a particular edge case might cause your code to break, then you should `raise <https://docs.python.org/3/tutorial/errors.html#raising-exceptions>`_ an appropriate exception to describe what has gone wrong and terminate the program's execution. An example of this in practice (taken from `here <https://github.com/tardis-sn/tardis/blob/7d7c4bc4f99c909ff45070ae9576390d96734014/tardis/widgets/kromer_plot.py#L447-L451>`_) is featured below:

.. code-block:: python

    def _calculate_plotting_data(self, packets_mode, packet_wvl_range, distance):
        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )
        # Rest of the code ...

Here, the ``packets_mode`` parameter can only be string "virtual" or "real". Thus, in case of an invalid value passed, it raises a specific exception (``ValueError``) and informs the user of what the error is and how it can be resolved. Aside from raising an exception, this method safely returns without executing code with an invalid parameter value which would have caused unwanted behaviour.

An important caveat to the above is that as research software, we expect the user
to correctly input parameters and respond to failures. Thus, we do not want to raise exceptions for every possible edge case, as this would make the code unnecessarily verbose and difficult to read. 
Instead, we should only raise exceptions for edge cases that are likely to occur and can be easily handled by the user. 
For example, in the above code snippet, it is likely that a user might accidentally pass an invalid value to ``packets_mode``, and it is easy for them to fix this error by simply passing a valid value. 
Hence, it is appropriate to raise an exception in this case. However, if there is an edge case that is unlikely to occur or cannot be easily fixed by the user, then it may be better to simply let the code fail without raising an exception.


Writing Unit Tests
==================

Unit tests are a type of software testing where individual units or components of a software are tested in isolation to ensure that they are working as expected. Writing unit tests is an important aspect of code quality as it helps to identify and fix bugs early in the development process, ensures that new changes do not break existing functionality, and provides documentation for how the code is supposed to work.

For TARDIS, we use `pytest <https://docs.pytest.org/en/7.2.x/>`_ as our testing framework. If you are adding a new feature or fixing a bug, it is important to write unit tests for your code to ensure that it works correctly and does not break existing functionality.

Unit tests should test meaningful functionality of the code, and not just test trivial things like whether a function runs without throwing an error. For example, if you are writing a function that adds two numbers, a unit test should check that the function returns the correct sum for different pairs of numbers. Use pytest parametrize to test multiple input/output pairs where individual or small input/output pairs are sufficient to test the functionality of the code.

Example:

.. code-block:: python

    @pytest.mark.parametrize(
        ["compton_opacity", "photoabsorption_opacity", "total_opacity", "expected"],
        [
            (1, 0, 1, GXPacketStatus.COMPTON_SCATTER),
            (0, 1, 1, GXPacketStatus.PHOTOABSORPTION),
            (0, 0, 1, GXPacketStatus.PAIR_CREATION),
        ],
    )
    def test_scatter_type(
        compton_opacity, photoabsorption_opacity, total_opacity, expected
    ):
        """Test the scattering type

        Parameters
        ----------
        compton_opacity : float
        photoabsorption_opacity : float
        total_opacity : float
        expected : list
            Expected parameters
        """
        actual = scatter_type(compton_opacity, photoabsorption_opacity, total_opacity)
        assert actual == expected


Tests should make use of pytest fixtures where possible to avoid code duplication and make the tests easier to read and maintain. For example, if you are testing a function that requires a certain input data structure, you can create a fixture that sets up that data structure and use it in your tests. Use existing fixture where possible to reduce code duplication and make the tests faster to run. Search for fixtures that are already defined using your IDE.

Example:

.. code-block:: python

    @pytest.fixture(scope="session")
    def simulation_verysimple_vpacket_tracking(config_verysimple, atomic_dataset):
        atomic_data = deepcopy(atomic_dataset)
        sim = Simulation.from_config(
            config_verysimple, atom_data=atomic_data, virtual_packet_logging=True
        )
        sim.last_no_of_packets = 4000
        sim.run_final()
        return sim


Using Regression Data
=====================

When writing tests, it is often helpful to use regression data to ensure that your code produces the expected results. Regression data is a set of output values that are used to test the functionality of a piece of code. By using regression data, you can ensure that your code produces the same output for the same input, which helps to catch any unintended changes or bugs that may have been introduced during development.

TARDIS has a `regression data repository <https://github.com/tardis-sn/tardis-regression-data>`_ that contains a variety of output data for different features of the code. When writing unit tests, you can use this regression data to test your code against known inputs and outputs. 

We use a regression data class defined in `tardisbase <https://github.com/tardis-sn/tardisbase>`_ to access the regression data in our unit tests. This class provides methods to load the regression data from numpy arrays or pandas dataframes stored as HDF files.

.. note:: Try to re-use existing regression data where possible so that we minimize the amount of regression data storage required, and reduce duplication of data and thus confusion.

Example:

.. code-block:: python  
    
    def test_transition_probabilities_regression(
        self, continuum_macro_atom_state, regression_data
    ):
        """Test transition probabilities using regression data.

        This test stores the transition probabilities in regression data
        for comparison across runs.
        """
        actual = continuum_macro_atom_state.transition_probabilities
        expected = regression_data.sync_dataframe(actual)
        pdt.assert_frame_equal(actual, expected)

Note that ``regression_data`` is loaded as a fixture. The ``sync_dataframe`` method is used to synchronize the actual dataframe with the expected dataframe from the regression data, ensuring that they have the same structure and index before comparison. The `assert_frame_equal` function from `pandas.testing` is then used to compare the two dataframes and assert that they are equal.

The ``sync_dataframe`` method allows the regression data to be updated or read
depending on the pytest arguments used. ``pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference`` will update the regression data with the actual dataframe, while ``pytest tardis --tardis-regression-data=/path/to/tardis-regression-data`` will read the expected dataframe from the regression data for comparison. For more information on how to update the regression data, please refer to the `Update the Regression Data <https://tardis-sn.github.io/tardis/contributing/development/update_regression_data.html>`_ section of the documentation.

.. important:: Updating the regression data should be done with caution and only when necessary, as it can potentially mask underlying issues in the code. Always ensure that the changes you are making to the regression data are intentional and justified, and that they do not introduce any unintended consequences or bugs into the code. It is safest to add new regression data rather than update existing regression data.