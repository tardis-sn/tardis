******************
Testing Guidelines
******************

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

.. note:: Try to reuse existing regression data where possible so that we minimize the amount of regression data storage required, and reduce duplication of data and thus confusion.

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
depending on the pytest arguments used. ``pytest tardis --tardis-regression-data=/path/to/tardis-regression-data --generate-reference`` will update the regression data with the actual dataframe, while ``pytest tardis --tardis-regression-data=/path/to/tardis-regression-data`` will read the expected dataframe from the regression data for comparison. for more information on how to update the regression data, please refer to the `Update the Regression Data <https://tardis-sn.github.io/tardis/contributing/development/update_regression_data.html>`_ section of the documentation.

.. important:: Updating the regression data should be done with caution and only when necessary, as it can potentially mask underlying issues in the code. Always ensure that the changes you are making to the regression data are intentional and justified, and that they do not introduce any unintended consequences or bugs into the code. It is safest to add new regression data rather than update existing regression data.