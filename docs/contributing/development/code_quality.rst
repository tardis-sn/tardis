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

Black
-----
`Black <https://black.readthedocs.io/en/stable/index.html>`_ is a PEP 8 compliant opinionated code formatter. At TARDIS. we use Black to automatically conform to PEP 8. It is already installed in the TARDIS conda environment, so all you have to do is to run Black before commiting your changes: ::

    black {source_file_or_directory}

A better method is to run Black automatically - first `integrate it within the code editor <https://black.readthedocs.io/en/stable/editor_integration.html>`_ you use and then enable the "format on save" or "format on type" option in your editor settings.

.. warning :: If your code doesn't follow the Black code style, then the Black-check action on your PR will fail.

Ruff
----
`Ruff <https://docs.astral.sh/ruff/>`_ is a code linter and formatter that checks for common mistakes and automatically fixes them. It is currently not installed in the TARDIS conda environment, so you will have to install it manually: ::

    conda install -c conda-forge ruff

To run Ruff, use the following command: ::

    ruff check <source_file_or_directory> # Lints the code
    ruff check <source_file_or_directory> --fix # Lints and fixes any fixable errors

Currently, Ruff is not integrated with the TARDIS CI and is not a requirement for merging a PR. However, it is recommended to run Ruff on your code before commiting it to ensure that new code already follows these rules.

.. note :: We adopt the linting rules utilized by astropy. Permanent rules are defined in the ``pyproject.toml``, non-permanent rules are defined in the ``.ruff.toml`` file. If you want to add a new rule, please add it to the ``.ruff.toml`` file. If you want to add a permanent rule, please open a PR to the ``pyproject.toml``.

.. note :: Ruff can also be used for formatting code, but for now we recommend using Black for this purpose as the CI is configured to run Black on all PRs.

Pre-commit (Optional)
----
`Pre-commit <https://pre-commit.com/>`_ hooks are tools that help enforce quality standards by running checks on your code before you commit. If you choose to use pre-commit on your local machine, please follow these steps:

Install pre-commit by running: ::

    pip install pre-commit

Set up the pre-commit hooks with: ::

    pre-commit install

This needs to be done only once per repository. The pre-commit hooks will now automatically run on each commit to ensure your changes meet our code quality standards.

Naming Conventions
------------------

While Black automatically conforms your code to a majority of the PEP 8 style guide (like whitespace usage, string quotes, code layout, etc.), your code still needs to conform to the `PEP 8 naming conventions <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_. The main things to keep in mind are:

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

Here, the ``packets_mode`` parameter can only be string "virtual" or "real". Thus, in case of an invalid value passed, it raises a specific exception (``ValueError``) and informs the user of what the error is and how it can be resolved. Aside from raising exception an exception, this method safely returns without executing code with an invalid parameter value which would have caused unwanted behaviour.
