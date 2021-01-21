***********************
Code Quality Guidelines
***********************

Code quality ensures that new developers will have an easier time understanding what previous developers have written. Hence, in an open-source software like TARDIS, writing quality code is very essential. Quoting from `this RealPython article <https://realpython.com/python-code-quality>`_, a high-quality code is identified by:

- **It does what it is supposed to do** - code should obviously perform the functionality it is written for.

- **It does not contain defects or problems** - things shouldn't break on edge cases and defects should throw exceptions instead of causing unwanted behavior.

- **It is easy to read, maintain, and extend** - code should be easy to comprehend and it should be easy to add a new feature in it without disrupting previous features.


Code Style Conventions
======================

TARDIS follows the `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ - style guide written by the author of the Python programming language itself. It defines a consistent way to write your code making it easier to read and maintain.

Black
-----
`Black <https://black.readthedocs.io/en/stable/index.html>`_ is a PEP 8 compliant opinionated code formatter. At TARDIS, we use Black that automatically takes care of PEP 8 compilance. It is already installed in the TARDIS conda environment you have created, so all you have to do is to run black before commiting your changes: ::
    
    black {source_file_or_directory}

A better way instead is to run black automatically - first `integrate it within the code editor <https://black.readthedocs.io/en/stable/editor_integration.html>`_ you use and then enable "format on save" or "format on type" option in your editor settings.

.. warning :: If your code doesn't follow black code style, black-check action on your PR will fail.

Naming Conventions
------------------

While black will take care of the most aspects of PEP 8 (like whitespace usage, string quotes, code layout, etc.), you still yourself need to take care of the `PEP 8 naming conventions <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_. The main things here to keep in mind are:

- Function names should be lowercase, with words separated by underscores as necessary to improve readability (i.e. snake_case).

- Variable names follow the same convention as function names. 

- Class names should normally use the CapWords convention.


.. _docstrings:

Docstrings
==========

Docstring (short for documentation string) is a string that describes a module, function, class, or method definition. The docstring is a special attribute of the object (``object.__doc__``) and, for consistency, is surrounded by triple double quotes. Besides helping others understand what your code does, docstrings are used by Sphinx to auto-generate `API documentation <https://tardis-sn.github.io/tardis/api/modules.html>`_.

At TARDIS, we use `Numpy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`_ - please go through this format guide before writing docstrings. Following is an example of a properly formatted docstring from `model_reader <https://github.com/tardis-sn/tardis/blob/master/tardis/io/model_reader.py>`_ module of TARDIS:

.. code-block:: python

    def read_density_file(filename, filetype):
        """
        read different density file formats

        Parameters
        ----------
        filename : str
            filename or path of the density file
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

Some of the important formatting conventions to note here, are:

- The docstring should have no leading or trailing carriage returns, and there should be a carriage return between each segment. 

- At the start of the docstring there is a summary explaining the purpose of the function/class/module. This summary should follow the standard English syntax, starting with a capitalized letter and ending with appropriate punctuation.

- The docstring summary should not explain individual lines or the returns, it should summarize the purpose of the function/class/module. Comments on how individual lines work should be written using inline comments (``# comment``).

- Variable, module, function, and class names should be written between single back-ticks \` \`.

- In the above example the return variable and type is specified. For the "Returns" section, the type must always be stated, even if the variable is not. The "Returns" section should follow the format of: ::

    Returns
    -------
    (`optional variable name` : )type
        (optional descriptor)

- The "Returns" section should not be included if the function/module/class does not have a return value(s).

- Always list the full path for a variable type if it is not a built-in type, like in above example it is shown for ``time_of_model``.


Edge Cases and Exception Handling
=================================

Code should be written with a bit of foresight to handle errors that can occur during the execution. If you know that an `exception <https://docs.python.org/3/tutorial/errors.html>`_ is likely to occur in a certain case and can be dealt accordingly, then your code should `handle <https://docs.python.org/3/tutorial/errors.html#handling-exceptions>`_ that exception. In another scenario, you may know that a particular edge case might cause your code to break, then you should `raise <https://docs.python.org/3/tutorial/errors.html#raising-exceptions>`_ an appropriate exception to inform what has gone wrong and to terminate the program execution. An example of this in practice (taken from `here <https://github.com/tardis-sn/tardis/blob/7d7c4bc4f99c909ff45070ae9576390d96734014/tardis/widgets/kromer_plot.py#L447-L451>`_):

.. code-block:: python

    def _calculate_plotting_data(self, packets_mode, packet_wvl_range, distance):
        if packets_mode not in ["virtual", "real"]:
            raise ValueError(
                "Invalid value passed to packets_mode. Only "
                "allowed values are 'virtual' or 'real'"
            )
        # Rest of the code ...

Here ``packets_mode`` parameter can only be string "virtual" or "real". Thus in case of an invalid value passed, it raises a specific exception (``ValueError``) and informs the user about what is the error and how it can be resolved. Besides, by raising exception, the method safely returns without executing code with wrong parameter value which would have caused unwanted behaviour.
