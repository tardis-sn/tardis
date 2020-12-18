##################
TARDIS Coding Guide
##################


TARDIS follows the `Black <https://black.readthedocs.io/en/stable/>`_ coding style for files. To install Black, do::

    pip install black
    
For more information on Black formatting, please refer to `Black <https://black.readthedocs.io/en/stable/>`_. 

For docstrings, TARDIS follows the docstring formatting of `Numpydocs <https://numpydoc.readthedocs.io/en/latest/format.html>`_. 
This is an example of a properly formatted docstring::
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
The docstring should have no leading or trailing carriage returns, and there should be a carriage return between each segment. In this examples the return variable and type is specified. For the Returns, the type must always be stated. It should follow the format of::
    Returns
    -------
    type
        (optional descriptor)
The return section is not necessary if the function has no return parameter. In the first example, the full path can be shown for the return time_of_model. Always list the full path for a variable type if it is not a built-in type. 

For the naming conventions in TARDIS, please follow the naming conventions of `PEP8 <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_. Taken from PEP8, "Function names should be lowercase, with words separated by underscores as necessary to improve readability. Variable names follow the same convention as function names. Class names should normally use the CapWords convention." For more detailed information on naming conventions, please refer to the `PEP8 <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_ guide. If you have questions, don't be afraid to ask. We would rather have someone ask then someone be too scared to ask. Remember, your work `is important <http://localhost:8888/files/tardis/docs/_build/html/CONTRIBUTING.html#imposter-syndrome-disclaimer>`_.