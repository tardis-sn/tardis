##################
Tardis Coding Guide
##################

TARDIS follows the `Black<https://black.readthedocs.io/en/stable/>`_coding sytle for files. To install Black, do::
    pip install black
For more information on running black, please refer to `Black<https://black.readthedocs.io/en/stable/>`_. 

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
The docstring should have no leading or trailing carriage returns, and there should be a carriage return between each segment. If the function does not have a return statement, then the Returns section is not necessary.  

For the naming conventions in TARDIS, please follow the naming conventions of `PEP8 <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_. Taken from PEP8, "Function names should be lowercase, with words separated by underscores as necessary to improve readability. Variable names follow the same convention as function names. Class names should normally use the CapWords convention." For more detailed information on naming conventions, please refer to the `PEP8 <https://www.python.org/dev/peps/pep-0008/#naming-conventions>`_ guide.