Building local docs is done through compiling docs locally. In order to do this, go into the tardis/docs/conf.py file. Create a new line before the line::
    # -- General configuration ----------------------------------------------------
and set it to be::
    nbsphinx_allow_errors = True
    # -- General configuration ----------------------------------------------------
Once this is done, in the terminal inside the main tardis folder, type "python setup.py build_sphinx". This should create a new folder, _build, inside tardis/docs. Go through this folder and find the local docs that you are changing/creating to view how they will be compiled once a PR is approved. This method is good to check the formatting of docstrings in the API, as well as any other aspects that will be loaded in the tardis docs. 

For information on how docstrings should be formatted, refer to 'here <link for the coding_guide.rst file>`_. If you are creating new documentation that will have a link inside the `TARDIS pages <https://tardis-sn.github.io/tardis>`_, you will need to edit the locations index.rst file to include your own .rst file. For example, when the documentation was made for papers, the `file <https://github.com/tardis-sn/tardis/blobl/master/docs/research/index.rst>`_ was editted to include the path to the .rst file that would be generated when the docs are compiled and create a notebook.

When making or adding to the functionality of an aspect of TARDIS, an example notebook file should be made to show the changes/function that is created using that file. If a file is being converted from one language to another, a notebook should be included to show that the same values are outputted. 