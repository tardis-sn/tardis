########################
Documentation Guidelines
########################


Building local docs is done via compiling docs locally. In order to do this, in the terminal type::

    python setup.py develop
    python setup.py build_docs
    
The first command only needs to be done once. This command allows you to create a link in the deployment directory so you can edit code and see changes without having to reinstall every time you make a change. 

The second command needs to be used every time you want to compile the docs to see your changes. If you don't type this command after making changes, your local docs will not recompile. Upon running this command, a new folder, _build, shuold be created inside of tardis/docs. Go through this folder and find the local docs that you are changing/creating to view how they will be compiled once a PR is approved. This method is good to check the formatting of docstrings in the API and the formatting of display files such as this one, as well as any other aspects that will be loaded in the tardis docs. For information on how docstrings should be formatted, refer `to the Tardis Coding Guide <https://tardis-sn.github.io/tarids/Code_Quality_Guidelines.html>`_ for more information and examples. 

To create new documentation that will have a link inside the `TARDIS pages <https://tardis-sn.github.io/tardis>`_, you will need to edit the corresponding index.rst file to include your file inside of the toctree. For example, when the documentation was made for `research papers <https://tardis-sn.github.io/tardis/research/research_done_using_TARDIS/research_papers.html>`_ that used TARDIS, this location's `index.rst <https://github.com/tardis-sn/tardis/blob/master/docs/research/index.rst>`_ toctree was editted to include the path to the .rst file that would be generated when the docs are compiled and create a notebook.

When making or adding to the functionality of an aspect of TARDIS, an example notebook file should be made to show the changes/function that is created using that file. If a file is being converted from one language to another, a notebook should be included to show that the same values are outputted. 

Another way to check documentation before it is submitted in a pull request is to preview the documentation online. This steps to do this is described `here <https://tardis-sn.github.io/tardis/development/documentation_preview.html>`_. This is an alternative way to build your own docs to check how they look before they are submitted if building locally is not available for any reason. 