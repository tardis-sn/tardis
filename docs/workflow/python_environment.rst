*****************************
Python Environment for TARDIS
*****************************


There are many ways to install the right Python packages for TARDIS. In this
cocument we outline an easy way (using `virtualenv`s) that ensures that the
installed packages are of the right version for TARDIS and not messing with
your python installation on your computer.

Using virtualenv
================

`virtualenv`_ is a tool for creating and activating isolated Python
environments that allow installing and experimenting with Python packages
without disrupting your production Python environment.  When using commands
such as ``setup.py develop``, for example, it is strong recommended to do
so within a virtualenv.

We won't provide a full tutorial on using virtualenv here--the virtualenv
documentation linked to above is a better place to start.  But here is a quick
overview on how to set up a virtualenv for TARDIS development with your
default Python version:

1. Install virtualenv::

       $ pip install virtualenv

   or::

       $ easy_install virtualenv

   or (on Debian/Ubuntu)::

       $ sudo apt-get install python-virtualenv

   etc.

2. (Recommended) Create a root directory for all your virtualenvs under a path
   you have write access to.  For example::

       $ mkdir ~/.virtualenvs

3. Create the TARDIS virtualenv::

       $ virtualenv --distribute --system-site-packages ~/.virtualenvs/tardis

   The ``--system-site-packages`` option inherits all packages already
   installed in your system site-packages directory; this frees you from having
   to reinstall packages like Numpy and Scipy in the virtualenv.  However, if
   you would like your virtualenv to use a development version of Numpy, for
   example, you can still install Numpy into the virtualenv and it will take
   precedence over the version installed in site-packages.

4. Activate the virtualenv::

       $ source ~/.virtualenvs/tardis/bin/activate

   or if you're using a csh-variant::

       $ source ~/.virtualenvs/tardis/bin/activate.csh

   virtualenv works on Windows too--see the documentation for details.

5. If the virtualenv successfully activated its name should appear in your
   shell prompt::

       (tardis) $

   The virtualenv can be disabled at any time by entering::

       (tardis) $ deactivate

6. Now as long as the virtualenv is activated packages you install with
   ``pip``, ``easy_install``, or by manually running ``python setup.py
   install`` will automatically install into your virtualenv instead of the
   system site-packages. To develop it is a useful idea to install TARDIS using
   `python setup.py develop` - which will allow you to develop and immediately
   test it out.