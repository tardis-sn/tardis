**************************
Outdated Installation FAQs
**************************

We highly encourage with any installation problems to try the recommended install
method because this often fixes problems. Here are some common problems when
installing and their fixes:

**Problem:** While building TARDIS via ``python setup.py`` build you
may encounter the following error::

    error: tardis/montecarlo/montecarlo.c: Could not find C file tardis/montecarlo/montecarlo.c for Cython file tardis/montecarlo/montecarlo.pyx when building extension tardis.montecarlo.montecarlo. Cython must be installed to build from a git checkout.


**Solution:** There are several solutions to this problem. A clean checkout will
help. To clean up your repository please try ``python setup.py clean`` and
then ``git clean -dfx`` (**WARNING** will delete any non-TARDIS file in that directory)
This will often clean this problem. If it still persists:

Go into the tardis/montecarlo directory and build montecarlo.c by hand::

    cython montecarlo.pyx

Then, ``python setup.py build`` should run without problems.


**Problem:** when trying to set up CC=gcc python setup.py develop --with-openmp the following error popped up: 
from tardis/_compiler.c:1: /Users/yssavo/miniconda2/envs/tardis-show2/lib/gcc/x86_64-apple-darwin13.4.0/5.2.0/include-fixed/limits.h:168:61: fatal error: limits.h: No such file or directory 
        
**Solution:** Run on terminal: 

    open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg

**Problem:** Symbol not found: _GOMP_parallel when compiling with `--with-openmp`

**Solution:** Install gcc8 from macports and then install with these flags: `link_args = ['-fopenmp','-Wl,-rpath,/opt/local/lib/gcc8/']`

**Problem:** While building TARDIS (via python 2.7) via ``python setup.py`` build you
may encounter the following error::

     TypeError: super() argument 1 must be type, not None
    
    ----------------------------------------
    Command "python setup.py egg_info" failed with error code 1 in /tmp/pip-req-build-wPB39p/


**Solution:** The cause for this problem is Sphinx or Sphinx version. It can be easily solved by installing Sphinx 1.5.6.
              The command for the same is :

    pip install sphinx==1.5.6
    
    or
    
    conda install sphinx==1.5.6

Then, ``python setup.py build install`` should run without problems.
