#setting the right include
#setting the right include
from Cython.Distutils import  Extension
import numpy as np
import os
import sys
import tempfile
import shutil
from distutils.ccompiler import new_compiler
from glob import glob
import subprocess
from astropy_helpers.setup_helpers import get_distutils_build_option


def get_extensions():
    sources = ['tardis/montecarlo/montecarlo.pyx']
    print(sources)
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src', '*.c'))]
    print(sources)
    sources += [os.path.relpath(fname) for fname in glob(
        os.path.join(os.path.dirname(__file__), 'src/randomkit', '*.c'))]
    print(sources)
    print('>>>>>>>>>>>>')
    (extra_compile_args, extra_link_args,
    compile_openmp) = _get_compile_link_args()

    return [Extension('tardis.montecarlo.montecarlo', sources,
                      include_dirs=['tardis/montecarlo/src',
                                    'tardis/montecarlo/src/randomkit',
                                    np.get_include()],
                      extra_compile_args= extra_compile_args,
                      extra_link_args= extra_link_args,
                      cython_compile_time_env = {'OPENMP': compile_openmp},
                      )]


##OpenMP

def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/7018879/
    # disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            f = open(fname, 'w')
            f.write('#include <omp.h>\n')
            f.write('#include <stdio.h>\n')
            f.write('int main(void) {\n')
            f.write('    %s();\n' % funcname)
            f.write('}\n')
            f.close()
            f = open(fname, 'r')
            f.read()
            f.close()
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            if 'CC' in os.environ:
                compiler = os.environ['CC']
                command = compiler + ' ' + fname +' -fopenmp' + ' -lgomp' '+O0'
                process = subprocess.Popen(command,shell=True, stdout=subprocess.PIPE)
                process.wait()
                if process.returncode == 0:
                    return True
                else:
                    return False
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir, extra_postargs=['-fopenmp'] )
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"), extra_postargs=['-lgomp'] )
        except:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
            shutil.rmtree(tmpdir)


def detect_openmp():
    "Does this compiler support OpenMP parallelization?"
    compiler = new_compiler()
    print "Attempting to autodetect OpenMP support... ",
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads')
    if not hasopenmp:
        print('No openMP found')
    else:
        print("OpenMP found")
    return hasopenmp


def _get_compile_link_args():
    no_openmp = get_distutils_build_option('no_openmp')
    has_openmp = detect_openmp()
    if (no_openmp is None or no_openmp==False) and has_openmp:
        extra_compile_args = ['-fopenmp']
        extra_link_args = ['-lgomp']
        compile_openmp = True
    else:
        extra_compile_args = []
        extra_link_args = []
        compile_openmp = False
    return extra_compile_args, extra_link_args, compile_openmp