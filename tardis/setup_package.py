#setting the right include

from setuptools import Extension

randomkit_files = ['tardis/randomkit/rk_isaac.c', 'tardis/randomkit/rk_mt.c', 'tardis/randomkit/rk_primitive.c',
                   'tardis/randomkit/rk_sobol.c']

def get_extensions():
    return [Extension('tardis.montecarlo',
                        ['tardis/montecarlo.pyx', 'tardis/cmontecarlo.c'] + randomkit_files, include_dirs=['tardis/randomkit']), Extension('tardis.tests.montecarlo_test_wrappers', ['tardis/tests/montecarlo_test_wrappers.pyx', 'tardis/cmontecarlo.c'] + randomkit_files, include_dirs=['tardis/randomkit'])]
    #return {'tardis.montecarlo_multizone':['randomkit/*.c']}

def get_package_data():
    return {
        'tardis.tests': ['coveragerc', 'data/*.h5']}
