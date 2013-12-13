#setting the right include

from setuptools import Extension

randomkit_files = ['tardis/randomkit/rk_isaac.c', 'tardis/randomkit/rk_mt.c', 'tardis/randomkit/rk_primitive.c',
                   'tardis/randomkit/rk_sobol.c']

def get_extensions():
    return [Extension('tardis.montecarlo',
                        ['tardis/montecarlo.pyx'] + randomkit_files)]
    #return {'tardis.montecarlo_multizone':['randomkit/*.c']}