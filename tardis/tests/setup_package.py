def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.tests': ['coveragerc', 'data/*.h5',
                                            'data/*.dat', 'data/*.npy',
                                            'tests_slow/w7/*.yml',
                                            'tests_slow/w7/*.dat',
                                            'tests_slow/at/*.yml',
                                            'tests_slow/at/*.dat']}
