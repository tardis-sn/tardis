def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_
        + ".tests": [
            "coveragerc",
            "data/*.h5",
            "data/*.dat",
            "data/*.npy",
            "integration_tests/*/*.yml",
            "integration_tests/*/*.dat",
        ]
    }
