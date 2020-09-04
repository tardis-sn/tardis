# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    return {
        "tardis.io.tests": [
            "data/*.dat",
            "data/*.yml",
            "data/*.csv",
            "data/*.csvy",
        ],
        "tardis.model.tests": [
            "data/*.dat",
            "data/*.yml",
            "data/*.csv",
            "data/*.csvy",
        ],
        "tardis.plasma.tests": [
            "data/*.dat",
            "data/*.yml",
            "data/*.csv",
            "data/*.txt",
        ],
    }
