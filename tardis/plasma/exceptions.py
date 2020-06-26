class PlasmaException(Exception):
    pass


class IncompleteAtomicData(PlasmaException):
    def __init__(self, atomic_data_name):
        message = (
            "The current plasma calculation requires {0}, "
            "which is not provided by the given atomic data".format(
                atomic_data_name
            )
        )
        super(PlasmaException, self).__init__(message)


class PlasmaMissingModule(PlasmaException):
    pass


class PlasmaIsolatedModule(PlasmaException):
    pass


class NotInitializedModule(PlasmaException):
    pass


class PlasmaIonizationError(PlasmaException):
    pass


class PlasmaConfigError(PlasmaException):
    pass
