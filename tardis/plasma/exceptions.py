class PlasmaException(Exception):
    pass


class IncompleteAtomicData(PlasmaException):
    def __init__(self, atomic_data_name):
        message = (
<<<<<<< HEAD
            "The current plasma calculation requires {0}, "
            "which is not provided by the given atomic data".format(
                atomic_data_name
            )
=======
            f"The current plasma calculation requires {atomic_data_name}, "
            f"which is not provided by the given atomic data"
>>>>>>> 5c7f60f3... all string formatting done
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
