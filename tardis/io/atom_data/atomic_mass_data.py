from tardis import constants as const


class AtomicMassData:
    def __init__(self, atomic_mass_data):
        # Convert atomic masses to CGS
        # We have to use constants.u because astropy uses
        # different values for the unit u and the constant.
        # This is changed in later versions of astropy (
        # the value of constants.u is used in all cases)
        atomic_mass_data.loc[:, "mass"] = (
            atomic_mass_data["mass"].values * const.u.cgs.value
        )
        self.data = atomic_mass_data
