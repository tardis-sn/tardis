from astropy.units import Quantity


class IonizationData:
    def __init__(self, ionization_data):
        # Convert ionization energies to CGS
        ionization_data = ionization_data.squeeze()
        ionization_data[:] = Quantity(ionization_data[:], "eV").cgs.value
        self.data = ionization_data

    def active_data(self):
        pass


class ZetaData:
    def __init__(self, zeta_data):
        self.data = zeta_data

    def active_data(self):
        pass


class PhotoIonizationData:
    def __init__(self, photoionization_data):
        self.data = photoionization_data

    def active_data(self):
        pass
