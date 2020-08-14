from astropy import units as u
import numpy as np
import pandas as pd
from tardis.io.decay import IsotopeAbundances
from pyne import nucname


class GenericModel:
    """A generic model that binds all the properties together,
    check the consistency of time and lengths, and allow easy
    access for the user"""

    def __init__(self, *properties, time=None):
        properties = [prop for prop in tproperties if prop is not None]
        try:
            max_property_time = time.cgs
        except (ValueError):
            print(
                "No time for the model was provided, time of properties will be"
                " used"
            )
            max_property_time = 0 * u.s
        self.names = []
        for tproperty in properties:
            setattr(self, tproperty.name, tproperty)
            self.names.append(tproperty.name)
            max_property_time = u.Quantity(
                max(max_property_time, tproperty.time)
            )

        try:
            time.to("s")
        except u.UnitsError:
            raise ValueError(
                f'"time" needs to be a quantity with units time (days, seconds,'
                f" ...) "
                + "by looking at other quantities, the earliest time you could"
                " set is"
                + "{max_property_time}"
            )

        self.number_of_cells = self.velocity.number_of_cells
        self.time = max(time, max_property_time).to("s")

        for tproperty in properties:
            tproperty.time = self.time
            tproperty.check_number_of_cells(self.number_of_cells)

        print(f"Model time set to {self.time} or {self.time.to('d')}")

    def __iter__(self):
        return iter(self.names)

    def to_dataframe(self):
        dataframe = pd.DataFrame(dtype=np.float64)
        for name in self:
            dataframe = getattr(self, name).to_dataframe(dataframe)
        dataframe.index.name = "cell"
        return dataframe


class BaseProperty:
    """Generic Class from which all properties inherit.
    Parameters
    ----------
    time_0 : astropy.units.Quantity time of the initialization
    name : str name of the property
    quantity : astropy.units.Quantity quantity to be checked
    units : str units of quantity

    Methods
    -------
    to_dataframe(dataframe) return a dataframe with all attributes in
                            dataframe_column_names. Abundances is ignored
                            for now because they are their own dataframes
    """

    def __init__(self, time_0, time=None):
        try:
            self.time_0 = time_0.to("s")
        except u.UnitConversionError:
            raise ValueError(
                '"time" needs to be a quantity with time units(days, seconds,'
                " ...)"
            )
        try:
            self.time = time.to("s")
        except AttributeError:
            # to logger
            print("no `time` provided, `time_0` will be used")
            self.time_0 = time_0.to("s")

    def to_dataframe(self, dataframe):

        names = self.dataframe_column_names
        for prop in names:
            try:
                dataframe[names[prop]] = pd.Series(getattr(self, prop).value)
            except:
                print(f"\n Warning: {prop} might be empty\n")
        return dataframe

    @u.quantity_input
    def cgs_units(self):
        print("This method is not implemented")

    def check_number_of_cells(self, number_of_cells):
        if self.number_of_cells != number_of_cells:
            raise ValueError(
                f"{self.name} number of cells is inconsistent with that of the"
                " model.                             It should be"
                f" {number_of_cells} but is {self.number_of_cells}"
            )


class Velocity(BaseProperty):
    """A class that holds the velocity arrays and derived quantities
    Parameters
    ----------
    velocity_inner : astropy.units.Quantity
    velocity_outer : astropy.units.Quantity
    time_0 : astropy.unit.Quantity

    Attributes
    ----------
    name            class name for handling purposes.

    number_of_cells return number of cells in the ejecta model.

    inner           return the inner velocity array in cgs.

    outer           return the outer velocity array in cgs.

    middle          return an array with the average of the inner and outer arrays.

    time_0          time assigned at initilization. Defaults to 0.

    time            time used for calculations. Can be updated manually or initializing GenericModel.

    radius          return array with outer radii in cgs. Uses time and velocity_outer.

    volume          return array with volume of every cells in cgs.

    dataframe       dictionary of properties that will be used in `to_dataframe()`.
    _dictionary
    """

    name = "velocity"

    def __init__(
        self, velocity_inner, velocity_outer, time_0=0 * u.d, time=None
    ):
        self._inner = self.cgs_units(velocity_inner)
        self._outer = self.cgs_units(velocity_outer)

        BaseProperty.__init__(self, time_0, time)
        self.dataframe_column_names = {
            "inner": "velocity_inner",
            "middle": "velocity_middle",
            "outer": "velocity_outer",
            "radius": "radius",
            "volume": "volume",
        }
        self.number_of_cells = len(self._inner)

    @property
    def number_of_cells(self):
        return self.__number_of_cells

    @number_of_cells.setter
    def number_of_cells(self, number_of_cells):
        if len(self._inner) != len(self._outer):
            raise ValueError(
                "cannot set number of cells, quantities lengths don't match"
            )
        self.__number_of_cells = number_of_cells

    @property
    def inner(self):
        return self._inner

    @property
    def outer(self):
        return self._outer

    @property
    def middle(self):
        return 0.5 * (self.inner + self.outer)

    @property
    def outer_radius(self):
        return (self.time - self.time_0) * self.outer

    @property
    def inner_radius(self):
        return (self.time - self.time_0) * self.inner

    @property
    def middle_radius(self):
        return (self.time - self.time_0) * self.middle

    @property
    def volume(self):
        volume = (
            4.0
            / 3
            * np.pi
            * (self.outer ** 3 - self.inner ** 3)
            * self.time ** 3
        )
        return volume.cgs

    @u.quantity_input
    def cgs_units(self, velocity: u.km / u.s):
        return velocity.cgs


class Density(BaseProperty):
    """A class that holds an initial density and time
    Parameters
    ----------

    density : astropy.units.Quantity Density at time_0

    time_0 : astropy.units.Quantity Time for initial density. Defaults to 1s to avoid `inf`

    electron :   astropy.units.Quantity electron densities

    dataframe    dictionary of properties that will be used in `to_dataframe()`.
    _dictionary

    Attributes
    ----------
    name     class name for handling purposes

    mass return density at `time`

    initial return density at `time_0`

    electron  return electron density at `time`

    Methods
    -------
    from_exponential() initializes Density class with exponential profile.

        Parameters
        ----------
        rho_0 : float64. Multiplicative constant for exponential profile
        velocity_object : Velocity object. Object from which velocities will be extracted
        exponent : exponent for the exponential profile
        electron_density : astropy.units.Quantity electron_density
        time_0 : astropy.units.Quantity time at initialization
    """

    name = "density"

    def __init__(
        self, density, electron_density=None, time_0=1 * u.s, time=None
    ):
        BaseProperty.__init__(self, time_0, time)
        self.mass_0 = self.cgs_units(density)
        if electron_density is not None:
            self.electron = self.cgs_electron(electron_density)
        self.number_of_cells = len(density)
        self.dataframe_column_names = {"mass": "density"}

    @property
    def number_of_cells(self):
        return self.__number_of_cells

    @number_of_cells.setter
    def number_of_cells(self, number_of_cells):
        if hasattr(self, "electron"):
            if len(self.mass_0) != len(self.electron):
                raise ValueError(
                    "Cannot set number of cells, quantities lengths don't match"
                )
        self.__number_of_cells = number_of_cells

    @property
    def mass(self):
        return (self.mass_0 * (self.time / self.time_0) ** -3).cgs

    @classmethod
    def from_exponential(
        cls, rho_0, velocity_object, exponent, electron_density, time_0
    ):
        density = rho_0 * np.exp(
            -(velocity_object.middle / velocity_object.inner[0])
        )
        return cls(density, electron_density, time_0)

    @u.quantity_input
    def cgs_units(self, density: u.g / u.cm ** 3):
        return density.cgs

    @u.quantity_input
    def cgs_electron(self, electron_density: 1.0 / u.cm ** 3):
        return electron_density.cgs


class Abundances(BaseProperty):
    """A class for abundances and isotope abundances.
    It takes dataframes for elemental abundance and isotopic.
    See tardis.io.decay.IsotopicAbundances for more info.

    Parameters
    ----------
    elemental : pandas.DataFrame contains elemental abundances at `time_0`

    isotope : pandas.DataFrame contains isotope abundances at `time_0`

    time_0 : astropy.units.Quantity

    Attributes:
    ----------
    elemental_0 : return DataFrame with elemental abundance at `time_0`

    isotope_0  return IsotopeAbundance object at `time_0`

    elemental return DataFrame with elemental abundances at `time`

    isotope return IsotopeAbundance object decayed at `time`
    """

    name = "abundance"

    def __init__(self, elemental, isotope, time_0=0 * u.s, time=None):
        BaseProperty.__init__(self, time_0, time)
        self.elemental_0 = elemental.sort_values(by="atomic_number", axis=0)
        self.isotope_0 = IsotopeAbundances(isotope)
        self.dataframe_column_names = {}
        self.number_of_cells = len(self.elemental.index)

    @property
    def _isotope(self):
        return self.isotope_0.decay(self.time)

    @property
    def _elemental(self):
        return self.isotope_0.decay(self.time).merge(self.elemental_0)

    @property
    def isotope(self):
        isotope_buffer = self._isotope
        atomic_index = isotope_buffer.index.get_level_values("atomic_number")
        mass_index = isotope_buffer.index.get_level_values("mass_number")
        isotope_buffer.index = self.get_isotope_index(atomic_index, mass_index)
        isotope_buffer.name = "isotope"
        isotope_buffer.columns.name = "cell"
        return isotope_buffer.T

    @property
    def elemental(self):
        elemental_buffer = self._elemental
        elemental_buffer.index = self.get_elemental_index(
            elemental_buffer.index
        )
        elemental_buffer.columns.name = "cell"
        elemental_buffer.index.name = "element"
        return elemental_buffer.T

    def get_elemental_index(self, number_index):
        names_index = []
        for number in number_index:
            names_index.append(nucname.name(number))
        return names_index

    def get_isotope_index(self, atomic_number_index, mass_number_index):
        names_index = []
        for atomic_number, mass_number in zip(
            atomic_number_index, mass_number_index
        ):
            names_index.append(
                "{0:s}{1:d}".format(nucname.name(atomic_number), mass_number)
            )
        return names_index


class RadiationField(BaseProperty):
    """A class for radiative temperature and dilution factor

    Parameters
    ----------
    radiative_tempertature : astropy.units.Quantity

    dilution_factor : numpy.ndarray

    time_0 : astropy.units.Quantity
    """

    name = "radiation_field"

    def __init__(
        self,
        radiative_temperature=None,
        dilution_factor=None,
        time_0=0 * u.d,
        time=None,
    ):
        BaseProperty.__init__(self, time_0, time)
        if radiative_temperature is not None:
            self.radiative_temperature = self.cgs_units(radiative_temperature)
        if dilution_factor is not None:
            self.dilution_factor = dilution_factor
        self.number_of_cells = self.set_number_of_cells(
            radiative_temperature, dilution_factor
        )
        self.dataframe_column_names = {
            "radiative_temperature": "radiative_temperature",
            "dilution_factor": "dilution_factor",
        }

    # FIX NUMBER OF CELLS
    @property
    def number_of_cells(self):
        return self.__number_of_cells

    @number_of_cells.setter
    def number_of_cells(self, number_of_cells):
        self.__number_of_cells = number_of_cells

    def set_number_of_cells(self, radiative_temperature, dilution_factor):
        if radiative_temperature is not None and dilution_factor is not None:
            if len(radiative_temperature) == len(dilution_factor):
                number_of_cells = len(radiative_temperature)
            else:
                raise ValueError(
                    "radiative_temperature and "
                    + "dilution_factor lengths must be the same"
                )
        elif radiative_temperature is not None:
            number_of_cells = len(radiative_temperature)
        elif dilution_factor is not None:
            number_of_cells = len(dilution_factor)
        else:
            number_of_cells = None

    @u.quantity_input
    def cgs_units(self, temperature: u.K):
        return temperature.cgs
