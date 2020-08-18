from astropy import units as u
import numpy as np
from tardis.io.decay import IsotopeAbundances
from pyne import nucname


class GenericModel:
    """A generic ejecta model that uses all properties.

    The model is intended to be used with BaseProperty children classes.
    Using the time of each property it will sync all properties and
    check the number of cells in each one.
    Allows easy access and modification of quantities
    contained in the properties.
    """

    def __init__(self, *properties, time=None):
        properties = [prop for prop in properties if prop is not None]
        try:
            max_property_time = time.cgs
        except AttributeError:
            print(
                "No time for the model was provided, "
                + "time of properties will be used"
            )
            max_property_time = 0 * u.s
            time = 0 * u.s
        except u.UnitsError:
            raise ValueError(
                '"time" needs to be a quantity with units time (days, seconds,'
                " ...) "
            )

        self.names = []
        for tproperty in properties:
            setattr(self, tproperty.name, tproperty)
            self.names.append(tproperty.name)
            if tproperty.time is not None:
                max_property_time = max(max_property_time, tproperty.time)

        self.number_of_cells = self.velocity.number_of_cells
        self.time = max(time, max_property_time)

        for tproperty in properties:
            tproperty.time = self.time
            tproperty.check_number_of_cells(self.number_of_cells)

        print(f"Model time set to {self.time} or {self.time.to('d')}")

    def __iter__(self):
        return iter(self.names)

    def to_dataframe(self):
        """Not implemented"""


class BaseProperty:
    """Generic Class from which all properties inherit.

    Attributes
    ----------
    time_0 : astropy.units.Quantity time of the initialization
    time : astropy.units.Quantity current time

    """

    def __init__(self, time_0, time=None):
        """Initialize `BaseProperty`.

        It is intended to handle the initialization of time and initial time.
        If no time is provided, it will fallback to time_0.

        Parameters:
        time_0: astropy.units.Quantity initial time of property
        time: astropy.units.Quantity current time of property. Default None.
        """

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
            print("no `time` provided, `time_0` will be used")
            self.time = time_0.to("s")

    @u.quantity_input
    def cgs_units(self):
        print("This method is not implemented")

    def check_number_of_cells(self, number_of_cells):
        """Check Property's number of cells."""

        if self.number_of_cells != number_of_cells:
            raise ValueError(
                f"{self.name} number of cells is inconsistent with that of the"
                " model.                             It should be"
                f" {number_of_cells} but is {self.number_of_cells}"
            )


class Velocity(BaseProperty):
    """A class that holds the velocity arrays and derived quantities.

    The Velocity Class is the most important to define an ejecta model.
    It provides the grid upon which all other Properties are referenced to.
    Inner and outer velocity arrays must have the same length and have the
    same internal values (velocity_inner[1:] == velocity_outer[:-1])

    Attributes
    ----------
    name        class name for handling purposes.

    number      int.
    _of_cells   Return number of cells in the ejecta model.

    inner       astropy.units.Quantity.
                Return the inner velocity array in cgs.

    outer       astropy.units.Quantity.
                Return the outer velocity array in cgs.

    middle      astropy.units.Quantity.
                Return an array with the average of the inner and outer arrays.

    time_0      astropy.units.Quantity.
                Time assigned at initilization. Default 0.

    time        astropy.units.Quantity. Time used for calculations.
                Can be updated manually or initializing GenericModel.

    radius      return array with outer radii in cgs.
                Uses time and velocity_outer.

    volume      return array with volume of every cells in cgs.

    """

    name = "velocity"

    def __init__(
        self, velocity_inner, velocity_outer, time_0=0 * u.d, time=None
    ):
        """Initialize Velocity class. Extend BaseProperty.

        Using inner velocity and outer velocity array check their lengths,
        and units. Expands `BaseProperty`.

        Arguments:
        ---------
        velocity_inner: astropy.units.Quantity. Cells' inner velocity array.
        velocity_outer: astropy.units.Quantity. Cells' outer velocity array.
        time_0: astropy.units.Quantity. Initial time of class. Default 0 day.
        time: astropy.units.Quantity. Current time of class. Default None.
        """

        self._inner = self.cgs_units(velocity_inner)
        self._outer = self.cgs_units(velocity_outer)

        BaseProperty.__init__(self, time_0, time)
        self.number_of_cells = len(self._inner)

    @property
    def number_of_cells(self):
        return self.__number_of_cells

    @number_of_cells.setter
    def number_of_cells(self, number_of_cells):
        """Check if arrays length match and return number of cells."""
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
        """Check input velocities units and retun CGS arrays"""
        return velocity.cgs


class Density(BaseProperty):
    """A class that holds an initial mass and electron density and holds
     density related quantities.

    Attributes
    ----------
    name:     str.
              Class name for handling purposes

    mass:     astropy.units.Quantity.
              Return density at current time

    initial:  astropy.units.Quantity.
              Return density at initial time

    electron: astropy.units.Quantity.
              Return electron density at `time`

    time_0:   astropy.units.Quantity.
              Time assigned at initilization. Default 0.

    time:     astropy.units.Quantity. Time used for calculations.
              Can be updated manually or initializing GenericModel.

    Methods
    -------
    from_exponential() initializes Density class with exponential profile.
    Return Density object.
    """

    name = "density"

    def __init__(
        self, density, electron_density=None, time_0=1 * u.s, time=None
    ):
        """Initialize Density class. Extend BaseProperty.

        Density array(and electron density if provided) are checked for units,
        and number of cells is defined by their lengths.

        Arguments:
        ---------
        density:  astropy.units.Quantity.Initial mass density array.
        electron
        _density: astropy.units.Quantity. Initial electron density array.
        time_0:   astropy.units.Quantity.
                  Initial time of Density. Default 1 second.
        time:     astropy.units.Quantity. Current time of Density.
                  Default None.
        """
        BaseProperty.__init__(self, time_0, time)
        self.mass_0 = self.cgs_units(density)
        if electron_density is not None:
            self.electron = self.cgs_electron(electron_density)
        self.number_of_cells = len(density)

    @property
    def number_of_cells(self):
        return self.__number_of_cells

    @number_of_cells.setter
    def number_of_cells(self, number_of_cells):
        if hasattr(self, "electron"):
            if len(self.mass_0) != len(self.electron):
                raise ValueError(
                    "Cannot set number of cells, quantities lengths don't"
                    " match."
                )
        self.__number_of_cells = number_of_cells

    @property
    def mass(self):
        """Return mass density of cells."""
        return (self.mass_0 * (self.time / self.time_0) ** -3).cgs

    @classmethod
    def from_exponential(
        cls,
        rho_0,
        velocity_middle,
        velocity_inner,
        exponent,
        electron_density,
        time_0,
        time,
    ):
        """Instantiate a Density Class with and exponential profile.

        Using the velocity array it generates an exponential density profile.

        Argument
        ----------
        rho_0:           float64.
                         Multiplicative constant for exponential profile.

        velocity_middle: astropy.unit.Quantity.
                         Middle velocity will be used to calculate the profile.

        velocity_inner:  astropy.unit.Quantity.
                         Only the first value is used to calculate the profile.

        exponent :       int.
                         exponent for the exponential profile.

        electron_density: astropy.units.Quantity electron_density.

        time_0 :         astropy.units.Quantity time at initialization.

        time:            astropy.units.Quantity. Current time of Density.
                         Default None.

        """

        density = rho_0 * np.exp(-(velocity_middle / velocity_inner[0]))
        return cls(density, electron_density, time_0=1 * u.s, time=None)

    @u.quantity_input
    def cgs_units(self, density: u.g / u.cm ** 3):
        """Check input mass density units and return CGS arrays"""

        return density.cgs

    @u.quantity_input
    def cgs_electron(self, electron_density: 1.0 / u.cm ** 3):
        """Check input electron density units and return CGS arrays"""

        return electron_density.cgs


class Abundances(BaseProperty):
    """A class for abundances and isotope abundances.

    It takes dataframes for elemental abundance and isotopic and
    instantiate a IsotopeAbundance object.
    See tardis.io.decay.IsotopicAbundances for more info.


    Attributes:
    ----------
    name:         str. Class name for handling purposes.

    elemental_0 : return DataFrame with elemental abundance at `time_0`.

    isotope_0:    return IsotopeAbundance object at `time_0`.

    elemental:    return DataFrame with elemental abundances at `time`

    isotope:      return IsotopeAbundance object decayed at `time`.

    time_0:       astropy.units.Quantity.
                  Time assigned at initilization. Default 0.

    time:         astropy.units.Quantity. Time used for calculations.
                  Can be updated manually or initializing GenericModel.
    """

    name = "abundance"

    def __init__(self, elemental, isotope, time_0=0 * u.s, time=None):
        """Instantiate Abundaces Class using IsotopeAbundance class.

        Arguments
        ----------
        elemental: pandas.DataFrame contains elemental abundances at `time_0`
        isotope: pandas.DataFrame contains isotope abundances at `time_0`
        time_0: astropy.units.Quantity.
        time: astropy.units.Quantity.
        """
        BaseProperty.__init__(self, time_0, time)
        self.elemental_0 = elemental.sort_values(by="atomic_number", axis=0)
        self.isotope_0 = IsotopeAbundances(isotope)
        self.number_of_cells = len(self.elemental.index)

    @property
    def _isotope(self):
        return self.isotope_0.decay(self.time)

    @property
    def _elemental(self):
        return self.isotope_0.decay(self.time).merge(self.elemental_0)

    @property
    def isotope(self):
        """Return dataframe of isotope abundances decayed by time - time_0"""
        isotope_buffer = self._isotope
        atomic_index = isotope_buffer.index.get_level_values("atomic_number")
        mass_index = isotope_buffer.index.get_level_values("mass_number")
        isotope_buffer.index = self.get_isotope_index(atomic_index, mass_index)
        isotope_buffer.name = "isotope"
        isotope_buffer.columns.name = "cell"
        return isotope_buffer.T

    @property
    def elemental(self):
        """Return dataframe of abundances decayed by time - time_0"""
        elemental_buffer = self._elemental
        elemental_buffer.index = self.get_elemental_index(
            elemental_buffer.index
        )
        elemental_buffer.columns.name = "cell"
        elemental_buffer.index.name = "element"
        return elemental_buffer.T

    def get_elemental_index(self, number_index):
        """Return list of element names from their atomic number"""
        names_index = []
        for number in number_index:
            names_index.append(nucname.name(number))
        return names_index

    def get_isotope_index(self, atomic_number_index, mass_number_index):
        """Return list of element names from their atomic number
        and mass number with format.
        e.g. Ni56 for isotope 56 of Nickel
        """

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
    radiative_tempertature: astropy.units.Quantity.
                            Radiative temperature of cells.

    dilution_factor: numpy.ndarray. Dilution factoy of cells

    time_0:          astropy.units.Quantity.
                     Time assigned at initilization. Default 0.

    time:            astropy.units.Quantity. Time used for calculations.
                     Can be updated manually or initializing GenericModel.
    """

    name = "radiation_field"

    def __init__(
        self,
        radiative_temperature=None,
        dilution_factor=None,
        time_0=0 * u.d,
        time=None,
    ):
        """Instantiate RadiationField Class using IsotopeAbundance class.
           Expand BaseProperty

        Arguments
        ----------
        radiative_temperature:  astropy.units.Quantity.
                                Radiative temperature of cells.
        isotope: numpy.ndarray. Dilution factoy of cells
        time_0: astropy.units.Quantity.
        time: astropy.units.Quantity.
        """
        BaseProperty.__init__(self, time_0, time)
        if radiative_temperature is not None:
            self.radiative_temperature = self.cgs_units(radiative_temperature)
        if dilution_factor is not None:
            self.dilution_factor = dilution_factor
        self.number_of_cells = self.set_number_of_cells(
            radiative_temperature, dilution_factor
        )

    @property
    def number_of_cells(self):
        return self.__number_of_cells

    @number_of_cells.setter
    def number_of_cells(self, number_of_cells):
        self.__number_of_cells = number_of_cells

    def set_number_of_cells(self, radiative_temperature, dilution_factor):
        """Set number of cells checking if both radiative_temperature and
        dilution_factor exists and compare length of cells if both exist"""

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
        return number_of_cells

    @u.quantity_input
    def cgs_units(self, temperature: u.K):
        """Check input radiative_temperatures units and return CGS array"""
        return temperature.cgs
