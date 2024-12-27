from typing import Union

import numpy as np
from astropy import units as u

from tardis.util.base import intensity_black_body


class DilutePlanckianRadiationField:
    """
    Represents the state of a dilute planckian radiation field.


    Parameters
    ----------
    temperature : u.Quantity
        Radiative temperature in each shell
    dilution_factor : numpy.ndarray
        Dilution Factors in each shell
    geometry: tardis.model.Radial1DModel
        The geometry of the model that uses to constrains the active shells
    """

    def __init__(
        self,
        temperature: u.Quantity,
        dilution_factor: np.ndarray,
        geometry=None,
    ):
        print("\n=== DilutePlanckianRadiationField.__init__ ===")
        print("Temperature array:")
        print(f"  First 10 elements: {temperature[:10]}")
        print(f"  Shape: {temperature.shape}")
        print(f"  dtype: {temperature.dtype}")
        print(f"  unit: {temperature.unit}")
        
        print("\nDilution factor array:")
        print(f"  First 10 elements: {dilution_factor[:10]}")
        print(f"  Shape: {dilution_factor.shape}")
        print(f"  dtype: {dilution_factor.dtype}")
        
        if geometry is not None:
            print("\nGeometry boundaries:")
            print(f"  v_inner_boundary_index: {geometry.v_inner_boundary_index}")
            print(f"  v_outer_boundary_index: {geometry.v_outer_boundary_index}")
        print("="*50)

        # ensuring that the radiation_field has both
        # dilution_factor and t_radiative equal length
        assert len(temperature) == len(dilution_factor)
        if (
            geometry is not None
        ):  # check the active shells only (this is used when setting up the radiation_field_state)
            assert np.all(
                temperature[
                    geometry.v_inner_boundary_index : geometry.v_outer_boundary_index
                ]
                > 0 * u.K
            )
            assert np.all(
                dilution_factor[
                    geometry.v_inner_boundary_index : geometry.v_outer_boundary_index
                ]
                > 0
            )
        else:
            assert np.all(temperature > 0 * u.K)
            assert np.all(dilution_factor >= 0)
        self.temperature = temperature
        self.dilution_factor = dilution_factor

    @property
    def temperature_kelvin(self):
        # print("\n=== DilutePlanckianRadiationField.temperature_kelvin ===")
        result = self.temperature.to(u.K).value
        # print(f"First 10 elements: {result[:10]}")
        # print(f"Shape: {result.shape}")
        # print(f"dtype: {result.dtype}")
        # print("="*50)
        return result

    def calculate_mean_intensity(self, nu: Union[u.Quantity, np.ndarray]):
        """
        Calculate the intensity of the radiation field at a given frequency.

        Parameters
        ----------
        nu : u.Quantity
            Frequency at which the intensity is to be calculated

        Returns
        -------
        intensity : u.Quantity
            Intensity of the radiation field at the given frequency
        """
        print("\n=== DilutePlanckianRadiationField.calculate_mean_intensity ===")
        print("Input frequency (nu):")
        print(f"  First 10 elements: {nu[:10]}")
        print(f"  Shape: {nu.shape}")
        print(f"  dtype: {nu.dtype}")
        print(f"self.temperature = {self.temperature}")
        if isinstance(nu, u.Quantity):
            print(f"  unit: {nu.unit}")
        
        result = self.dilution_factor * intensity_black_body(
            nu[np.newaxis].T, self.temperature
        )
        print("\nOutput intensity:")
        print(f"  First 10 elements: {result[:10]}")
        print(f"  Shape: {result.shape}")
        print(f"  dtype: {result.dtype}")
        if isinstance(result, u.Quantity):
            print(f"  unit: {result.unit}")
        print("="*50)
        return result

    def to_planckian_radiation_field(self):
        # print("\n=== DilutePlanckianRadiationField.to_planckian_radiation_field ===")
        # print("Converting with temperature:")
        # print(f"  First 10 elements: {self.temperature[:10]}")
        # print(f"  Shape: {self.temperature.shape}")
        # print(f"  dtype: {self.temperature.dtype}")
        # print(f"  unit: {self.temperature.unit}")
        # print("="*50)
        return PlanckianRadiationField(self.temperature)


class PlanckianRadiationField:
    """
    Represents the state of a planckian radiation field.


    Parameters
    ----------
    temperature : u.Quantity
        Radiative temperature in each shell
    dilution_factor : numpy.ndarray
        Dilution Factors in each shell
    geometry: tardis.model.Radial1DModel
        The geometry of the model that uses to constrains the active shells
    """

    def __init__(
        self,
        temperature: u.Quantity,
        geometry=None,
    ):
        if (
            geometry is not None
        ):  # check the active shells only (this is used when setting up the radiation_field_state)
            assert np.all(
                temperature[
                    geometry.v_inner_boundary_index : geometry.v_outer_boundary_index
                ]
                > 0 * u.K
            )
        else:
            assert np.all(temperature > 0 * u.K)
        self.temperature = temperature

    @property
    def temperature_kelvin(self):
        return self.temperature.to(u.K).value

    def calculate_mean_intensity(self, nu: Union[u.Quantity, np.ndarray]):
        """
        Calculate the intensity of the radiation field at a given frequency.

        Parameters
        ----------
        nu : u.Quantity
            Frequency at which the intensity is to be calculated

        Returns
        -------
        intensity : u.Quantity
            Intensity of the radiation field at the given frequency
        """
        return intensity_black_body(nu[np.newaxis].T, self.temperature)
