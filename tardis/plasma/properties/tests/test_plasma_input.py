import pandas as pd

from tardis.plasma.base import BasePlasma
from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.properties.plasma_input import (
    ElectronDensitiesInput,
    IonNumberDensityInput,
    LevelNumberDensityInput,
)


class EquilibriumStateObserver(ProcessingPlasmaProperty):
    outputs = ("equilibrium_state_observation",)

    def __init__(self, plasma_parent: BasePlasma) -> None:
        self.calculation_count = 0
        super().__init__(plasma_parent)

    def calculate(
        self,
        electron_densities: pd.Series,
        ion_number_density: pd.DataFrame,
        level_number_density: pd.DataFrame,
    ) -> tuple[
        pd.Series,
        pd.DataFrame,
        pd.DataFrame,
    ]:
        self.calculation_count += 1
        return (
            electron_densities,
            ion_number_density,
            level_number_density,
        )


def test_equilibrium_inputs_are_published_before_downstream_update() -> None:
    input_properties = (
        ElectronDensitiesInput,
        IonNumberDensityInput,
        LevelNumberDensityInput,
    )
    initial_values = {
        "electron_densities": pd.Series([1.0]),
        "ion_number_density": pd.DataFrame([[2.0]]),
        "level_number_density": pd.DataFrame([[3.0]]),
    }
    plasma = BasePlasma(
        [*input_properties, EquilibriumStateObserver], **initial_values
    )
    observer = plasma.plasma_properties_dict["EquilibriumStateObserver"]
    published_values = {
        "electron_densities": pd.Series([10.0]),
        "ion_number_density": pd.DataFrame([[20.0]]),
        "level_number_density": pd.DataFrame([[30.0]]),
    }

    plasma.update(**published_values)

    assert observer.calculation_count == 2
    assert all(
        actual is expected
        for actual, expected in zip(
            plasma.equilibrium_state_observation,
            published_values.values(),
            strict=True,
        )
    )
    assert {
        type(plasma.outputs_dict[output_name])
        for output_name in published_values
    } == set(input_properties)
