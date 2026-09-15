"""Plasma properties for evaluator-owned Type IIP equilibrium states."""

from tardis.plasma.properties import (
    AtomicData,
    BetaElectron,
    BetaRadiation,
    ContinuumInteractionSpecies,
    DilutePlanckianRadField,
    DilutionFactor,
    ElectronDensitiesInput,
    ElectronTemperature,
    GElectron,
    HeliumTreatment,
    IonizationData,
    IonNumberDensityInput,
    JBlues,
    LevelBoltzmannFactorDiluteLTE,
    LevelBoltzmannFactorLTE,
    LevelBoltzmannFactorNoNLTE,
    LevelNumberDensityInput,
    Levels,
    Lines,
    LinesLowerLevelIndex,
    LinesUpperLevelIndex,
    LinkTRadTElectron,
    NumberDensity,
    PartitionFunction,
    PhiSahaLTE,
    PhiSahaNebular,
    PhotoIonizationData,
    RadiationFieldCorrection,
    SahaFactor,
    SelectedAtoms,
    StimulatedEmissionFactor,
    ThermalGElectron,
    ThermalLevelBoltzmannFactorLTE,
    ThermalLTEPartitionFunction,
    ThermalPhiSahaLTE,
    TimeExplosion,
    TRadiative,
    ZetaData,
)


class PlasmaPropertyCollection(list):
    """Group properties assembled for the Type IIP workflow."""


basic_inputs = PlasmaPropertyCollection(
    [
        DilutePlanckianRadField,
        NumberDensity,
        TimeExplosion,
        AtomicData,
        JBlues,
        LinkTRadTElectron,
        HeliumTreatment,
        ContinuumInteractionSpecies,
        ElectronDensitiesInput,
        IonNumberDensityInput,
        LevelNumberDensityInput,
    ]
)
basic_properties = PlasmaPropertyCollection(
    [
        TRadiative,
        DilutionFactor,
        BetaRadiation,
        Levels,
        Lines,
        GElectron,
        IonizationData,
        LinesLowerLevelIndex,
        LinesUpperLevelIndex,
        SelectedAtoms,
        ElectronTemperature,
        ThermalLevelBoltzmannFactorLTE,
        ThermalLTEPartitionFunction,
        BetaElectron,
        ThermalGElectron,
        ThermalPhiSahaLTE,
        SahaFactor,
        StimulatedEmissionFactor,
        PartitionFunction,
    ]
)
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
nebular_ionization_properties = PlasmaPropertyCollection(
    [PhiSahaNebular, ZetaData, RadiationFieldCorrection]
)
dilute_lte_excitation_properties = PlasmaPropertyCollection(
    [LevelBoltzmannFactorDiluteLTE]
)
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
nlte_properties = PlasmaPropertyCollection()
continuum_properties = PlasmaPropertyCollection([PhotoIonizationData])
macro_atom_properties = PlasmaPropertyCollection()
helium_nlte_properties = PlasmaPropertyCollection()
helium_lte_properties = PlasmaPropertyCollection()
helium_numerical_nlte_properties = PlasmaPropertyCollection()
