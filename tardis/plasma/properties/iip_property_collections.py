"""Property collections used by the equilibrium-backed IIP plasma."""

from tardis.opacities.tau_sobolev import BetaSobolev, TauSobolev
from tardis.plasma.properties import (
    AtomicData,
    BetaElectron,
    BetaRadiation,
    BoundFreeHeatingEstimator,
    CollisionalDeexcitationRateCoefficient,
    CollisionalExcitationRateCoefficient,
    CollisionalIonizationRateCoefficient,
    ContinuumInteractionSpecies,
    DilutePlanckianRadField,
    DilutionFactor,
    ElectronTemperature,
    FreeFreeHeatingEstimator,
    GElectron,
    HydrogenContinuumFractionalHeating,
    HydrogenContinuumIonRatio,
    HydrogenContinuumLTEPopulations,
    IIPContinuumPopulations,
    IonizationData,
    Iteration,
    JBlues,
    LevelBoltzmannFactorDiluteLTE,
    LevelBoltzmannFactorNoNLTE,
    Levels,
    Lines,
    LinesLowerLevelIndex,
    LinesUpperLevelIndex,
    LinkTRadTElectron,
    MacroAtomData,
    NLTEData,
    NumberDensity,
    PartitionFunction,
    PhiSahaNebular,
    PhotoIonizationData,
    PhotoionizationRateEstimator,
    PreviousBetaSobolev,
    PreviousElectronDensities,
    PreviousIonNumberDensity,
    PreviousLevelNumberDensity,
    RadiationFieldCorrection,
    SelectedAtoms,
    StimulatedEmissionFactor,
    StimulatedRecombinationCoolingEstimator,
    StimulatedRecombinationRateEstimator,
    ThermalGElectron,
    ThermalLevelBoltzmannFactorLTE,
    ThermalLTEPartitionFunction,
    TimeExplosion,
    TRadiative,
    ZetaData,
)


class PlasmaPropertyCollection(list):
    """Group plasma properties used by the IIP configuration."""


basic_inputs = PlasmaPropertyCollection(
    [
        DilutePlanckianRadField,
        NumberDensity,
        TimeExplosion,
        AtomicData,
        JBlues,
        LinkTRadTElectron,
        ContinuumInteractionSpecies,
    ]
)
basic_properties = PlasmaPropertyCollection(
    [
        TRadiative,
        DilutionFactor,
        BetaRadiation,
        Levels,
        Lines,
        PartitionFunction,
        GElectron,
        IonizationData,
        LinesLowerLevelIndex,
        LinesUpperLevelIndex,
        StimulatedEmissionFactor,
        SelectedAtoms,
        ElectronTemperature,
        ThermalLevelBoltzmannFactorLTE,
        ThermalLTEPartitionFunction,
        BetaElectron,
        ThermalGElectron,
    ]
)
dilute_lte_excitation_properties = PlasmaPropertyCollection(
    [LevelBoltzmannFactorDiluteLTE]
)
nebular_ionization_properties = PlasmaPropertyCollection(
    [PhiSahaNebular, ZetaData, RadiationFieldCorrection]
)
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
macro_atom_properties = PlasmaPropertyCollection(
    [
        TauSobolev,
        BetaSobolev,
        MacroAtomData,
    ]
)

hydrogen_continuum_inputs = PlasmaPropertyCollection(
    [
        Iteration,
        PhotoIonizationData,
        PreviousElectronDensities,
        PreviousBetaSobolev,
        PreviousIonNumberDensity,
        PreviousLevelNumberDensity,
        CollisionalIonizationRateCoefficient,
        CollisionalExcitationRateCoefficient,
        CollisionalDeexcitationRateCoefficient,
        FreeFreeHeatingEstimator,
        BoundFreeHeatingEstimator,
        StimulatedRecombinationCoolingEstimator,
        PhotoionizationRateEstimator,
        StimulatedRecombinationRateEstimator,
    ]
)
hydrogen_continuum_properties = PlasmaPropertyCollection(
    [
        NLTEData,
        IIPContinuumPopulations,
        HydrogenContinuumLTEPopulations,
        HydrogenContinuumFractionalHeating,
        HydrogenContinuumIonRatio,
    ]
)
