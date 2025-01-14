from tardis.plasma.properties import *


class PlasmaPropertyCollection(list):
    pass


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
        NLTEIonizationSpecies,
        NLTEExcitationSpecies,
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
    ]
)
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = []
nebular_ionization_properties = PlasmaPropertyCollection(
    [PhiSahaNebular, ZetaData, BetaElectron, RadiationFieldCorrection]
)
dilute_lte_excitation_properties = PlasmaPropertyCollection(
    [LevelBoltzmannFactorDiluteLTE]
)
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
nlte_properties = PlasmaPropertyCollection(
    [
        LevelBoltzmannFactorNLTE,
        NLTEData,
        PreviousElectronDensities,
    ]
)
nlte_root_solver_properties = PlasmaPropertyCollection(
    [NLTEIndexHelper, NLTEPopulationSolverRoot]
)
nlte_lu_solver_properties = PlasmaPropertyCollection(
    [NLTEIndexHelper, NLTEPopulationSolverLU]
)
helium_nlte_properties = PlasmaPropertyCollection(
    [
        HeliumNLTE,
        RadiationFieldCorrection,
        ZetaData,
        BetaElectron,
        LevelNumberDensityHeNLTE,
        IonNumberDensityHeNLTE,
    ]
)
helium_lte_properties = PlasmaPropertyCollection(
    [LevelNumberDensity, IonNumberDensity]
)
helium_numerical_nlte_properties = PlasmaPropertyCollection(
    [HeliumNumericalNLTE]
)
continuum_interaction_inputs = PlasmaPropertyCollection(
    [
        PhotoIonRateCoeff,
        StimRecombRateFactor,
        BfHeatingRateCoeffEstimator,
        StimRecombCoolingRateCoeffEstimator,
        YgData,
    ]
)
continuum_interaction_properties = PlasmaPropertyCollection(
    [
        StimRecombRateCoeff,
        PhotoIonizationData,
        SpontRecombRateCoeff,
        ThermalLevelBoltzmannFactorLTE,
        ThermalLTEPartitionFunction,
        BetaElectron,
        ThermalGElectron,
        ThermalPhiSahaLTE,
        SahaFactor,
        CorrPhotoIonRateCoeff,
        SpontRecombCoolingRateCoeff,
        YgInterpolator,
        CollExcRateCoeff,
        CollDeexcRateCoeff,
        RawCollisionTransProbs,
        MarkovChainIndex,
        FreeFreeCoolingRate,
        FreeBoundCoolingRate,
        LevelNumberDensityLTE,
        PhotoIonBoltzmannFactor,
        FreeBoundEmissionCDF,
        LevelIdxs2LineIdx,
        LevelIdxs2TransitionIdx,
        CollIonRateCoeffSeaton,
        CollRecombRateCoeff,
        ContinuumInteractionHandler,
    ]
)
adiabatic_cooling_properties = PlasmaPropertyCollection([AdiabaticCoolingRate])
two_photon_properties = PlasmaPropertyCollection(
    [
        RawTwoPhotonTransProbs,
        TwoPhotonData,
        TwoPhotonEmissionCDF,
        TwoPhotonFrequencySampler,
    ]
)
