from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactorLTE,
    Levels, Lines, AtomicMass, PartitionFunction,
    LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion, BetaSobolev, JBlues,
    TransitionProbabilities, StimulatedEmissionFactor, SelectedAtoms,
    PhiSahaNebular, LevelBoltzmannFactorDiluteLTE, DilutionFactor,
    ZetaData, ElectronTemperature, LinkTRadTElectron, BetaElectron,
    RadiationFieldCorrection, RadiationFieldCorrectionInput,
    LevelBoltzmannFactorNoNLTE, LevelBoltzmannFactorNLTE, NLTEExcitationData,
    NLTEExcitationSpecies, PreviousBetaSobolevs, LTEJBlues,
    PreviousElectronDensities, Chi0, HeliumNLTE, NLTEIonizationData,
    PhiSahaNLTE, PhiSahaNoNLTE, NLTEIonizationSpecies)

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData, JBlues, DilutionFactor, LinkTRadTElectron,
    RadiationFieldCorrectionInput, NLTEExcitationSpecies,
    NLTEIonizationSpecies, PreviousBetaSobolevs, PreviousElectronDensities])
basic_properties = PlasmaPropertyCollection([BetaRadiation,
    Levels, Lines, AtomicMass, PartitionFunction,
    GElectron, IonizationData, NumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, LevelNumberDensity, IonNumberDensity,
    StimulatedEmissionFactor, SelectedAtoms, ElectronTemperature])
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = PlasmaPropertyCollection([BetaSobolev,
    TransitionProbabilities])
nebular_ionization_properties = PlasmaPropertyCollection([PhiSahaNebular,
    ZetaData, BetaElectron, RadiationFieldCorrection, Chi0])
dilute_lte_excitation_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorDiluteLTE])
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
nlte_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorNLTE, NLTEExcitationData, NLTEExcitationSpecies,
    LTEJBlues])
helium_nlte_properties = PlasmaPropertyCollection([HeliumNLTE,
    RadiationFieldCorrection, ZetaData,
    BetaElectron, Chi0])
nlte_ionization_properties = PlasmaPropertyCollection([NLTEIonizationData,
    PhiSahaNLTE])
non_nlte_ionization_properties = PlasmaPropertyCollection([PhiSahaNoNLTE])