from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactorLTE,
    Levels, Lines, AtomicMass, PartitionFunction,
    LevelPopulation, LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion, BetaSobolev, JBlues,
    TransitionProbabilities, StimulatedEmissionFactor, SelectedAtoms,
    PhiGeneral, PhiSahaNebular, LevelBoltzmannFactorDiluteLTE, DilutionFactor,
    ZetaData, ElectronTemperature, LinkTRadTElectron, BetaElectron,
    RadiationFieldCorrection, RadiationFieldCorrectionInput,
    LevelBoltzmannFactorNoNLTE, LevelBoltzmannFactorNLTECoronal, NLTEData,
    NLTESpecies, PreviousIterationBetaSobolevs, LTEJBlues,
    LevelBoltzmannFactorNLTEGeneral)

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData, JBlues, DilutionFactor, LinkTRadTElectron,
    RadiationFieldCorrectionInput, NLTESpecies, PreviousIterationBetaSobolevs])
basic_properties = PlasmaPropertyCollection([BetaRadiation,
    Levels, Lines, AtomicMass, LevelPopulation, PartitionFunction,
    GElectron, IonizationData, NumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, LevelNumberDensity, IonNumberDensity,
    StimulatedEmissionFactor, SelectedAtoms, PhiGeneral])
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = PlasmaPropertyCollection([BetaSobolev,
    TransitionProbabilities])
nebular_ionization_properties = PlasmaPropertyCollection([PhiSahaNebular,
    ZetaData, ElectronTemperature, BetaElectron, RadiationFieldCorrection])
dilute_lte_excitation_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorDiluteLTE])
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
nlte_coronal_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorNLTECoronal, NLTEData, NLTESpecies])
nlte_general_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorNLTEGeneral, NLTEData, NLTESpecies,
    PreviousIterationBetaSobolevs, LTEJBlues])
