from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactorLTE,
    Levels, Lines, AtomicMass, PartitionFunction,
    LevelPopulation, LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion, ElectronDensity, BetaSobolev, JBlues,
    TransitionProbabilities, StimulatedEmissionFactor, SelectedAtoms,
    PhiGeneral, PhiSahaNebular, LevelBoltzmannFactorDiluteLTE, DilutionFactor,
    ZetaData, ElectronTemperature, LinkTRadTElectron, BetaElectron,
    RadiationFieldCorrection, RadiationFieldCorrectionInput, NLTEInput,
    LevelBoltzmannFactorNLTE, NLTEData, LevelBoltzmannFactorWithNLTE,
    LevelBoltzmannFactorWithoutNLTE)

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData, JBlues, DilutionFactor, LinkTRadTElectron,
    RadiationFieldCorrectionInput])
basic_properties = PlasmaPropertyCollection([BetaRadiation,
    Levels, Lines, AtomicMass, LevelPopulation, PartitionFunction,
    GElectron, IonizationData, NumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, LevelNumberDensity, IonNumberDensity,
    ElectronDensity, StimulatedEmissionFactor, SelectedAtoms, PhiGeneral,
    ElectronTemperature])
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = PlasmaPropertyCollection([BetaSobolev,
    TransitionProbabilities])
nebular_ionization_properties = PlasmaPropertyCollection([PhiSahaNebular,
    ZetaData, BetaElectron, RadiationFieldCorrection])
dilute_lte_excitation_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorDiluteLTE])
nlte_properties = PlasmaPropertyCollection([NLTEInput,
    LevelBoltzmannFactorNLTE, NLTEData, LevelBoltzmannFactorWithNLTE])
non_nlte_properties = \
    PlasmaPropertyCollection([LevelBoltzmannFactorWithoutNLTE])