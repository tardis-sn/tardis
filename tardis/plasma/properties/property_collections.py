from tardis.plasma.properties import *

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData, JBlues, DilutionFactor, LinkTRadTElectron,
    HeliumTreatment])
basic_properties = PlasmaPropertyCollection([BetaRadiation,
    Levels, Lines, AtomicMass, PartitionFunction,
    GElectron, IonizationData, NumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, IonNumberDensity,
    StimulatedEmissionFactor, SelectedAtoms, ElectronTemperature,
    RadiationFieldCorrection])
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = PlasmaPropertyCollection([BetaSobolev,
    TransitionProbabilities])
nebular_ionization_properties = PlasmaPropertyCollection([PhiSahaNebular,
    ZetaData, BetaElectron])
dilute_lte_excitation_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorDiluteLTE])
non_nlte_properties = PlasmaPropertyCollection([LevelBoltzmannFactorNoNLTE])
nlte_properties = PlasmaPropertyCollection([
    LevelBoltzmannFactorNLTE, NLTEData, LTEJBlues, PreviousElectronDensities,
    PreviousBetaSobolev, BetaSobolev])
helium_nlte_properties = PlasmaPropertyCollection([HeliumNLTE,
    ZetaData, BetaElectron, LevelNumberDensityHeNLTE])
helium_lte_properties = PlasmaPropertyCollection([LevelNumberDensity])
helium_numerical_nlte_properties = PlasmaPropertyCollection([
    HeliumNumericalNLTE])