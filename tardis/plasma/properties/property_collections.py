from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactor,
    Levels, Lines, SelectedAtoms, AtomicMass, LTEPartitionFunction,
    LevelPopulationLTE, LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion, ElectronDensity)

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData])
basic_properties = PlasmaPropertyCollection([BetaRadiation,
    LevelBoltzmannFactor, Levels, Lines, SelectedAtoms, AtomicMass,
    GElectron, IonizationData, NumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, LevelNumberDensity, IonNumberDensity,
    ElectronDensity])
lte_ionization_properties = PlasmaPropertyCollection([LTEPartitionFunction,
    PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LTEPartitionFunction,
    LevelPopulationLTE])