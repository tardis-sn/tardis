from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactor,
    Levels, Lines, SelectedAtoms, AtomicMass, LTEPartitionFunction,
    LevelPopulationLTE, LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion)

class PlasmaPropertyCollection(list):
    pass

basic_inputs = PlasmaPropertyCollection([TRadiative, Abundance, Density,
    TimeExplosion, AtomicData])
lte_processing_properties = PlasmaPropertyCollection([BetaRadiation,
    LevelBoltzmannFactor, Levels, Lines, SelectedAtoms, AtomicMass,
    LTEPartitionFunction, LevelPopulationLTE, LevelNumberDensity, PhiSahaLTE,
    GElectron, IonizationData, NumberDensity, IonNumberDensity,
    LinesLowerLevelIndex, LinesUpperLevelIndex, TauSobolev])
