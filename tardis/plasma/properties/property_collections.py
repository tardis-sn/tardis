from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactor,
    Levels, Lines, SelectedAtoms, AtomicMass, LTEPartitionFunction,
    LevelPopulationLTE, LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion)

class LTEInputs(list):
    def __init__(self):
        self.append = [TRadiative, BetaRadiation, LevelBoltzmannFactor,
            Levels, Lines, AtomicData, Abundance, SelectedAtoms, AtomicMass,
            LTEPartitionFunction, LevelPopulationLTE, PhiSahaLTE, GElectron,
            IonizationData, Density, NumberDensity, IonNumberDensity,
            LevelNumberDensity, LinesLowerLevelIndex, LinesUpperLevelIndex,
            TauSobolev, TimeExplosion]