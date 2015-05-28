from tardis.plasma.properties import (BetaRadiation, LevelBoltzmannFactor,
    Levels, Lines, SelectedAtoms, AtomicMass, LTEPartitionFunction,
    LevelPopulationLTE, LevelNumberDensity, PhiSahaLTE, GElectron,
    IonizationData, NumberDensity, IonNumberDensity, LinesLowerLevelIndex,
    LinesUpperLevelIndex, TauSobolev, TRadiative, AtomicData, Abundance,
    Density, TimeExplosion)

class LTEProperties(list):
    def __init__(self):
        properties = (BetaRadiation, LevelBoltzmannFactor, Levels, Lines, 
            SelectedAtoms, AtomicMass, LTEPartitionFunction, 
            LevelPopulationLTE, PhiSahaLTE, GElectron, IonizationData, 
            NumberDensity, IonNumberDensity, LevelNumberDensity, 
            LinesLowerLevelIndex, LinesUpperLevelIndex, TauSobolev,
            TRadiative, Abundance, Density, TimeExplosion, AtomicData)
        for module in properties:
            self.append(module)
