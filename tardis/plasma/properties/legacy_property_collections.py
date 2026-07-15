from tardis.opacities.macro_atom.base import (
    NonMarkovChainTransitionProbabilities,
    TransitionProbabilities,
)
from tardis.opacities.tau_sobolev import BetaSobolev, TauSobolev
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
        TauSobolev,
        StimulatedEmissionFactor,
        SelectedAtoms,
        ElectronTemperature,
    ]
)
lte_ionization_properties = PlasmaPropertyCollection([PhiSahaLTE])
lte_excitation_properties = PlasmaPropertyCollection([LevelBoltzmannFactorLTE])
macro_atom_properties = PlasmaPropertyCollection(
    [BetaSobolev, TransitionProbabilities, MacroAtomData]
)
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
        PreviousBetaSobolev,
        BetaSobolev,
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
