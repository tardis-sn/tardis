

class PlasmaProperty(object):
    def __init__(self, plasma_parent):
        self.plasma_parent = plasma_parent

    def get_dependencies(self):
        return self.dependencies

    def get_fullfilled_dependencies(self):
        return X

    def get_required_dependencies(self):
        return Y

    def register_dependency(self, name):
        pass

    def f1(self):
        self.plasma_parent[]

class PartitionFunction(PlasmaProperty):
    """
    Calculate partition functions for the ions using the following formula, where
    :math:`i` is the atomic_number, :math:`j` is the ion_number and :math:`k` is the level number.

    .. math::
        Z_{i,j} = \\sum_{k=0}^{max(k)_{i,j}} g_k \\times e^{-E_k / (k_\\textrm{b} T)}



    if self.initialize is True set the first time the partition functions are initialized.
    This will set a self.partition_functions and initialize with LTE conditions.


    Returns
    -------

    partition_functions : `~astropy.table.Table`
        with fields atomic_number, ion_number, partition_function

    """

    dependency = ['number_density', 'levels']


    def evaluate(self, *args, **kwargs):
        levels = self.atom_data.levels

        level_population_proportional_array = levels.g.values[np.newaxis].T *\
                                              np.exp(np.outer(levels.energy.values, -self.beta_rads))
        level_population_proportionalities = pd.DataFrame(level_population_proportional_array,
                                                               index=self.atom_data.levels.index,
                                                               columns=np.arange(len(self.t_rads)), dtype=np.float64)


        #level_props = self.level_population_proportionalities

        partition_functions = level_population_proportionalities[self.atom_data.levels.metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()
        partition_functions_non_meta = self.ws * level_population_proportionalities[~self.atom_data.levels.metastable].groupby(
            level=['atomic_number', 'ion_number']).sum()
        partition_functions.ix[partition_functions_non_meta.index] += partition_functions_non_meta
        if self.nlte_config is not None and self.nlte_config.species != [] and not initialize_nlte:
            for species in self.nlte_config.species:
                partition_functions.ix[species] = self.atom_data.levels.g.ix[species].ix[0] * \
                                                       (self.level_populations.ix[species] /
                                                        self.level_populations.ix[species].ix[0]).sum()

        return level_population_proportionalities, partition_functions


pf = PartitionFunction()