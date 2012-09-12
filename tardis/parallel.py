#parallel launcher scripts
import tardis.simulation
from IPython.parallel import require

@require(tardis.simulation)
def simulation_launcher(config_dict):
    return tardis.simulation.run_multizone(config_dict, globals()['amodel'])


def parallel_macro_multizone(lbv, default_dict, **kwargs):
    #lbv = client.load_balanced_view()

    config_dicts = []

    for i in xrange(len(kwargs.items()[0])):
        tmp_config_dict = default_dict.copy()
        for key in kwargs:
            tmp_config_dict[key] = kwargs[key][i]
        config_dicts.append(tmp_config_dict)

    specs = lbv.map_async(simulation_launcher, config_dicts)
    return specs