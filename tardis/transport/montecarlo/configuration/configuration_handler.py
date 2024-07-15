from tardis.transport.montecarlo.configuration import montecarlo_globals


def update_mc_global(global_str: str, value):

    current_value = getattr(montecarlo_globals, global_str)
    if value != current_value:
        montecarlo_globals.RECOMPILE_FLAG = True
        setattr(montecarlo_globals, global_str, value)

def dynamic_update_jit(jitted_function):

    
    if montecarlo_globals.RECOMPILE_FLAG:
        jitted_function.recompile()
    return jitted_function