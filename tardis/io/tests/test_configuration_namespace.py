from tardis.io.config_reader import TARDISConfigurationNameSpace

simple_config_dict = {'a' : {'b' : {'param1' : 1}}}

def test_simple_configuration_namespace():
    config_ns = TARDISConfigurationNameSpace(simple_config_dict)
    assert config_ns.a.b.param1 == 1



