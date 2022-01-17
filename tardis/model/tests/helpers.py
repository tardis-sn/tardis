import pytest
import os
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.io.decay import IsotopeAbundances

def data_path(filename):
    return os.path.abspath(os.path.join("tardis/io/tests/data/", filename))


class TestSetup:
    def edit_config(self, config, options):
        for item in options:
            if isinstance(options[item], dict):
                self.edit_config(config[item], options[item])
            else:
                config[item] = options[item]
        return config

    @pytest.fixture(autouse=True, scope="class")
    def auto_inject(self, request):
        # class attributes need to be defined explicitly
        cls = type(self)
        filename = self.scenario["filename"]

        config = Configuration.from_yaml(data_path(filename))
        if "config_options" in self.scenario:
            config = self.edit_config(
                config, options=self.scenario["config_options"]
            )
        model = Radial1DModel.from_config(config)

        setattr(cls, "config", config)
        setattr(cls, "model", model)
