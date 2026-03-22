# Interactive Widgets Demo

An Interactive demonstration of Shell Info Widget and Line Info Widget.


First create and run a simulation that we can use to generate widgets (more details about running simulation in [Quickstart](https://tardis-sn.github.io/tardis/quickstart/quickstart.html) section):

```{pyodide}
:skip-embed:
sim = None
```

```python
from tardis import run_tardis
from tardis.io.atom_data import download_atom_data
from tardis.io.configuration.config_reader import Configuration

# We download the atomic data needed to run the simulation
download_atom_data('kurucz_cd23_chianti_H_He_latest')

config = Configuration.from_yaml("tardis_example.yml")
config.montecarlo.tracking.track_rpacket = True

sim = run_tardis(config, virtual_packet_logging=True)
```

Now, import functions \& class to create widgets from `visualization` subpackage:

```{pyodide}
:skip-embed:
from tardis.visualization import (
    shell_info_from_simulation,
    LineInfoWidget,
    )
```

```{note}
Click run icon in any cell to view interactive widgets!
```


## Shell Info Widget

This widget allows you to explore chemical abundances of each shell \- all the way from elements to ions to levels \- by just clicking on the rows you want to explore!

### Using a Simulation object

We will use the simulation object we created in the beginning, `sim` to generate shell info widget. Then simply display it to start using.

```{pyodide}
:skip-embed:
shell_info_widget = shell_info_from_simulation(sim)
shell_info_widget.display()
```

## Line Info Widget

This widget lets you explore the atomic lines responsible for producing features in the simulated spectrum.

You can select any wavelength range in the spectrum interactively to display a table giving the fraction of packets that experienced their last interaction with each species. Using toggle buttons, you can specify whether to filter the selected range by the emitted or absorbed wavelengths of packets. Clicking on a row in the species table, shows packet counts for each last line interaction of the selected species, which can be grouped in several ways.


To generate line info widget, we will again use the simulation object `sim` and then display the widget:

```{pyodide}
:skip-embed:
line_info_widget = LineInfoWidget.from_simulation(sim)
line_info_widget.display()
```