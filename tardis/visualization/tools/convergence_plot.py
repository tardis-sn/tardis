from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors as clr
import plotly.graph_objects as go
import matplotlib as mpl


def transistion_colors(name="jet", iterations=20):
    cmap = mpl.cm.get_cmap(name, iterations)
    colors = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3]
        colors.append(mpl.colors.rgb2hex(rgb))


class ConvergencePlots:
    def __init__(self):
        self.iterable_data = {}
        self.value_data = defaultdict(list)
        self.rows = 4
        self.cols = 2
        self.specs = [
            [{}, {}],
            [{"colspan": 2}, None],
            [{"colspan": 2}, None],
            [{"colspan": 2}, None],
        ]
        self.row_heights = [0.45, 0.1, 0.4, 0.1]
        self.vertical_spacing = 0.07

    def fetch_data(self, name=None, value=None):
        """
        This allows user to fetch data from the Simulation class.
        This data is stored and used when an iteration is completed.
        Returns:
            iterable_data: iterable data from the simulation, values like radiation temperature
            value_data: values like luminosity
        """
        self.iterable_data[name] = value
        self.value_data[name].append(value)

    def get_data(self):
        """
        Returns data for convergence plots

        Returns:
            iterable_data: iterable data from the simulation, values like radiation temperature
            value_data: values like luminosity
        """
        return self.iterable_data, self.value_data


class BuildCplots(ConvergencePlots):
    def __init__(self):
        super().__init__()
        self.use_vbox = True


class UpdateCplots(BuildCplots):
    def __init__(self):
        super().__init__()
