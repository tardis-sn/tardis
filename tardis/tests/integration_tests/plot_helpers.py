import os
import tempfile

from pytest_html import extras
from tardis import __githash__ as tardis_githash


thumbnail_html_remote = """
<div class="image" style="float: left">
    <a href="#">
        <img src= "{dokuwiki_url}lib/exe/fetch.php?media=reports:{githash}:{name}.png" />
    </a>
</div>
"""


class PlotUploader(object):
    def __init__(self, request, dokuwiki_url):
        """A helper class to collect plots from integration tests and upload
        them to DokuWiki.

        Parameters
        ----------
        request : _pytest.python.RequestObject

        """
        self.request = request
        self._plots = list()
        self.plot_html = list()
        self.dokuwiki_url = dokuwiki_url

    def add(self, plot, name):
        """Accept a plot figure and add it to ``self._plots``.

        Parameters
        ----------
        plot : matplotlib.pyplot.figure
        name : str

        """
        self._plots.append((plot, name))

    def upload(self, report):
        """Upload the content of self._plots to dokuwiki.

        Parameters
        ----------
        report : _pytest.runner.TestReport

        """

        for plot, name in self._plots:
            plot_file = tempfile.NamedTemporaryFile(suffix=".png")
            axes = plot.axes[0]

            if report.passed:
                axes.text(0.8, 0.8, 'passed', transform=axes.transAxes,
                            bbox={'facecolor': 'green', 'alpha': 0.5, 'pad': 10})
            else:
                axes.text(0.8, 0.8, 'failed', transform=axes.transAxes,
                            bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})

            plot.savefig(plot_file.name)

            self.request.config.dokureport.doku_conn.medias.add(
                "reports:{0}:{1}.png".format(tardis_githash[:7], name),
                plot_file.name
            )

            self.plot_html.append(extras.html(
                thumbnail_html_remote.format(
                    dokuwiki_url=self.dokuwiki_url,
                    githash=tardis_githash[:7],
                    name=name)
                )
            )
            plot_file.close()

    def get_extras(self):
        """Return ``self.plot_html`` which is further added into html report.

        Returns
        -------
        list
             List of strings containing raw html snippet to embed images.
        """
        return self.plot_html


thumbnail_html_local = """
<div class="image" style="float: left">
    <a href="#">
        <img src= "assets/{name}.png" />
    </a>
</div>
"""


class LocalPlotSaver(object):
    def __init__(self, request, assets_dirpath):
        """A helper class to collect plots from integration tests and save
        them in a particular directory locally.

        Parameters
        ----------
        request : _pytest.python.RequestObject

        """
        self.request = request
        self._plots = list()
        self.plot_html = list()
        self.assets_dirpath = assets_dirpath

    def add(self, plot, name):
        """Accept a plot figure and add it to ``self._plots``.

        Parameters
        ----------
        plot : matplotlib.pyplot.figure
        name : str

        """
        self._plots.append((plot, name))

    def upload(self, report):
        """Upload (save) content of ``self._plots`` to ``self.assets_dirpath``.

        Parameters
        ----------
        report : _pytest.runner.TestReport

        """

        for plot, name in self._plots:
            axes = plot.axes[0]

            if report.passed:
                axes.text(0.8, 0.8, 'passed', transform=axes.transAxes,
                          bbox={'facecolor': 'green', 'alpha': 0.5, 'pad': 10})
            else:
                axes.text(0.8, 0.8, 'failed', transform=axes.transAxes,
                          bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10})

            plot.savefig(os.path.join(self.assets_dirpath, "{0}.png".format(name)))

            self.plot_html.append(extras.html(
                thumbnail_html_local.format(name=name))
            )

    def get_extras(self):
        """Return ``self.plot_html`` which is further added into html report.

        Returns
        -------
        list
             List of strings containing raw html snippet to embed images.
        """
        return self.plot_html
