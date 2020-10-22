import os
import tempfile

from pytest_html import extras
from tardis import __githash__ as tardis_githash


class BasePlotSaver(object):
    def __init__(self, request, dokuwiki_url=None, assets_dirpath=None):
        """Base Class for RemotePlotSaver and LocalPlotSaver classes.
        These help in uploading plots to DokuWiki instance or saving them
        locally in a directory.

        Parameters
        ----------
        request : _pytest.python.RequestObject
        dokuwiki_url : str
        assets_dirpath : str
        """
        self.request = request
        self._plots = list()
        self.plot_html = list()
        self.dokuwiki_url = dokuwiki_url
        self.assets_dirpath = assets_dirpath

    def add(self, plot, name):
        """Accept a plot figure and add it to ``self._plots``.

        Parameters
        ----------
        plot : matplotlib.pyplot.figure
        name : str

        """
        self._plots.append((plot, name))

    def save(self, plot, filepath, report):
        """Mark pass / fail and save a plot with ``name`` to ``filepath``.

        Parameters
        ----------
        plot : matplotlib.pyplot.figure
        filepath : str
        report : _pytest.runner.TestReport
        """
        axes = plot.axes[0]

        if report.passed:
            axes.text(
                0.8,
                0.8,
                "passed",
                transform=axes.transAxes,
                bbox={"facecolor": "green", "alpha": 0.5, "pad": 10},
            )
        else:
            axes.text(
                0.8,
                0.8,
                "failed",
                transform=axes.transAxes,
                bbox={"facecolor": "red", "alpha": 0.5, "pad": 10},
            )

        plot.savefig(filepath)

    def get_extras(self):
        """Return ``self.plot_html`` which is further added into html report.

        Returns
        -------
        list
             List of strings containing raw html snippet to embed images.
        """
        return self.plot_html


thumbnail_html_remote = """
<div class="image" style="float: left">
    <a href="#">
        <img src= "{dokuwiki_url}/lib/exe/fetch.php?media=reports:{githash}:{name}.png" />
    </a>
</div>
"""


class RemotePlotSaver(BasePlotSaver):
    def __init__(self, request, dokuwiki_url):
        super(RemotePlotSaver, self).__init__(
            request, dokuwiki_url=dokuwiki_url
        )

    def upload(self, report):
        """Upload content of ``self._plots`` to ``self.dokuwiki_url``.

        Parameters
        ----------
        report : _pytest.runner.TestReport

        """

        for plot, name in self._plots:
            plot_file = tempfile.NamedTemporaryFile(suffix=".png")
            self.save(plot, plot_file.name, report)

            self.request.config.dokureport.doku_conn.medias.add(
                "reports:{0}:{1}.png".format(tardis_githash[:7], name),
                plot_file.name,
            )

            self.plot_html.append(
                extras.html(
                    thumbnail_html_remote.format(
                        dokuwiki_url=self.dokuwiki_url,
                        githash=tardis_githash[:7],
                        name=name,
                    )
                )
            )
            plot_file.close()


thumbnail_html_local = """
<div class="image" style="float: left">
    <a href="#">
        <img src= "assets/{name}.png" />
    </a>
</div>
"""


class LocalPlotSaver(BasePlotSaver):
    def __init__(self, request, assets_dirpath):
        super(LocalPlotSaver, self).__init__(
            request, assets_dirpath=assets_dirpath
        )

    def upload(self, report):
        """Save content of ``self._plots`` to ``self.assets_dirpath``.

        Parameters
        ----------
        report : _pytest.runner.TestReport

        """

        for plot, name in self._plots:
            self.save(
                plot,
                os.path.join(self.assets_dirpath, "{0}.png".format(name)),
                report,
            )

            self.plot_html.append(
                extras.html(thumbnail_html_local.format(name=name))
            )
