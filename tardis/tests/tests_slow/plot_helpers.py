import tempfile

from pytest_html import extras
import tardis


thumbnail_html = """
<div class="image" style="float: left">
    <a href="#">
        <img src= "{dokuwiki_url}lib/exe/fetch.php?media=reports:{githash}:{name}.png" />
    </a>
</div>
"""


class PlotUploader(object):
    def __init__(self, request):
        """A helper class to collect plots from integration tests and upload
        them to DokuWiki.

        Parameters
        ----------
        request : _pytest.python.RequestObject

        """
        self.request = request
        self._plots = list()
        self.plot_html = list()
        self.dokuwiki_url = self.request.config.dokureport.dokuwiki_url

    def add(self, plot, name):
        """Accept a plot figure and add it to ``self._plots``.

        Parameters
        ----------
        plot : matplotlib.pyplot.figure
        name : str

        """
        self._plots.append((plot, name))

    def upload(self, report):
        """Upload the content in self._plots to dokuwiki.

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
                "reports:{0}:{1}.png".format(tardis.__githash__[0:7], name),
                plot_file.name
            )

            self.plot_html.append(extras.html(
                thumbnail_html.format(
                    dokuwiki_url=self.dokuwiki_url,
                    githash=tardis.__githash__[0:7],
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
