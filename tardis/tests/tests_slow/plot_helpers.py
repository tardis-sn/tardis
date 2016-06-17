import tempfile

from pytest_html import extras
import tardis


class PlotUploader(object):
    def __init__(self, request):
        self.request = request
        self._plots = list()
        self.plot_html = list()
        self.dokuwiki_url = self.request.config.dokureport.dokuwiki_url

    def add(self, plot, name):
        """
        Accept a `matplotlib.pyplot.figure` and add it to self._plots.
        """
        self._plots.append((plot, name))

    def upload(self, report):
        """
        Upload the content in self._plots to dokuwiki.
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
                "plots:{0}_{1}.png".format(tardis.__githash__[0:7], name),
                plot_file.name
            )

            thumbnail_html = """
            <div class="image" style="float: left">
                <a href="#">
                    <img src= "{0}lib/exe/fetch.php?media=plots:{1}_{2}.png" />
                </a>
            </div>
            """.format(self.dokuwiki_url, tardis.__githash__[0:7], name)

            self.plot_html.append(extras.html(thumbnail_html))
            plot_file.close()

    def get_extras(self):
        return self.plot_html
