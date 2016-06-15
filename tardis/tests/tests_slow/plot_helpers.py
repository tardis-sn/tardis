import tempfile
from pytest_html import extras
import tardis


class PlotUploader(object):
    def __init__(self):
        self._plots = list()
        self.plot_links = list()

    def add(self, plot, name):
        """
        Accept a `matplotlib.pyplot.figure` and add it to self._plots.
        """
        self._plots.append((plot, name))

    def upload(self, request):
        """
        Upload the content in self._plots to dokuwiki.
        """
        dokuwiki_url = request.config.dokureport.dokuwiki_url
        for plot, name in self._plots:
            plot_file = tempfile.NamedTemporaryFile(suffix=".png")
            plot.savefig(plot_file.name)

            request.config.dokureport.doku_conn.medias.add(
                "plots:{0}_{1}.png".format(tardis.__githash__[0:7], name),
                plot_file.name
            )
            self.plot_links.append(extras.url(
                "{0}lib/exe/fetch.php?media=plots:{1}_{2}.png".format(
                    dokuwiki_url, tardis.__githash__[0:7], name
                ), name)
            )
            plot_file.close()

    def get_extras(self):
        return self.plot_links
