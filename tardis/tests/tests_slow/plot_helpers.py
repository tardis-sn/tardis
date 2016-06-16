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

    def upload(self):
        """
        Upload the content in self._plots to dokuwiki.
        """
        for plot, name in self._plots:
            plot_file = tempfile.NamedTemporaryFile(suffix=".png")
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
