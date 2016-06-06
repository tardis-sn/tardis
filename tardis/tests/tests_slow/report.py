import datetime
import pkg_resources
import os
import time

# For specifying error while exception handling
from socket import gaierror

from py.xml import html, raw
from pytest_html.plugin import HTMLReport
import tardis

try:
    import dokuwiki
except ImportError:
    dokuwiki = None


class DokuReport(HTMLReport):

    def __init__(self, logfile, dokuwiki_details):
        super(DokuReport, self).__init__(logfile.name)
        self.logfile = logfile

        if dokuwiki is not None:
            try:
                self.doku_conn = dokuwiki.DokuWiki(
                    url=dokuwiki_details["url"],
                    user=dokuwiki_details["username"],
                    password=dokuwiki_details["password"])
            except gaierror, dokuwiki.DokuWikiError:
                self.doku_conn = None
                print "Dokuwiki connection could not be established!"
        else:
            self.doku_conn = None

    def _generate_report(self):
        suite_stop_time = time.time()
        suite_time_delta = suite_stop_time - self.suite_start_time
        numtests = self.passed + self.failed + self.xpassed + self.xfailed
        generated = datetime.datetime.now()

        style_css = pkg_resources.resource_string(
            __name__, os.path.join('resources', 'style.css'))

        head = html.head(
            html.meta(charset='utf-8'),
            html.title('Test Report'),
            html.style(raw(style_css)))

        summary = [html.h2('Summary'), html.p(
            '{0} tests ran in {1:.2f} seconds.'.format(
                numtests, suite_time_delta),
            html.br(),
            html.span('{0} passed'.format(
                self.passed), class_='passed'), ', ',
            html.span('{0} skipped'.format(
                self.skipped), class_='skipped'), ', ',
            html.span('{0} failed'.format(
                self.failed), class_='failed'), ', ',
            html.span('{0} errors'.format(
                self.errors), class_='error'), '.',
            html.br(),
            html.span('{0} expected failures'.format(
                self.xfailed), class_='skipped'), ', ',
            html.span('{0} unexpected passes'.format(
                self.xpassed), class_='failed'), '.')]

        results = [html.h2('Results'), html.table([html.thead(
            html.tr([
                html.th('Result',
                        class_='sortable initial-sort result',
                        col='result'),
                html.th('Test', class_='sortable', col='name'),
                html.th('Duration',
                        class_='sortable numeric',
                        col='duration'),
                html.th('Links')]), id='results-table-head'),
            html.tbody(*self.test_logs, id='results-table-body')],
            id='results-table')]

        main_js = pkg_resources.resource_string(
            __name__, os.path.join('resources', 'main.js'))

        body = html.body(
            html.script(raw(main_js)),
            html.p('Report generated on {0} at {1}'.format(
                generated.strftime('%d-%b-%Y'),
                generated.strftime('%H:%M:%S'))))

        body.extend(summary)
        body.extend(results)

        doc = html.html(head, body)

        self.logfile.write(
            "Test executed on commit "
            "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{0}]]\n\n".format(
                tardis.__githash__
            )
        )
        unicode_doc = doc.unicode(indent=2)

        self.logfile.write(unicode_doc)

    def _save_report(self):
        # Go back to the beginning to read everything
        self.logfile.seek(0)

        if self.doku_conn is not None:
            self.doku_conn.pages.set("reports:{0}".format(
                    tardis.__githash__[:7]), self.logfile.read())
            print "Uploaded report on Dokuwiki."

        self.logfile.close()

    def pytest_sessionfinish(self, session):
        self._generate_report()
        self._save_report()
