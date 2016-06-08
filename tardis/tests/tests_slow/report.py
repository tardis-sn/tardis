"""
A helper class which works as a plugin to generate the test report and upload it
to the group server's dokuwiki. It inheirts from the class `HTMLReport` of
the `pytest-html` plugin. The test report contains the following details:

* The git commit hash on which test run was executed.
* The time of generation of test report.
* Number of passes, fails, errors, skips etc.
* Tabular representation of each method - name, result, duration.
* Embedded image of plot(s) and error log below a particular method (if any).

As a subclass, this class serves as a plugin and hence, `pytest-html` has to be
unregistered during the test run for tis plugin to function.

When the integration tests are selected for a particular test run, this class
is registered as a plugin in `pytest_configure` and subsequently unregistered in
`pytest_unconfigure`. As a plugin, it implements several "hook" functions
specified in pytest's official documentation.


References
==========
1. "Writing Plugins" ( https://pytest.org/latest/writing_plugins.html )
2. "Hookspec Source" ( https://pytest.org/latest/_modules/_pytest/hookspec.html )
3. "pytest-html" ( https://www.github.com/davehunt/pytest-html )
"""
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

    def __init__(self, dokuwiki_details):
        """
        Initialization of a DokuReport object and registration as a plugin
        occurs in `pytest_configure`, where a dict containing url, username and
        password of dokuwiki is passed through `dokuwiki_details`.
        """
        # Base class accepts a file path to save the report, but we pass an
        # empty string and then delete it anyhow.
        super(DokuReport, self).__init__(" ")
        del self.logfile

        try:
            self.doku_conn = dokuwiki.DokuWiki(
                url=dokuwiki_details["url"],
                user=dokuwiki_details["username"],
                password=dokuwiki_details["password"])
        except (TypeError, gaierror, dokuwiki.DokuWikiError):
            self.doku_conn = None

    def _generate_report(self, session):
        """
        The method writes HTML report to a temporary logfile.
        """
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

        if session.config._environment:
            environment = set(session.config._environment)
            body.append(html.h2('Environment'))
            body.append(html.table(
                [html.tr(html.td(e[0]), html.td(e[1])) for e in sorted(
                    environment, key=lambda e: e[0]) if e[1]],
                id='environment'))

        body.extend(summary)
        body.extend(results)

        doc = html.html(head, body)

        # A string which holds the complete report.
        report_content = (
            "Test executed on commit "
            "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{0}]]\n\n".format(
                tardis.__githash__
            )
        )
        report_content += doc.unicode(indent=2)
        return report_content

    def _save_report(self, report_content):
        """
        The method uploads the report and closes the temporary file. Temporary
        file is made using `tempfile` built-in module, it gets deleted upon
        closing.
        """
        try:
            self.doku_conn.pages.set("reports:{0}".format(
                tardis.__githash__[:7]), report_content)
        except (gaierror, TypeError):
            pass

    def pytest_sessionfinish(self, session):
        """
        This hook function is called by pytest when whole test run is completed.
        It calls the two helper methods `_generate_report` and `_save_report`.
        """
        report_content = self._generate_report(session)
        self._save_report(report_content)

    def pytest_terminal_summary(self, terminalreporter):
        """
        This hook is called by pytest after session ends, and it adds an extra
        summary at the end. Here, the success / failure of upload of report
        to dokuwiki is logged.
        """
        try:
            uploaded_report = self.doku_conn.pages.get(
                "reports:{0}".format(tardis.__githash__[0:7]))
        except (gaierror, TypeError):
            uploaded_report = ""

        if len(uploaded_report) > 0:
            terminalreporter.write_sep(
                "-", "Successfully uploaded report to Dokuwiki")
        else:
            terminalreporter.write_sep(
                "-", "Connection not established, upload failed.")
