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
import json
import time

# For specifying error while exception handling
from socket import gaierror

import tardis

try:
    from pytest_html import __name__ as pytest_html_path
    from pytest_html.plugin import HTMLReport
    import dokuwiki
    import requests
except ImportError:
    pytest_html = None
    dokuwiki = None
    requests = None


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
            self.dokuwiki_url = ""
        else:
            self.dokuwiki_url = dokuwiki_details["url"]

    def _generate_report(self, session):
        """Writes HTML report to a temporary logfile."""
        # Little hack to include suite_time_delta in wiki overview page.
        suite_stop_time = time.time()
        self.suite_time_delta = suite_stop_time - self.suite_start_time

        report_content = super(DokuReport, self)._generate_report(session)

        # A string which holds the complete report.
        report_content = (
            "Test executed on commit "
            "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{0}]]\n\n".format(
                tardis.__githash__
            )
        ) + report_content

        # Quick hack for preventing log to be placed in narrow left out space
        report_content = report_content.replace(
            u'class="log"', u'class="log" style="clear: both"'
        )
        # It was displayed raw on wiki pages, but not needed.
        report_content = report_content.replace(u'<!DOCTYPE html>', u'')
        return report_content

    def _save_report(self, report_content):
        """
        Uploads the report and closes the temporary file. Temporary file is
        made using `tempfile` built-in module, it gets deleted upon closing.
        """
        try:
            self.doku_conn.pages.set("reports:{0}".format(
                tardis.__githash__[:7]), report_content)
        except (gaierror, TypeError):
            pass

    def _wiki_overview_entry(self):
        """Makes an entry of current test run on overview page of dokuwiki."""
        if self.errors == 0:
            if self.failed + self.xpassed == 0:
                status = "Passed"
            else:
                status = "Failed"
        else:
            status = "Errored"

        suite_start_datetime = datetime.datetime.utcfromtimestamp(self.suite_start_time)

        # Fetch commit message from github.
        gh_request = requests.get(
            "https://api.github.com/repos/tardis-sn/tardis/git/commits/{0}".format(
                tardis.__githash__
            )
        )
        gh_commit_data = json.loads(gh_request.content)
        # Pick only first line of commit message
        gh_commit_message = gh_commit_data['message'].split('\n')[0]

        # Truncate long commit messages
        if len(gh_commit_message) > 60:
            gh_commit_message = "{0}...".format(gh_commit_message[:57])
        row = "|  "
        # Append hash
        row += "[[reports:{0}|{0}]]  | ".format(tardis.__githash__[:7])
        # Append commit message
        row += "[[https://www.github.com/tardis-sn/tardis/commit/{0}|{1}]] |  ".format(
            tardis.__githash__, gh_commit_message
        )
        # Append start time
        row += "{0}  |  ".format(suite_start_datetime.strftime('%d %b %H:%M:%S'))
        # Append time elapsed
        row += "{0:.2f} sec  |  ".format(self.suite_time_delta)
        # Append status
        row += "{0}  |\n".format(status)
        try:
            self.doku_conn.pages.append('/', row)
        except (gaierror, TypeError):
            pass

    def pytest_sessionfinish(self, session):
        """
        This hook function is called by pytest when whole test run is completed.
        It calls the two helper methods `_generate_report` and `_save_report`.
        """
        report_content = self._generate_report(session)
        self._save_report(report_content)
        self._wiki_overview_entry()

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
            terminalreporter.write_sep(
                "-", "URL: {0}doku.php?id=reports:{1}".format(
                    self.dokuwiki_url, tardis.__githash__[0:7]
                )
            )
        else:
            terminalreporter.write_sep(
                "-", "Connection not established, upload failed.")
