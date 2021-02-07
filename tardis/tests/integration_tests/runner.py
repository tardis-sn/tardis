import argparse
import datetime
import json
import logging
import subprocess
import time
import yaml

import dokuwiki
import requests


logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Run slow integration tests")
parser.add_argument(
    "--integration-tests",
    dest="yaml_filepath",
    help="Path to YAML config file for integration tests.",
)
parser.add_argument(
    "--tardis-refdata",
    dest="tardis_refdata",
    help="Path to Tardis Reference Data.",
)
parser.add_argument(
    "--less-packets",
    action="store_true",
    default=False,
    help="Run integration tests with less packets.",
)


def run_tests():
    args = parser.parse_args()

    integration_tests_config = yaml.load(
        open(args.yaml_filepath), Loader=yaml.CLoader
    )
    doku_conn = dokuwiki.DokuWiki(
        url=integration_tests_config["dokuwiki"]["url"],
        user=integration_tests_config["dokuwiki"]["username"],
        password=integration_tests_config["dokuwiki"]["password"],
    )
    less_packets = "--less-packets" if args.less_packets else ""
    test_command = [
        "python",
        "setup.py",
        "test",
        "--test-path=tardis/tests/integration_tests/test_integration.py",
        "--args",
<<<<<<< HEAD
        "--capture=no --integration-tests={0} --tardis-refdata={1} --remote-data "
        "{2}".format(args.yaml_filepath, args.tardis_refdata, less_packets),
=======
        f"--capture=no --integration-tests={args.yaml_filepath} --tardis-refdata={args.tardis_refdata} --remote-data "
        f"{less_packets}",
>>>>>>> 5c7f60f3... all string formatting done
    ]
    subprocess.call(test_command)

    while True:
        # Request Github API and get githash of master on Github.
        gh_request = requests.get(
            "https://api.github.com/repos/tardis-sn/tardis/branches/master"
        )
        gh_master_head_data = json.loads(gh_request.content)
        gh_tardis_githash = gh_master_head_data["commit"]["sha"][:7]

        # Check whether a report of this githash is uploaded on dokuwiki.
        # If not, then this is a new commit and tests should be executed.
        dokuwiki_report = doku_conn.pages.get(
<<<<<<< HEAD
            "reports:{0}".format(gh_tardis_githash)
=======
            f"reports:{gh_tardis_githash}"
>>>>>>> 5c7f60f3... all string formatting done
        )

        # If dokuwiki returns empty string, then it means that report has not
        # been created yet.
        if len(dokuwiki_report) == 0:
            subprocess.call(
                [
                    "git",
                    "pull",
                    "https://www.github.com/tardis-sn/tardis",
                    "master",
                ]
            )
            subprocess.call(test_command)
        else:
            checked = datetime.datetime.now()
            logger.info(
<<<<<<< HEAD
                "Up-to-date. Checked on {0} {1}".format(
                    checked.strftime("%d-%b-%Y"), checked.strftime("%H:%M:%S")
                )
=======
                f'Up-to-date. Checked on {checked.strftime("%d-%b-%Y")} {checked.strftime("%H:%M:%S")}'
>>>>>>> 5c7f60f3... all string formatting done
            )
            time.sleep(600)
