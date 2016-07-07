import argparse
import datetime
import json
import logging
import subprocess
import time

import requests
from tardis import __githash__ as tardis_githash


logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Run slow integration tests")
parser.add_argument("--integration-tests", dest="yaml_filepath",
                    help="Path to YAML config file for integration tests.")
parser.add_argument("--atomic-dataset", dest="atomic_dataset",
                    help="Path to atomic dataset.")


if __name__ == "__main__":
    args = parser.parse_args()
    while True:
        gh_request = requests.get(
            "https://api.github.com/repos/tardis-sn/tardis/branches/master"
        )
        gh_master_head_data = json.loads(gh_request.content)
        gh_tardis_githash = gh_master_head_data['commit']['sha']

        if gh_tardis_githash != tardis_githash:
            subprocess.call([
                "git", "pull", "https://www.github.com/tardis-sn/tardis", "master"
            ])

            subprocess.call([
                "python", "setup.py", "test",
                "--test-path=tardis/tests/tests_slow/test_integration.py", "--args",
                "-rs --integration-tests={0} --atomic-dataset={1} --remote-data".format(
                    args.yaml_filepath, args.atomic_dataset
                )
            ])
        else:
            checked = datetime.datetime.now()
            logger.info("Up-to-date. Checked on {0} {1}".format(
                checked.strftime("%d-%b-%Y"), checked.strftime("%H:%M:%S")
            ))
            time.sleep(600)
