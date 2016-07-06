import argparse
import json
import os
import time

import requests
from tardis import __githash__ as tardis_githash


parser = argparse.ArgumentParser(description="Run slow integration tests")
parser.add_argument("--yaml", dest="yaml_filepath",
                    help="Path to YAML config file for integration tests.")
parser.add_argument("--atomic-dataset", dest="atomic_dataset",
                    help="Path to atomic dataset.")

test_command = (
    "python setup.py test --test-path=tardis/tests/tests_slow/test_integration.py "
    "--args=\"-rs --integration-tests={0} --atomic-dataset={1} --remote-data\""
)


if __name__ == "__main__":
    args = parser.parse_args()
    while True:
        gh_request = requests.get(
            "https://api.github.com/repos/tardis-sn/tardis/branches/master"
        )
        gh_master_head_data = json.loads(gh_request.content)
        gh_tardis_githash = gh_master_head_data['commit']['sha']

        if gh_tardis_githash != tardis_githash:
            os.system("git pull origin master")
            os.system(test_command.format(args.yaml_filepath,
                                          args.atomic_dataset))
        else:
            time.sleep(600)
