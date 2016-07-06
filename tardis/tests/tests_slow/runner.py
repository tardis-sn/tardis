import argparse
import os
import time


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
        os.system(test_command.format(args.yaml_filepath, args.atomic_dataset))
        time.sleep(20)
