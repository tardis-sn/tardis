import argparse
import logging
import subprocess

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
    less_packets = "--less-packets" if args.less_packets else ""
    test_command = [
        "python",
        "setup.py",
        "test",
        "--test-path=tardis/tests/integration_tests/test_integration.py",
        "--args",
        f"--capture=no --integration-tests={args.yaml_filepath} --tardis-refdata={args.tardis_refdata} --remote-data "
        f"{less_packets}",
    ]
    subprocess.call(test_command)
