import time
import subprocess


if __name__ == "__main__":
    while True:
        subprocess.call([
            "python", "setup.py", "test", "--test-path=tardis/tests/test_util.py",
        ])
        time.sleep(20)
