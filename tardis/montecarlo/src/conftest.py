import subprocess
import pytest
import os


def pytest_collect_file(parent, path):
    if path.ext == ".c" and path.basename.startswith("test_"):
        return CTestFile(path, parent)


class CTestFile(pytest.File):

    def collect(self):
        test_exe = os.path.splitext(str(self.fspath))[0]
        test_output = subprocess.check_output(test_exe)
        lines = test_output.split("\n")
        lines = [line.strip() for line in lines]
        lines = [line for line in lines if line.startswith("[")]
        test_results = []
        for line in lines:
            token, data = line.split(" ", 1)
            token = token[1:-1]

            if token in ("PASS", "FAIL"):
                file_name, function_name, line_number = data.split(":")
                test_results.append({"condition": token,
                                     "file_name": file_name,
                                     "function_name": function_name,
                                     "line_number": int(line_number),
                                     "EXP": 'EXP(no data found)',
                                     "GOT": 'GOT(no data found)',
                                     "TST": 'TST(no data found)',
                                     })
            elif token in ("EXP", "GOT", "TST"):
                test_results[-1][token] = data

        for test_result in test_results:
            yield CTestItem(test_result["function_name"], self, test_result)


class CTestItem(pytest.Item):

    def __init__(self, name, parent, test_result):
        super(CTestItem, self).__init__(name, parent)
        self.test_result = test_result

    def runtest(self):
        if self.test_result["condition"] == "FAIL":
            raise CTestException(self, self.name)

    def repr_failure(self, exception):
        if isinstance(exception.value, CTestException):
            return ("Test failed : {TST} at {file_name}:{line_number}\n"
                    "         got: {GOT}\n"
                    "    expected: {EXP}\n".format(**self.test_result))

    def reportinfo(self):
        return self.fspath, self.test_result["line_number"] - 1, self.name


class CTestException(Exception):
    pass

