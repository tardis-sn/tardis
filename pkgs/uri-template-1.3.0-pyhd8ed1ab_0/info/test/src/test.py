#!/usr/bin/env python3
"""Run comprehensive test suite."""

import collections
import glob
import json
import os
from typing import List, Optional, Union

from uri_template import URITemplate


def _pass(test: str, result: Optional[str]) -> None:
    print('  PASS: "{test}" == {result}'.format(test=test, result=str(result)))


def _fail(test: str, result: Optional[str], expected: Optional[Union[str, List[str]]]) -> None:
    if (isinstance(expected, list)):
        expected = '\n' + ' or\n'.join(['    "{result}"'.format(result=acceptable) for acceptable in expected])
    print('* FAIL: "{test}" got: "{result}", expected: "{expected}"'.format(test=test, result=str(result),
                                                                            expected=expected))


def _check_result(test: str, result: Optional[str], expected_result: Optional[Union[str, List[str]]]) -> int:
    if (isinstance(expected_result, str)):
        if (expected_result != result):
            _fail(test, result, expected_result)
            return 1
    elif (isinstance(expected_result, list)):
        for possible_result in expected_result:
            if (possible_result == result):
                break
        else:
            _fail(test, result, expected_result)
            return 1
    elif (not expected_result):
        if (result):
            _fail(test, result, expected_result)
            return 1
    else:
        print('** Unknown expected result type: ' + repr(expected_result))
        return 1
    _pass(test, result)
    return 0


def run_tests(test_file_search: str) -> int:
    """Load all tests matching search and run them."""
    fail_count: int = 0
    for test_file_path in sorted(glob.glob(test_file_search)):
        print('Running tests from: ' + test_file_path)
        with open(test_file_path, encoding='utf-8') as test_file:
            test_data = json.load(test_file, object_pairs_hook=collections.OrderedDict)
            for test_set_name in test_data:
                print(test_set_name + ':')
                test_set = test_data[test_set_name]
                for test in test_set.get('testcases', []):
                    expected_result = test[1]
                    try:
                        template = URITemplate(test[0])
                        if (str(template) != test[0]):
                            _fail(test[0], str(template), test[0])
                            fail_count += 1

                        result = template.expand(**test_set['variables'])
                        fail_count += _check_result(test[0], result, expected_result)
                    except Exception:
                        if (expected_result):
                            _fail(test[0], 'Exception', expected_result)
                            fail_count += 1
                        else:
                            _pass(test[0], 'Exception')

                for test in test_set.get('partial_testcases', []):
                    expected_result = test[1]
                    try:
                        template = URITemplate(test[0])
                        if (str(template) != test[0]):
                            _fail(test[0], str(template), test[0])
                            fail_count += 1

                        partial = template.partial(**test_set['partial_variables'])
                        fail_count += _check_result(test[0] + ' partial', str(partial), expected_result)
                        if (2 < len(test)):
                            fail_count += _check_result(test[0] + ' expanded partial == expanded', partial.expand(), template.expand(**test_set['partial_variables']))
                            result = str(partial.expand(**test_set['variables']))
                            fail_count += _check_result(test[0] + ' completed partial', result, test[2])
                            fail_count += _check_result(test[0] + ' completed partial == fully expanded', result, str(template.expand(**test_set['variables'])))
                    except Exception:
                        if (expected_result):
                            _fail(test[0], 'Exception', expected_result)
                            fail_count += 1
                        else:
                            _pass(test[0], 'Exception')

    return fail_count


if ('__main__' == __name__):      # called from the command line
    fail_count = 0
    fail_count += run_tests(os.path.join('tests', '*.json'))
    print('{count} failures'.format(count=fail_count))
    if (fail_count):
        exit(1)
