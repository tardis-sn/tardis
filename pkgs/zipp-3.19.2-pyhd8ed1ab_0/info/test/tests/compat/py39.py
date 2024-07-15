import sys

from jaraco.test.cpython import from_test_support, try_import


os_helper = try_import('os_helper') or from_test_support(
    'FakePath',
    'temp_dir',
)

sys.modules[__name__ + '.os_helper'] = os_helper
