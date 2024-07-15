#  tests for llvmlite-0.42.0-py312hb06c811_1 (this is a generated file);
print('===== testing package: llvmlite-0.42.0-py312hb06c811_1 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import os
from llvmlite.tests import main

# Disable tests for distribution only
# These check for static linkage, which we don't do.
os.environ['LLVMLITE_DIST_TEST'] = ''
main()
#  --- run_test.py (end) ---

print('===== llvmlite-0.42.0-py312hb06c811_1 OK =====');
print("import: 'llvmlite'")
import llvmlite

print("import: 'llvmlite.binding'")
import llvmlite.binding

