#  tests for nomkl-1.0-h5ca1d4c_0 (this is a generated file);
print('===== testing package: nomkl-1.0-h5ca1d4c_0 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import subprocess
import json

args = ("conda", "list", "mkl", "--json")

mkl_list = json.loads(subprocess.check_output(args))
assert "mkl" not in mkl_list#  --- run_test.py (end) ---

print('===== nomkl-1.0-h5ca1d4c_0 OK =====');
