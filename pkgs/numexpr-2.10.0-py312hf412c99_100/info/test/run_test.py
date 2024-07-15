#  tests for numexpr-2.10.0-py312hf412c99_100 (this is a generated file);
print('===== testing package: numexpr-2.10.0-py312hf412c99_100 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import numexpr
import numexpr.interpreter

from multiprocessing import freeze_support

if __name__ == "__main__":
    freeze_support()
    numexpr.test()

#  --- run_test.py (end) ---

print('===== numexpr-2.10.0-py312hf412c99_100 OK =====');
print("import: 'numexpr'")
import numexpr

print("import: 'numexpr.interpreter'")
import numexpr.interpreter

