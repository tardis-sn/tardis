#  tests for harfbuzz-8.5.0-hfac3d4d_0 (this is a generated file);
print('===== testing package: harfbuzz-8.5.0-hfac3d4d_0 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import gi
gi.require_version('HarfBuzz', '0.0')
from gi.repository import HarfBuzz as hb
import sys

if hb.buffer_create () is None:
    sys.exit(1)
#  --- run_test.py (end) ---

print('===== harfbuzz-8.5.0-hfac3d4d_0 OK =====');
