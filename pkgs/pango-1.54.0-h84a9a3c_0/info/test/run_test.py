#  tests for pango-1.54.0-h84a9a3c_0 (this is a generated file);
print('===== testing package: pango-1.54.0-h84a9a3c_0 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import gi
gi.require_version('PangoCairo', '1.0')
from gi.repository import PangoCairo
import sys

fontmap = PangoCairo.FontMap.get_default()
fontmap.list_families()
#  --- run_test.py (end) ---

print('===== pango-1.54.0-h84a9a3c_0 OK =====');
