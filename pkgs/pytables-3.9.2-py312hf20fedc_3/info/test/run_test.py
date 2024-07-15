#  tests for pytables-3.9.2-py312hf20fedc_3 (this is a generated file);
print('===== testing package: pytables-3.9.2-py312hf20fedc_3 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import sys
import os
import tables
import tables._comp_bzip2
# We don't build this one on Windows.
if not sys.platform == "win32":
    import tables._comp_lzo
import tables.hdf5extension
import tables.indexesextension
import tables.linkextension
import tables.lrucacheextension
import tables.tableextension
import tables.utilsextension


if __name__ == "__main__":
    from multiprocessing import freeze_support
    freeze_support()
    tables.test()
#  --- run_test.py (end) ---

print('===== pytables-3.9.2-py312hf20fedc_3 OK =====');
