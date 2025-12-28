Convert CMFGEN files to TARDIS file format
==========================================

This script converts `CMFGEN
<http://kookaburra.phyast.pitt.edu/hillier/web/CMFGEN.htm>`_ input files so that
they can be used as input for TARDIS. The script expects the path to the CMFGEN
files and the destination path for the TARDIS files as parameters
(:code:`cmfgen2tardis /path/to/input_file path/to/output/`):

.. code-block:: bash

   $ cmfgen2tardis tardis_example/DDC15_SN_HYDRO_DATA_0.976d tardis_example/
