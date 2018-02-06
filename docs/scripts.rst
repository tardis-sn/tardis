***************
Running scripts
***************

Convert CMFGEN files to TARDIS file format
==========================================

This script takes a CMFGEN file as an input , and an output path to save converted files in new TARDIS file format. 
CSV file contains abundances of both elements and isotopes.
DAT file contains values of velocity, density, electron_density and temperature.  

Format of command - :code:`cmfgen2tardis /path/to/input_file path/to/output/`  

Example, for how to use-  


.. code-block:: bash

   $ cmfgen2tardis tardis_example/DDC15_SN_HYDRO_DATA_0.976d tardis_example/
