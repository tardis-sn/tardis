"""
A collection of parsers to read data of different formats. 

Accepted data formats are Arepo snapshots, the Blondin toy-model
format, CSVY files and STELLA files. These are parsed into 
either CSV files, YAML files or Pandas DataFrames, which are 
data formats interpretable by TARDIS. 
"""

from tardis.io.parsers.blondin_toymodel import (
    read_blondin_toymodel,
    convert_blondin_toymodel,
)
