import pytest
from scripts.parse_data_carsus.main import read_data, ParameterException

# Predefined Data
filename = 'si2_osc_kurucz'
directory_path = '/Users/sashmish/Documents/personal/carsus/atomic/'
    
METADATA = r'^([\w\-\.]+)\s+!([\w ]+)$'
PATTERN1 = r'^([\w\[\]\(\)\/]+)\s+([\w\.]+)\s+([+-]?[\d\.]+)\s+([+-]?[\d\.]+)\s+([+-]?[\d\.]+)\s+([+-]?[\dE\.\+\-]+)\s+([+-]?\d+)\s+([+-]?[\dE\.\+\-]+)\s+([+-]?[\dE\.\+\-]+)\s+([+-]?[\dE\.\+\-]+)'
PATTERN2 = r'^([\w\[\]\(\)\/]+)\s*-([\w\[\]\(\)\/]+)\s+([+-]?[\dE\.\+\-]+)\s+([+-]?[\dE\.\+\-]+)\s+([+-]?[\dE\.\+\-]+)\s+(\d+-\s*\d+)'

columns_meta = ['Value', 'Name']
columns1 = ['level','g','E(cm^-1)','10^15 Hz','eV','Lam(A)','ID','ARAD','C4','C6']
columns2 = ['Transition1','Transition2','f','A','Lam(A)','i-j']
    
@pytest.mark.parametrize("pattern, columns", [(PATTERN1, columns1), (PATTERN2, columns2)])
def test_table(pattern, columns):
    df = read_data(directory_path, filename, pattern, columns)
    assert all(df.columns == columns)

def test_metadata():
    df = read_data(directory_path, filename, METADATA, columns_meta)
    assert all(df.columns == columns_meta)

def test_incorrect_filename():
    with pytest.raises(ParameterException):
        read_data(directory_path, 'si2_osc_kurucz1', METADATA, columns_meta)
