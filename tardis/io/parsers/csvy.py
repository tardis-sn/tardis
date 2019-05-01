import yaml
import pandas as pd

YAML_DELIMITER = '---'
def load_csvy(fname):
    with open(fname) as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert line.strip() == YAML_DELIMITER, 'First line of csvy file is not \'---\''
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise Exception('End YAML_DELIMITER not found') 
        yaml_dict = yaml.load(''.join(yaml_lines[1:-1]))
        data = pd.read_csv(fname, skiprows=yaml_end_ind - 1)

    return yaml_dict, data
