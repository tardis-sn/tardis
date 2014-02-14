__author__ = 'michi'
import yaml

f = open('default_conf_test/conf_def.yaml')
default = yaml.safe_load(f)
f = open('default_conf_test/conf_tes.yaml')
config = yaml.safe_load(f)
