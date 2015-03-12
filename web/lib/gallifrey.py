import json
import tornado
import urllib
import yaml
from collections import OrderedDict

def ordered_dump(data, stream=None, Dumper=yaml.Dumper, **kwds):
    class OrderedDumper(Dumper):
        pass
    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())
    OrderedDumper.add_representer(OrderedDict, _dict_representer)
    return yaml.dump(data, stream, OrderedDumper, **kwds)

def make_tree(data):
    subtrees = {}
    details = OrderedDict()
    if isinstance(data,dict):
        for i in data:
            elem = i.split('__',1)
            if len(elem) > 1:
                if elem[0] in data and data[elem[0]]=='y':
                    if elem[0] in subtrees:
                        if data[i]:
                            subtrees[elem[0]][elem[1]] = data[i]
                    else:
                        subtrees[elem[0]] = {}
                        if data[i]:
                            subtrees[elem[0]][elem[1]] = data[i]
            else:
                if data[i]:
                    details[i] = data[i]
        for i in subtrees:
            details[i] = subtrees[i]
    else:
        return data
    return details

def clean(data):
    data = OrderedDict([('tardis_config_version','v1.0')]+data.items())
    try:
        data['model']['abundances']['uniform_abundances'] = data['model']['abundances']['uniform_abundances'].replace('\r\n','\n')
        uniform_abundances = data['model']['abundances']['uniform_abundances'].split('\n')
        for i in uniform_abundances:
            elem = i.split()
            data['model']['abundances'][elem[0]] = elem[1]
        del data['model']['abundances']['uniform_abundances']
    except:
        pass
    return data

def generate_yaml(post,form):
    details = ''
    details = OrderedDict()
    param_list = []
    for i in post.request.body.split('&'):
        param_list.append(i.split('=')[0])
    for f in param_list:
        elem = f.split('-',1)
        if len(elem)>1:
            if elem[0] not in details:
                details[elem[0]] = {}
                details[elem[0]][elem[1]] = post.get_argument(f, default=None, strip=False)
                if elem[0]==elem[1]:
                    details[elem[0]] = post.get_argument(f, default=None, strip=False)
            else:
                details[elem[0]][elem[1]] = post.get_argument(f, default=None, strip=False)
        else:
            details[f] = post.get_argument(f, default=None, strip=False)

    data = OrderedDict()
    for i in details:
        try:
            if i=='abundances_model' or i=='structure_model':
                if 'model' not in data:
                    data['model'] = {}
                data['model'][i.split('_')[0]] = make_tree(details[i])
            else:
                data[i] = make_tree(details[i])
        except:pass
    data = clean(data)
    print
    yaml_file = file('input_data.yml', 'w')
    yaml_file.write('#This input doesn\'t have quantity_range_sampled types broken into Start, Stop and Num, but as comma separated values so it should not work as the input. It would be added in the next push.\n\n')
    a = ordered_dump(data,yaml_file, Dumper=yaml.SafeDumper, default_flow_style=False)
    yaml_file.close()
    post.set_header("Content-Type", "text/yaml")
    post.set_status(201)
    with open('input_data.yml', 'rb') as f:
        post.write(f.read())
