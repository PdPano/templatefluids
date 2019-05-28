'''
Reads and processes an input file

Usage:
    read_input(input_file)
Returns:

'''

from default_values import DEFAULT_OPTIONS
from geometry import Geometry


def read_input(input_file):
    '''
    Takes an open input_file
    returns a dict with the processed options read from file
    and a geometry object
    '''
    local_options = DEFAULT_OPTIONS.copy()
    geometry = Geometry()
    for i, line in enumerate(input_file):
        tokens = _tokenize(line)
        if tokens[0] == '':
            continue

        if len(tokens) != 2:
            raise ValueError(
                "In line {}: Wrong format\nGot: {}".format(i, line)
            )
        key, value = tokens[0].upper(), tokens[1].upper()

        if key in local_options:
            local_options[key] = value
        elif key == "BOUNDARY":
            geometry.add_boundary(value)
        elif key == "BODY":
            geometry.add_body_from_file(value)
        elif key[:4] == "GEO_":
            geometry.add_geometry(key, value)
        else:
            raise ValueError(
                "In line {}: Invalid option {}".format(i, key)
            )
    return local_options, geometry


def _tokenize(line):
    '''
    Separates option name and values, while ignoring comments
    '''
    return [tkn.strip() for tkn in line.split("#")[0].split("=")]
