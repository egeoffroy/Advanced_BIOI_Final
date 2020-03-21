__author__ = "alvaro barbeira"

import pandas

from .. import Utilities

def try_parse(string, fail=None):
    try:
        return float(string)
    except Exception:
        return fail


def skip_na(key, value):
    skip = (not value or value == "NA")
    return skip

def load_data(path, key_name, value_name, white_list=None, numeric=False, should_skip=None, value_conversion=None, key_filter=None):
    data = {}
    c_key=None
    c_value=None
    for i, line in Utilities.iterate_file(path):
        if i==0:
            header = line.strip().split()
            c_key = header.index(key_name)
            c_value = header.index(value_name)
            continue

        comps = line.strip().split()
        key = comps[c_key]
        if white_list and not key in white_list:
            continue
        if key_filter and key_filter(key):
            continue

        value = comps[c_value]
        if should_skip and should_skip(key, value):
            continue

        if value_conversion:
            value = value_conversion(value)
        elif numeric:
            value = try_parse(value)

        if value:
            data[key] = value

    return data

def to_data_frame(data, keys, key_column, value_column, to_numeric=None):
    ids = [x for x in keys if x in data]
    data = [(x, data[x]) for x in ids]
    if len(data) == 0:
        return pandas.DataFrame({key_column: [], value_column: []})
    data = Utilities.to_dataframe(data, [key_column, value_column], to_numeric)
    return data