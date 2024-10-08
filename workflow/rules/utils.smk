import os
import pandas as pd
from typing import List, Dict
from collections import defaultdict
from snakemake.io import Namedlist


#################
### FUNCTIONS ###
#################

### General

def flatten(list_of_lists: List) -> List:
    """Flatten an irregular list of lists recursively

    https://stackoverflow.com/a/53778278

    :param list_of_lists: A list of lists
    :return result: A string that has been flattened from a list of lists
    """
    result = list()
    for i in list_of_lists:
        if isinstance(i, list):
            result.extend(flatten(i))
        else:
            result.append(str(i))
    return result


def expand_pandas(string: List, df: pd.DataFrame, allow_missing=False) -> List:
    """Expand string following columns in the dataframe"""
    return set(
        flatten(
            [
                expand(string, **row._asdict(), allow_missing=allow_missing)
                for row in df.itertuples(False)
            ]
        )
    )


def to_dict(keys: List, values: List) -> Dict:
    """Convert two lists (keys and values) to a dictionary

    Even though it is a simple function, it makes the code more readable.
    """
    if len(keys) != len(values):
        raise ValueError("")

    out_dict = defaultdict(list)
    for key, val in list(zip(keys, values)):
        out_dict[key].append(val)

    return out_dict


### Config

def _item_or_sample(row, item):
    i = getattr(row, item, None)
    if pd.isnull(i):
        return getattr(row, "sample")
    return i


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


### Samples

def get_rule_stats(rule_name):
    r = re.compile("^stats/")
    return set(filter(r.match, getattr(rules, rule_name).output))


def get_samples(seq_type=".", material=".", all=True):
    # Get state where all samples are TRUE
    bool_true = units["sample"].str.match(".")

    type_cond = (
        units.seq_type.str.match(seq_type, case=False)
        if "seq_type" in units.columns
        else bool_true
    )
    material_cond = (
        units.material.str.match(material, case=False)
        if "material" in units.columns
        else bool_true
    )
    tot_cond = (type_cond & material_cond).groupby(level="sample")
    if all:
        tot_cond = tot_cond.all()

    return [idx for idx, bool in tot_cond.items() if bool]


def get_groups():
    return samples["group"].unique()


def get_group_samples(group):
    return samples.loc[samples["group"] == group]["sample"]


### Units


def get_libraries():
    return units["library"].unique()


def get_lanes():
    return units["lane"].unique()


def is_units_align(units):
    return units.data.str.endswith(".cram")


def get_units_pe(reverse=False):
    units_pe = units.seq_type.isin(["pe", "hic"])

    if reverse:
        return units[~units_pe]
    else:
        return units[units_pe]


def get_units_se():
    return get_units_pe(reverse=True)


def get_adapters(wildcards):
    adapters = units.loc[
        (wildcards.sample, wildcards.library, wildcards.lane), "adapters"
    ]
    if isinstance(adapters, str):
        return adapters.split(",")

    # If no adapters found, return None
    return None


def get_sample_libraries(sample):
    return units.loc[sample]["library"].unique()


def get_sample_library_lanes(sample, library):
    return units.loc[(sample, library)].lane.unique()


def get_sample_library_lane_data(sample, library, lane):
    return expand(units.loc[(sample, library, lane)].data, Read=["1", "2"])


def is_lane_pe(sample, library, lane):
    return get_units_pe().index.isin([(sample, library, lane)]).any()


def get_read_type_trim(read_type_map):
    if read_type_map == "pe":
        return ["R1", "R2"]
    elif read_type_map == "se":
        return ["R"]
    else:
        return read_type_map


def get_read_type_map(sample, library, lane):
    read_type_map = list()

    if is_lane_pe(sample, library, lane):
        read_type_map.append("pe")
        if is_activated("reads/trim"):
            trimmer = config["reads"]["trim"]["tool"]

            if trimmer == "adapterremoval":
                read_type_map.append("singleton")
                if is_activated("reads/collapse"):
                    read_type_map.append("collapsed")
                    read_type_map.append("collapsedtrunc")
            elif trimmer == "fastp":
                read_type_map.append("singleton")
                if is_activated("reads/collapse"):
                    read_type_map.append("collapsed")
            elif trimmer == "bbduk":
                read_type_map.append("singleton")
            elif trimmer == "trimmomatic":
                read_type_map.append("singleton1")
                read_type_map.append("singleton2")
            elif trimmer == "cutadapt":
                fall_through = True
            else:
                raise ValueError("Invalid trimmer provided!")
    else:
        read_type_map.append("se")

    return read_type_map
