import pandas as pd
from typing import List


#################
### FUNCTIONS ###
#################

### General


def expand_pd(string: List, df: pd.DataFrame, allow_missing=False) -> List:
    return set(expand(string, zip, **df.to_dict("list"), allow_missing=allow_missing))


### Config


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


### Samples


def get_rule_stats(rule_name):
    r = re.compile("^stats/")
    return set(filter(r.match, getattr(rules, rule_name).output))


### Units


def is_units_align(units):
    return units.data.str.endswith(".cram")


def is_lane_pe(sample, library, lane):
    return units.loc[sample, library, lane].seq_type in ["pe", "hic"]


def get_sample_library_lane_data(sample, library, lane):
    return expand(units.loc[(sample, library, lane)].data, Read=["1", "2"])


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
        if is_activated("trim"):
            trimmer = config["trim"]["tool"]

            if trimmer == "adapterremoval":
                read_type_map.append("singleton")
                if is_activated("collapse"):
                    read_type_map.append("collapsed")
                    read_type_map.append("collapsedtrunc")
            elif trimmer == "fastp":
                read_type_map.append("singleton")
                if is_activated("collapse"):
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
