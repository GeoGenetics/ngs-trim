#################
### FUNCTIONS ###
#################


def allow_remote(url):
    # HTTP(s)
    if url.startswith("http"):
        return storage.http(url)

    # FTP(s)
    if url.startswith("ftp"):
        return storage.ftp(url)

    # If a file URL, just return it
    return url


# Add RAW read type to UNITS
units["read_type_raw"] = [
    ["R1", "R2"] if is_lane_pe(u.sample, u.library, u.lane) else ["R"]
    for u in units.itertuples()
]


#############
### RULES ###
#############
wildcard_constraints:
    read_type_raw="|".join(set(flatten(units.read_type_raw))),


def _get_input_data(wildcards):
    data = get_sample_library_lane_data(
        wildcards.sample, wildcards.library, wildcards.lane
    )
    if wildcards.read_type_raw == "R2":
        return allow_remote(data[1])
    else:
        return allow_remote(data[0])


rule get_fastq_raw:
    input:
        _get_input_data,
    output:
        "results/reads/raw/{sample}_{library}_{lane}_{read_type_raw}.fastq.gz",
    localrule: True
    threads: 1
    run:
        import gzip

        src = Path(input[0])
        dst = Path(output[0])
        if src.suffix != ".gz":
            raise ValueError("Input FASTQ files is not GZip'ed: {}!".format(input[0]))

        with gzip.open(src, "rb") as f:
            f.read(2)
        dst.symlink_to(src.absolute())


##########
### QC ###
##########
rule fastqc_raw:
    input:
        rules.get_fastq_raw.output,
    output:
        html="stats/reads/fastqc/raw/{sample}_{library}_{lane}_{read_type_raw}.html",
        zip="stats/reads/fastqc/raw/{sample}_{library}_{lane}_{read_type_raw}_fastqc.zip",
    log:
        "logs/reads/fastqc/raw/{sample}_{library}_{lane}_{read_type_raw}.log",
    benchmark:
        "benchmarks/reads/fastqc/raw/{sample}_{library}_{lane}_{read_type_raw}.jsonl"
    threads: 4
    resources:
        mem=lambda w, attempt: f"{5* attempt} GiB",
        runtime=lambda w, attempt: f"{2* attempt} h",
    wrapper:
        f"{wrapper_ver}/bio/fastqc"
