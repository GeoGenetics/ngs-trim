
#################
### FUNCTIONS ###
#################

def get_read_type_raw(sample, library, lane):
    return ["R1","R2"] if is_lane_pe(sample, library, lane) else ["R"]


# Add RAW read type to UNITS
units["read_type_raw"] = [get_read_type_raw(u.sample, u.library, u.lane) for u in units.itertuples()]



#############
### RULES ###
#############

wildcard_constraints:
    read_type_raw  = "|".join(set(flatten(units.read_type_raw)))


def _get_input(wildcards):
    data = get_sample_library_lane_data(wildcards.sample, wildcards.library, wildcards.lane)
    if wildcards.read_type_raw == "R2":
        return allow_remote(data[1])
    else:
        return allow_remote(data[0])


rule get_fastq_raw:
    input:
        _get_input
    output:
        temp("temp/reads/raw/{sample}_{library}_{lane}_{read_type_raw}.fastq.gz"),
    localrule: True
    threads: 1
    run:
        import gzip
        src = Path(input[0])
        dst = Path(output[0])
        if src.suffix != ".gz":
            raise ValueError("Input FASTQ files is not GZip'ed: {}!".format(input[0]))

        f = gzip.open(src, 'rb')
        f.read(2)
        f.close()
        create_symlink(src, dst)



##########
### QC ###
##########
rule fastqc_raw:
    input:
        rules.get_fastq_raw.output,
    output:
        html = "stats/reads/fastqc_raw/{sample}_{library}_{lane}_{read_type_raw}.html",
        zip = "stats/reads/fastqc_raw/{sample}_{library}_{lane}_{read_type_raw}_fastqc.zip"
    log:
        "logs/reads/fastqc_raw/{sample}_{library}_{lane}_{read_type_raw}.log"
    benchmark:
        "benchmarks/reads/fastqc_raw/{sample}_{library}_{lane}_{read_type_raw}.jsonl"
    params: ""
    threads: 2
    resources:
        # Memory is hard-coded to 250M per thread (https://github.com/bcbio/bcbio-nextgen/issues/2989)
        mem = lambda w, threads: f"{512 * threads} MiB",
        runtime = lambda w, attempt: f"{30 * attempt} m",
    wrapper:
        wrapper_ver + "/bio/fastqc"
