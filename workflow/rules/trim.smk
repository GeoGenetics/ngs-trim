
#################
### FUNCTIONS ###
#################

# Add TRIM read type
units["read_type_trim"] = [flatten([get_read_type_trim(read_type_map) for read_type_map in get_read_type_map(u.sample, u.library, u.lane)]) for u in units.itertuples()]


def get_fastqs_trimmed(format = "fastq"):
    return expand("results/reads/trim/{trimmer}/{sample}_{library}_{lane}_{read_type_trim}.{format}.gz",
                   trimmer = config["reads"]["trim"]["tool"],
                   format = format,
                   allow_missing = True)



#############
### RULES ###
#############

wildcard_constraints:
    read_type_trim = "|".join(set(flatten(units.read_type_trim))),
    trimmer = config["reads"]["trim"]["tool"]



rule adapters_to_file:
    output:
        fas = temp("temp/reads/trim/adapters/{sample}_{library}_{lane}.adapters.fas"),
    log:
        "logs/reads/trim/adapters/{sample}_{library}_{lane}.log"
    params:
        adapters = get_adapters,
    run:
        with open(output.fas, 'w') as fp:
            for i in range(len(params.adapters)):
                print(">adapter_R" + str(i+1), file=fp)
                print(params.adapters[i], file=fp)



rule cutadapt_fastq_pe:
    input:
        lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        fastq1 = "results/reads/trim/cutadapt/{sample}_{library}_{lane}_R1.fastq.gz",
        fastq2 = "results/reads/trim/cutadapt/{sample}_{library}_{lane}_R2.fastq.gz",
        qc = "stats/reads/trim/cutadapt/{sample}_{library}_{lane}_pe.qc.txt",
    log:
        "logs/reads/trim/cutadapt/{sample}_{library}_{lane}_pe.log"
    benchmark:
        "benchmarks/reads/trim/cutadapt/{sample}_{library}_{lane}_pe.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args = ["--version", "-j", "--cores", "-a", "-A", "-o", "-p"]) + (" -a {} -A {} ".format(*get_adapters(w)) if get_adapters(w) else " "),
    priority: 10
    threads: 20
    resources:
        mem = lambda w, attempt: f"{15 * attempt} GiB",
        runtime = lambda w, attempt: f"{5 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/cutadapt/pe"

rule cutadapt_fastq_se:
    input:
        lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        fastq = "results/reads/trim/cutadapt/{sample}_{library}_{lane}_R.fastq.gz",
        qc = "stats/reads/trim/cutadapt/{sample}_{library}_{lane}_se.qc.txt",
    log:
        "logs/reads/trim/cutadapt/{sample}_{library}_{lane}_se.log"
    benchmark:
        "benchmarks/reads/trim/cutadapt/{sample}_{library}_{lane}_se.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args = ["--version", "-j", "--cores", "-a", "-A", "-o", "-p"]) + (" -a {} ".format(*get_adapters(w)) if get_adapters(w) else " "),
    priority: 10
    threads: 20
    resources:
        mem = lambda w, attempt: f"{15 * attempt} GiB",
        runtime = lambda w, attempt: f"{5 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/cutadapt/se"



def _collapsed_input():
    return {"collapsed": "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_collapsed.fastq.gz",
            "collapsed_trunc": "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_collapsedtrunc.fastq.gz"}

rule adapterremoval_fastq_pe:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        fq1 = "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_R1.fastq.gz",
        fq2 = "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_R2.fastq.gz",
        singleton = "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_singleton.fastq.gz",
        **_collapsed_input() if is_activated("reads/collapse") else {},
        discarded = "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_discarded.fastq.gz",
        settings = "stats/reads/trim/adapterremoval/{sample}_{library}_{lane}_pe.settings",
    log:
        "logs/reads/trim/adapterremoval/{sample}_{library}_{lane}_pe.log"
    benchmark:
        "benchmarks/reads/trim/adapterremoval/{sample}_{library}_{lane}_pe.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args=["--adapter1", "--adapter2"]) + (" --adapter1 {} --adapter2 {} ".format(*get_adapters(w)) if get_adapters(w) else " ") + (config["reads"]["collapse"]["params"] if is_activated("reads/collapse") else " "),
    priority: 10
    threads: 10
    resources:
        mem = lambda w, attempt: f"{15 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/adapterremoval"

rule adapterremoval_fastq_se:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        fq = "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_R.fastq.gz",
        discarded = "results/reads/trim/adapterremoval/{sample}_{library}_{lane}_discarded.fastq.gz",
        settings = "stats/reads/trim/adapterremoval/{sample}_{library}_{lane}_se.settings",
    log:
        "logs/reads/trim/adapterremoval/{sample}_{library}_{lane}_se.log"
    benchmark:
        "benchmarks/reads/trim/adapterremoval/{sample}_{library}_{lane}_se.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args=["--adapter1", "--adapter2"]) + (" --adapter1 {} ".format(*get_adapters(w)) if get_adapters(w) else " "),
    priority: 10
    threads: 10
    resources:
        mem = lambda w, attempt: f"{15 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/adapterremoval"



rule fastp_fastq_pe:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        trimmed = expand("results/reads/trim/fastp/{sample}_{library}_{lane}_{R}.fastq.gz", R = ["R1","R2"], allow_missing=True),
        unpaired = "results/reads/trim/fastp/{sample}_{library}_{lane}_singleton.fastq.gz",
        merged = "results/reads/trim/fastp/{sample}_{library}_{lane}_collapsed.fastq.gz" if is_activated("reads/collapse") else [],
        failed = "results/reads/trim/fastp/{sample}_{library}_{lane}_discarded.fastq.gz",
        html = "stats/reads/trim/fastp/{sample}_{library}_{lane}_pe.html",
        json = "stats/reads/trim/fastp/{sample}_{library}_{lane}_pe.json",
    log:
        "logs/reads/trim/fastp/{sample}_{library}_{lane}_pe.log"
    benchmark:
        "benchmarks/reads/trim/fastp/{sample}_{library}_{lane}_pe.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args = ["--adapter_sequence", "--adapter_sequence_r2"]) + (" --adapter_sequence {} --adapter_sequence_r2 {} ".format(*get_adapters(w)) if get_adapters(w) else " ") + (config["reads"]["collapse"]["params"] if is_activated("reads/collapse") else " "),
    priority: 10
    threads: 10
    resources:
        mem = lambda w, attempt: f"{30 * attempt} GiB",
        runtime = lambda w, attempt: f"{30 * attempt} m",
    wrapper:
        wrapper_ver + "/bio/fastp"

rule fastp_fastq_se:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        trimmed = "results/reads/trim/fastp/{sample}_{library}_{lane}_R.fastq.gz",
        failed = "results/reads/trim/fastp/{sample}_{library}_{lane}_discarded.fastq.gz",
        html = "stats/reads/trim/fastp/{sample}_{library}_{lane}_se.html",
        json = "stats/reads/trim/fastp/{sample}_{library}_{lane}_se.json",
    log:
        "logs/reads/trim/fastp/{sample}_{library}_{lane}_se.log"
    benchmark:
        "benchmarks/reads/trim/fastp/{sample}_{library}_{lane}_se.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args = ["--adapter_sequence", "--adapter_sequence_r2"]) + (" --adapter_sequence {} ".format(*get_adapters(w)) if get_adapters(w) else " "),
    priority: 10
    threads: 10
    resources:
        mem = lambda w, attempt: f"{30 * attempt} GiB",
        runtime = lambda w, attempt: f"{30 * attempt} m",
    wrapper:
        wrapper_ver + "/bio/fastp"



rule trimmomatic_fastq_pe:
    input:
        unpack(lambda w: dict(zip(["r1","r2"], expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True)))),
        adapt = rules.adapters_to_file.output.fas,
    output:
        r1 = "results/reads/trim/trimmomatic/{sample}_{library}_{lane}_R1.fastq.gz",
        r2 = "results/reads/trim/trimmomatic/{sample}_{library}_{lane}_R2.fastq.gz",
        r1_unpaired = "results/reads/trim/trimmomatic/{sample}_{library}_{lane}_singleton1.fastq.gz",
        r2_unpaired = "results/reads/trim/trimmomatic/{sample}_{library}_{lane}_singleton2.fastq.gz",
        trim_log = "stats/reads/trim/trimmomatic/{sample}_{library}_{lane}_pe.log",
    log:
        "logs/reads/trim/trimmomatic/{sample}_{library}_{lane}_pe.log"
    benchmark:
        "benchmarks/reads/trim/trimmomatic/{sample}_{library}_{lane}_pe.tsv"
    params:
        extra = lambda w, output: f"-trimlog {output.trim_log}",
        trimmer = lambda w, input: [config["reads"]["trim"]["params"].format(ADAPTER_FASTA = input.adapt)],
    priority: 10
    threads: 20
    resources:
        mem = lambda w, attempt: f"{15 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/trimmomatic/pe"

rule trimmomatic_fastq_se:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
        adapt = rules.adapters_to_file.output.fas,
    output:
        fq = "results/reads/trim/trimmomatic/{sample}_{library}_{lane}_R.fastq.gz",
        trim_log = "stats/reads/trim/trimmomatic/{sample}_{library}_{lane}_se.log",
    log:
        "logs/reads/trim/trimmomatic/{sample}_{library}_{lane}_se.log"
    benchmark:
        "benchmarks/reads/trim/trimmomatic/{sample}_{library}_{lane}_se.tsv"
    params:
        extra = lambda w, output: f"-trimlog {output.trim_log}",
        trimmer = lambda w, input: [config["reads"]["trim"]["params"].format(ADAPTER_FASTA = input.adapt)],
    priority: 10
    threads: 20
    resources:
        mem = lambda w, attempt: f"{15 * attempt} GiB",
        runtime = lambda w, attempt: f"{2 * attempt} h",
    wrapper:
        wrapper_ver + "/bio/trimmomatic/se"



rule bbduk_fastq_pe:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        trimmed = expand("results/reads/trim/bbduk/{sample}_{library}_{lane}_{R}.fastq.gz", R = ["R1","R2"], allow_missing=True),
        singleton = "results/reads/trim/bbduk/{sample}_{library}_{lane}_singleton.fastq.gz",
        discarded = "results/reads/trim/bbduk/{sample}_{library}_{lane}_discarded.fastq.gz",
        stats = "stats/reads/trim/bbduk/{sample}_{library}_{lane}_pe.stats.txt",
    log:
        "logs/reads/trim/bbduk/{sample}_{library}_{lane}_pe.log"
    benchmark:
        "benchmarks/reads/trim/bbduk/{sample}_{library}_{lane}_pe.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args=["literal"]) + (" literal={},{} ".format(*get_adapters(w)) if get_adapters(w) else " "),
    priority: 10
    threads: 10
    resources:
        mem = lambda w, attempt: f"{20 * attempt} GiB",
        runtime = lambda w, attempt: f"{30 * attempt} m",
    wrapper:
        wrapper_ver + "/bio/bbtools/bbduk"

rule bbduk_fastq_se:
    input:
        sample = lambda w: expand(rules.get_fastq_raw.output, read_type_raw = get_read_type_raw(w.sample, w.library, w.lane), allow_missing=True),
    output:
        trimmed = "results/reads/trim/bbduk/{sample}_{library}_{lane}_R.fastq.gz",
        discarded = "results/reads/trim/bbduk/{sample}_{library}_{lane}_discarded.fastq.gz",
        stats = "stats/reads/trim/bbduk/{sample}_{library}_{lane}_se.stats.txt",
    log:
        "logs/reads/trim/bbduk/{sample}_{library}_{lane}_se.log"
    benchmark:
        "benchmarks/reads/trim/bbduk/{sample}_{library}_{lane}_se.tsv"
    params:
        extra = lambda w: check_cmd(config["reads"]["trim"]["params"], forbidden_args=["literal"]) + (" literal={} ".format(*get_adapters(w)) if get_adapters(w) else " "),
    priority: 10
    threads: 10
    resources:
        mem = lambda w, attempt: f"{20 * attempt} GiB",
        runtime = lambda w, attempt: f"{30 * attempt} m",
    wrapper:
        wrapper_ver + "/bio/bbtools/bbduk"



##########
### QC ###
##########

use rule fastqc_raw as fastqc_trim with:
    input:
        fq = get_fastqs_trimmed(),
    output:
        html = "stats/reads/fastqc_trim/{sample}_{library}_{lane}_{read_type_trim}.html",
        zip = "stats/reads/fastqc_trim/{sample}_{library}_{lane}_{read_type_trim}_fastqc.zip"
    log:
        "logs/reads/fastqc_trim/{sample}_{library}_{lane}_{read_type_trim}.log"
    benchmark:
        "benchmarks/reads/fastqc_trim/{sample}_{library}_{lane}_{read_type_trim}.tsv"
