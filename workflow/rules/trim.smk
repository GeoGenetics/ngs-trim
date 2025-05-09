#################
### FUNCTIONS ###
#################


# Add TRIM read type
units["read_type_trim"] = [
    flatten(
        [
            get_read_type_trim(read_type_map)
            for read_type_map in get_read_type_map(u.sample, u.library, u.lane)
        ]
    )
    for u in units.itertuples()
]


def get_adapters(wildcards):
    adapters = units.loc[
        (wildcards.sample, wildcards.library, wildcards.lane), "adapters"
    ]
    if isinstance(adapters, str):
        return adapters.split(",")

    # If no adapters found, return None
    return None


#############
### RULES ###
#############


wildcard_constraints:
    read_type_trim="|".join(set(flatten(units.read_type_trim))),


rule adapters_to_file:
    output:
        fas=temp("temp/reads/trim/adapters/{sample}_{library}_{lane}.adapters.fas"),
    log:
        "logs/reads/trim/adapters/{sample}_{library}_{lane}.log",
    params:
        adapters=get_adapters,
    run:
        with open(output.fas, "w") as fp:
            for i in range(len(params.adapters)):
                print(">adapter_R" + str(i + 1), file=fp)
                print(params.adapters[i], file=fp)


if config["trim"]["tool"] == "cutadapt":

    rule cutadapt_fastq_pe:
        input:
            lambda w: expand(
                rules.get_fastq_raw.output,
                read_type_raw=units.loc[w.sample, w.library, w.lane].read_type_raw,
                allow_missing=True,
            ),
        output:
            fastq1="results/reads/trim/{sample}_{library}_{lane}_R1.fastq.gz",
            fastq2="results/reads/trim/{sample}_{library}_{lane}_R2.fastq.gz",
            qc="stats/reads/trim/{sample}_{library}_{lane}_pe.qc.txt",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_pe.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_pe.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (" -a {} -A {} ".format(*get_adapters(w)) if get_adapters(w) else " "),
        priority: 10
        threads: 10
        resources:
            mem=lambda w, attempt: f"{15* attempt} GiB",
            runtime=lambda w, attempt: f"{5* attempt} h",
        wrapper:
            f"{wrapper_ver}/bio/cutadapt/pe"

    rule cutadapt_fastq_se:
        input:
            lambda w: expand(
                rules.get_fastq_raw.output,
                read_type_raw=units.loc[w.sample, w.library, w.lane].read_type_raw,
                allow_missing=True,
            ),
        output:
            fastq="results/reads/trim/{sample}_{library}_{lane}_R.fastq.gz",
            qc="stats/reads/trim/{sample}_{library}_{lane}_se.qc.txt",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_se.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_se.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (" -a {} ".format(*get_adapters(w)) if get_adapters(w) else " "),
        priority: 10
        threads: 10
        resources:
            mem=lambda w, attempt: f"{15* attempt} GiB",
            runtime=lambda w, attempt: f"{5* attempt} h",
        wrapper:
            f"{wrapper_ver}/bio/cutadapt/se"

elif config["trim"]["tool"] == "adapterremoval":

    def _collapsed_input():
        return {
            "collapsed": "results/reads/trim/{sample}_{library}_{lane}_collapsed.fastq.gz",
            "collapsed_trunc": "results/reads/trim/{sample}_{library}_{lane}_collapsedtrunc.fastq.gz",
        }

    rule adapterremoval_fastq_pe:
        input:
            sample=lambda w: expand(
                rules.get_fastq_raw.output,
                read_type_raw=units.loc[w.sample, w.library, w.lane].read_type_raw,
                allow_missing=True,
            ),
        output:
            **_collapsed_input() if is_activated("collapse") else {},
            fq1="results/reads/trim/{sample}_{library}_{lane}_R1.fastq.gz",
            fq2="results/reads/trim/{sample}_{library}_{lane}_R2.fastq.gz",
            singleton="results/reads/trim/{sample}_{library}_{lane}_singleton.fastq.gz",
            discarded="results/reads/trim/{sample}_{library}_{lane}_discarded.fastq.gz",
            settings="stats/reads/trim/{sample}_{library}_{lane}_pe.settings",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_pe.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_pe.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (
                " --adapter1 {} --adapter2 {} ".format(*get_adapters(w))
                if get_adapters(w)
                else " "
            )
            + (config["collapse"]["params"] if is_activated("collapse") else " "),
        priority: 10
        threads: 10
        resources:
            mem=lambda w, attempt: f"{15* attempt} GiB",
            runtime=lambda w, attempt: f"{3* attempt} h",
        wrapper:
            f"{wrapper_ver}/bio/adapterremoval"

    use rule adapterremoval_fastq_pe as adapterremoval_fastq_se with:
        output:
            fq="results/reads/trim/{sample}_{library}_{lane}_R.fastq.gz",
            discarded="results/reads/trim/{sample}_{library}_{lane}_discarded.fastq.gz",
            settings="stats/reads/trim/{sample}_{library}_{lane}_se.settings",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_se.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_se.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (" --adapter1 {} ".format(*get_adapters(w)) if get_adapters(w) else " "),

elif config["trim"]["tool"] == "fastp":

    rule fastp_fastq_pe:
        input:
            sample=lambda w: expand(
                rules.get_fastq_raw.output,
                read_type_raw=units.loc[w.sample, w.library, w.lane].read_type_raw,
                allow_missing=True,
            ),
        output:
            trimmed=expand(
                "results/reads/trim/{sample}_{library}_{lane}_{R}.fastq.gz",
                R=["R1", "R2"],
                allow_missing=True,
            ),
            unpaired="results/reads/trim/{sample}_{library}_{lane}_singleton.fastq.gz",
            merged=(
                "results/reads/trim/{sample}_{library}_{lane}_collapsed.fastq.gz"
                if is_activated("collapse")
                else []
            ),
            failed="results/reads/trim/{sample}_{library}_{lane}_discarded.fastq.gz",
            html="stats/reads/trim/{sample}_{library}_{lane}_pe.html",
            json="stats/reads/trim/{sample}_{library}_{lane}_pe.json",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_pe.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_pe.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (
                " --adapter_sequence {} --adapter_sequence_r2 {} ".format(
                        *get_adapters(w)
                    )
                    if get_adapters(w)
                else " "
            )
            + (config["collapse"]["params"] if is_activated("collapse") else " "),
        priority: 10
        threads: 10
        resources:
            mem=lambda w, attempt: f"{20* attempt} GiB",
            runtime=lambda w, attempt: f"{30* attempt} m",
        wrapper:
            f"{wrapper_ver}/bio/fastp"

    use rule fastp_fastq_pe as fastp_fastq_se with:
        output:
            trimmed="results/reads/trim/{sample}_{library}_{lane}_R.fastq.gz",
            failed="results/reads/trim/{sample}_{library}_{lane}_discarded.fastq.gz",
            html="stats/reads/trim/{sample}_{library}_{lane}_se.html",
            json="stats/reads/trim/{sample}_{library}_{lane}_se.json",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_se.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_se.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (
                " --adapter_sequence {} ".format(*get_adapters(w))
                if get_adapters(w)
                else " "
            ),

elif config["trim"]["tool"] == "trimmomatic":

    rule trimmomatic_fastq_pe:
        input:
            unpack(
                lambda w: dict(
                    zip(
                        ["r1", "r2"],
                        expand(
                            rules.get_fastq_raw.output,
                            read_type_raw=get_read_type_raw(
                                w.sample, w.library, w.lane
                            ),
                            allow_missing=True,
                        ),
                    )
                )
            ),
            adapt=rules.adapters_to_file.output.fas,
        output:
            r1="results/reads/trim/{sample}_{library}_{lane}_R1.fastq.gz",
            r2="results/reads/trim/{sample}_{library}_{lane}_R2.fastq.gz",
            r1_unpaired="results/reads/trim/{sample}_{library}_{lane}_singleton1.fastq.gz",
            r2_unpaired="results/reads/trim/{sample}_{library}_{lane}_singleton2.fastq.gz",
            trim_log="stats/reads/trim/{sample}_{library}_{lane}_pe.log",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_pe.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_pe.jsonl"
        params:
            extra=lambda w, output: f"-trimlog {output.trim_log}",
            trimmer=lambda w, input: [
                config["trim"]["params"].format(ADAPTER_FASTA=input.adapt)
            ],
        priority: 10
        threads: 10
        resources:
            mem=lambda w, attempt: f"{15* attempt} GiB",
            runtime=lambda w, attempt: f"{2* attempt} h",
        wrapper:
            f"{wrapper_ver}/bio/trimmomatic/pe"

    rule trimmomatic_fastq_se:
        input:
            sample=lambda w: expand(
                rules.get_fastq_raw.output,
                read_type_raw=units.loc[w.sample, w.library, w.lane].read_type_raw,
                allow_missing=True,
            ),
            adapt=rules.adapters_to_file.output.fas,
        output:
            fq="results/reads/trim/{sample}_{library}_{lane}_R.fastq.gz",
            trim_log="stats/reads/trim/{sample}_{library}_{lane}_se.log",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_se.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_se.jsonl"
        params:
            extra=lambda w, output: f"-trimlog {output.trim_log}",
            trimmer=lambda w, input: [
                config["trim"]["params"].format(ADAPTER_FASTA=input.adapt)
            ],
        priority: 10
        threads: 20
        resources:
            mem=lambda w, attempt: f"{15* attempt} GiB",
            runtime=lambda w, attempt: f"{2* attempt} h",
        wrapper:
            f"{wrapper_ver}/bio/trimmomatic/se"

elif config["trim"]["tool"] == "bbduk":

    rule bbduk_fastq_pe:
        input:
            sample=lambda w: expand(
                rules.get_fastq_raw.output,
                read_type_raw=units.loc[w.sample, w.library, w.lane].read_type_raw,
                allow_missing=True,
            ),
        output:
            trimmed=expand(
                "results/reads/trim/{sample}_{library}_{lane}_{R}.fastq.gz",
                R=["R1", "R2"],
                allow_missing=True,
            ),
            singleton="results/reads/trim/{sample}_{library}_{lane}_singleton.fastq.gz",
            discarded="results/reads/trim/{sample}_{library}_{lane}_discarded.fastq.gz",
            stats="stats/reads/trim/{sample}_{library}_{lane}_pe.stats.txt",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_pe.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_pe.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (" literal={},{} ".format(*get_adapters(w)) if get_adapters(w) else " "),
        priority: 10
        threads: 10
        resources:
            mem=lambda w, attempt: f"{20* attempt} GiB",
            runtime=lambda w, attempt: f"{30* attempt} m",
        wrapper:
            f"{wrapper_ver}/bio/bbtools/bbduk"

    use rule bbduk_fastq_pe as bbduk_fastq_se with:
        output:
            trimmed="results/reads/trim/{sample}_{library}_{lane}_R.fastq.gz",
            discarded="results/reads/trim/{sample}_{library}_{lane}_discarded.fastq.gz",
            stats="stats/reads/trim/{sample}_{library}_{lane}_se.stats.txt",
        log:
            "logs/reads/trim/{sample}_{library}_{lane}_se.log",
        benchmark:
            "benchmarks/reads/trim/{sample}_{library}_{lane}_se.jsonl"
        params:
            extra=lambda w: config["trim"]["params"]
            + (" literal={} ".format(*get_adapters(w)) if get_adapters(w) else " "),


##########
### QC ###
##########


use rule fastqc_raw as fastqc_trim with:
    input:
        fq="results/reads/trim/{sample}_{library}_{lane}_{read_type_trim}.fastq.gz",
    output:
        html="stats/reads/fastqc/trim/{sample}_{library}_{lane}_{read_type_trim}.html",
        zip="stats/reads/fastqc/trim/{sample}_{library}_{lane}_{read_type_trim}_fastqc.zip",
    log:
        "logs/reads/fastqc/trim/{sample}_{library}_{lane}_{read_type_trim}.log",
    benchmark:
        "benchmarks/reads/fastqc/trim/{sample}_{library}_{lane}_{read_type_trim}.jsonl"
