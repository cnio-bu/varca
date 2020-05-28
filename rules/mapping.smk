rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.fastq.gz")
    params:
        extra="",
        **config["params"]["trimmomatic"]["se"]
    log:
        f"{LOGDIR}/trimmomatic/{{sample}}-{{unit}}.log"
    threads: get_resource("trim_reads","threads")
    resources:
        mem = get_resource("trim_reads","mem")
    wrapper:
        "0.35.0/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.fastq.gz"),
        r2=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.fastq.gz"),
        r1_unpaired=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.unpaired.fastq.gz"),
        trimlog=f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        f"{LOGDIR}/trimmomatic/{{sample}}-{{unit}}.log"
    threads: get_resource("trim_reads","threads")
    resources:
        mem = get_resource("trim_reads","mem")
    wrapper:
        "0.35.0/bio/trimmomatic/pe"


rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        temp(f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam")
    log:
        f"{LOGDIR}/bwa_mem/{{sample}}-{{unit}}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: get_resource("map_reads","threads")
    resources:
        mem = get_resource("map_reads","mem")
    wrapper:
        "0.35.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam"
    output:
        bam=temp(f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"),
        metrics=f"{OUTDIR}/qc/dedup/{{sample}}-{{unit}}.metrics.txt"
    log:
        f"{LOGDIR}/picard/dedup/{{sample}}-{{unit}}.log"
    threads: get_resource("mark_duplicates","threads")
    resources:
        mem = get_resource("mark_duplicates","mem")
    params:
        config["params"]["picard"]["MarkDuplicates"] + " -Xmx{}m".format(get_resource("mark_duplicates","mem"))
    wrapper:
        "0.35.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"]
    output:
        bam=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam"
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        f"{LOGDIR}/gatk/bqsr/{{sample}}-{{unit}}.log"
    threads: get_resource("recalibrate_base_qualities","threads")
    resources:
        mem = get_resource("recalibrate_base_qualities","mem")
    wrapper:
        "0.35.0/bio/gatk/baserecalibrator"

rule samtools_index:
    input:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem = get_resource("samtools_index","mem")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    wrapper:
        "0.35.0/bio/samtools/index"
