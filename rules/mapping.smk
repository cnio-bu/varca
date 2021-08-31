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
        mem_mb = get_resource("trim_reads","mem"),
        walltime = get_resource("trim_reads","walltime")
    wrapper:
        "0.77.0/bio/trimmomatic/se"

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
        mem_mb = get_resource("trim_reads","mem"),
        walltime = get_resource("trim_reads","walltime")
    wrapper:
        "0.77.0/bio/trimmomatic/pe"

idx_cmd = "bwa index {input} > {log.out} 2> {log.err}"
rule bwa_idx_genome:
    shadow:"shallow"
    input:
        config["ref"]["genome"]
    output:
        f"{config['ref']['genome']}.amb",
        f"{config['ref']['genome']}.ann",
        f"{config['ref']['genome']}.bwt",
        f"{config['ref']['genome']}.pac",
        f"{config['ref']['genome']}.sa"
    threads: get_resource("bwa_idx_genome","threads")
    resources:
        mem = get_resource("bwa_idx_genome","mem"),
        walltime = get_resource("bwa_idx_genome","walltime")
    log:
        f"{LOGDIR}/bwa_idx_genome/bwa_idx_genome.log"
    benchmark:
        f"{LOGDIR}/bwa_idx_genome/bwa_idx_genome.bmk"
    wrapper:
        "0.77.0/bio/bwa/index"

rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=f"{config['ref']['genome']}.bwt"
    output:
        temp(f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam")
    log:
        f"{LOGDIR}/bwa_mem/{{sample}}-{{unit}}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate"
    threads: get_resource("map_reads","threads")
    resources:
        mem = get_resource("map_reads","mem"),
        walltime = get_resource("map_reads","walltime")
    wrapper:
        "0.77.0/bio/bwa/mem"

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
        mem_mb = get_resource("mark_duplicates","mem"),
        walltime = get_resource("mark_duplicates","walltime")
    params:
        config["params"]["picard"]["MarkDuplicates"] + " -Xmx{}m".format(get_resource("mark_duplicates","mem"))
    wrapper:
        "0.77.0/bio/picard/markduplicates"

rule genome_faidx:
    input:
        config["ref"]["genome"]
    output:
        f"{config['ref']['genome']}.fai"
    log:
        f"{LOGDIR}/genome_faidx/genome_faidx.log"
    threads: get_resource("genome_faidx","threads")
    resources:
        mem = get_resource("genome_faidx","mem"),
        walltime = get_resource("genome_faidx","walltime")
    wrapper:
        "0.77.0/bio/samtools/faidx"

rule obtain_recal_table:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        ref_idx=f"{config['ref']['genome']}.fai",
        known=config["ref"]["known-variants"],
        known_idx=f"{config['ref']['known-variants']}.tbi"
    output:
        recal_table=f"{OUTDIR}/recal/{{sample}}-{{unit}}.grp"
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        f"{LOGDIR}/gatk/bqsr/{{sample}}-{{unit}}.log"
    threads: get_resource("recalibrate_base_qualities","threads")
    resources:
        mem_mb = get_resource("recalibrate_base_qualities","mem"),
        walltime = get_resource("recalibrate_base_qualities","walltime")
    wrapper:
        "0.77.0/bio/gatk/baserecalibrator"

rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        ref_idx=f"{config['ref']['genome']}.fai",
        recal_table=f"{OUTDIR}/recal/{{sample}}-{{unit}}.grp"
    output:
        bam=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam",
        bai=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bai"
    params:
        extra=""
    log:
        f"{LOGDIR}/gatk/gatk_applybqsr/{{sample}}-{{unit}}.log"
    resources:
        mem_mb = get_resource("recalibrate_base_qualities","mem"),
        walltime = get_resource("recalibrate_base_qualities","walltime")
    wrapper:
        "0.77.0/bio/gatk/applybqsr"

rule samtools_index:
    input:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem = get_resource("samtools_index","mem"),
        walltime = get_resource("samtools_index","walltime")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    wrapper:
        "0.77.0/bio/samtools/index"

rule samtools_index_sorted:
    input:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam"
    output:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem = get_resource("samtools_index","mem"),
        walltime = get_resource("samtools_index","walltime")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    wrapper:
        "0.77.0/bio/samtools/index"

rule index_known_variants:
    input:
        file=config['ref']['known-variants']
    output:
        index=f"{config['ref']['known-variants']}.tbi"
    resources:
        mem = get_resource("index_known_variants","mem"),
        walltime = get_resource("index_known_variants","walltime")
    conda:
        "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/index_known_variants.log"
    shell:"""
        gatk IndexFeatureFile -F {input.file} -O {output.index}
    """
