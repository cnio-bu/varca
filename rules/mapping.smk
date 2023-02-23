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
        "0.79.0/bio/trimmomatic/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.fastq.gz"),
        r2=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.fastq.gz"),
        r1_unpaired=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.unpaired.fastq.gz"),
        trimlog=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.trimlog.txt")
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
        "0.79.0/bio/trimmomatic/pe"

idx_cmd = "bwa index {input} > {log.out} 2> {log.err}"
rule bwa_idx_genome:
    shadow:"shallow"
    input:
        config["ref"]["genome"]
    output:
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".0123",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".amb",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".ann",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".bwt.2bit.64",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".pac"
    threads: get_resource("bwa_idx_genome","threads")
    resources:
        mem_mb = get_resource("bwa_idx_genome","mem"),
        walltime = get_resource("bwa_idx_genome","walltime")
    log:
        f"{LOGDIR}/bwa_idx_genome/bwa_idx_genome.log"
    params:
        prefix=config['ref']['genome_idx']+os.path.basename(config['ref']['genome'])
    benchmark:
        f"{LOGDIR}/bwa_idx_genome/bwa_idx_genome.bmk"
    wrapper:
        "v1.21.6/bio/bwa-mem2/index"

rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=multiext(f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome']), ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        alt=f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".alt" if config['ref']['genome_alt'] != None else []
    output:
        temp(f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam")
    log:
        f"{LOGDIR}/bwa_mem/{{sample}}-{{unit}}.log"
    params:
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    shadow: "shallow"
    threads: get_resource("map_reads","threads")
    resources:
        mem_mb = get_resource("map_reads","mem"),
        walltime = get_resource("map_reads","walltime")
    wrapper:
        "v1.21.6/bio/bwa-mem2/mem"

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
        "0.79.0/bio/picard/markduplicates"

checkpoint genome_faidx:
    input:
        config["ref"]["genome"]
    output:
        f"{config['ref']['genome']}.fai"
    log:
        f"{LOGDIR}/genome_faidx/genome_faidx.log"
    threads: get_resource("genome_faidx","threads")
    resources:
        mem_mb = get_resource("genome_faidx","mem"),
        walltime = get_resource("genome_faidx","walltime")
    wrapper:
        "v1.21.6/bio/samtools/faidx"

rule obtain_recal_table:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        ref_idx=f"{config['ref']['genome']}.fai",
        ref_dict=os.path.splitext(config["ref"]["genome"])[0] + ".dict",
        known=config["ref"]["known_variants"],
        known_idx=f"{config['ref']['known_variants']}.tbi"
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
        "0.79.0/bio/gatk/baserecalibrator"

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
        "0.79.0/bio/gatk/applybqsr"

rule samtools_index:
    input:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem_mb = get_resource("samtools_index","mem"),
        walltime = get_resource("samtools_index","walltime")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    wrapper:
        "v1.21.6/bio/samtools/index"

rule samtools_index_sorted:
    input:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam"
    output:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem_mb = get_resource("samtools_index","mem"),
        walltime = get_resource("samtools_index","walltime")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    wrapper:
        "v1.21.6/bio/samtools/index"

rule index_known_variants:
    input:
        file=config['ref']['known_variants']
    output:
        index=f"{config['ref']['known_variants']}.tbi"
    resources:
        mem_mb = get_resource("index_known_variants","mem"),
        walltime = get_resource("index_known_variants","walltime")
    conda:
        "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/index_known_variants.log"
    shell:"""
        gatk IndexFeatureFile -I {input.file} -O {output.index}
    """

rule use_alt_contigs_file:
    input:
        file=f"{config['ref']['genome_alt']}"
    output:
        file=f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".alt"
    shell:
        "cp {input.file} {output.file}"
