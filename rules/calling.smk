if "restrict-regions" in config["processing"]:
    rule sort_bed:
        input:
            config["processing"]["restrict-regions"]
        output:
            temp("called/regions.sorted.bed")
        conda:
            "../envs/bedops.yaml"
        shell:
            "sort-bed {input} > {output}"

    rule compose_regions:
        input:
            "called/regions.sorted.bed"
        output:
            f"{OUTDIR}/called/{{contig}}.regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"

rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions=f"{OUTDIR}/called/{{contig}}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=f"{OUTDIR}/called/{{sample}}.{{contig}}.g.vcf.gz"
    log:
        f"{LOGDIR}/gatk/haplotypecaller/{{sample}}.{{contig}}.log"
    threads: get_resource("call_variants","threads")
    resources:
        mem = get_resource("call_variants","mem")
    params:
        extra=get_call_variants_params
    wrapper:
        "0.35.0/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("{OUTDIR}/called/{sample}.{{contig}}.g.vcf.gz", OUTDIR=OUTDIR, sample=samples.index)
    output:
        gvcf=f"{OUTDIR}/called/all.{{contig}}.g.vcf.gz"
    log:
        f"{LOGDIR}/gatk/combinegvcfs.{{contig}}.log"
    threads: get_resource("combine_calls","threads")
    resources:
        mem = get_resource("combine_calls","mem")
    wrapper:
        "0.35.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf=f"{OUTDIR}/called/all.{{contig}}.g.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/genotyped/all.{{contig}}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        f"{LOGDIR}/gatk/genotypegvcfs.{{contig}}.log"
    threads: get_resource("genotype_variants","threads")
    resources:
        mem = get_resource("genotype_variants","mem")
    wrapper:
        "0.35.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcf=expand(f"{OUTDIR}/genotyped/all.{{contig}}.vcf.gz", contig=contigs)
    output:
        vcf=f"{OUTDIR}/genotyped/all.vcf.gz"
    log:
        f"{LOGDIR}/picard/merge-genotyped.log"
    threads: get_resource("merge_variants","threads")
    resources:
        mem = get_resource("merge_variants","mem")
    wrapper:
        "0.35.0/bio/picard/mergevcfs"

rule merge_bams:
    input: get_sample_bams
    output:
        f"{OUTDIR}/merged_bams/{{sample}}.bam"
    threads: get_resource("merge_bams","threads")
    resources:
        mem = get_resource("merge_bams","mem")
    wrapper:
        "0.35.0/bio/samtools/merge"

rule samtools_index_recal:
    input:
        f"{OUTDIR}/recal/{{sample}}.bam"
    output:
        f"{OUTDIR}/recal/{{sample}}.bam.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem = get_resource("samtools_index","mem")
    log:
        f"{LOGDIR}/samtools/index_recal/{{sample}}.log"
    wrapper:
        "0.35.0/bio/samtools/index"

rule samtools_index_merged:
    input:
        f"{OUTDIR}/merged_bams/{{sample}}.bam"
    output:
        f"{OUTDIR}/merged_bams/{{sample}}.bam.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem = get_resource("samtools_index","mem")
    log:
        f"{LOGDIR}/samtools/index_merged/{{sample}}.log"
    wrapper:
        "0.35.0/bio/samtools/index"

rule mutect:
    input:
        bam=lambda wc: get_merged_bam(wc)[0],
        bai=lambda wc: get_merged_bam(wc)[1],
        ref=config["ref"]["genome"]
    output:
        f"{OUTDIR}/mutect/{{sample}}.vcf.gz"
    threads: get_resource("mutect","threads")
    resources:
        mem = get_resource("mutect","mem")
    conda: "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/mutect.{{sample}}.log"
    shell:"""
        gatk Mutect2 -R {input.ref} -I {input.bam} -tumor {wildcards.sample} -O {output}
    """
