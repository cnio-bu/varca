if "restrict-regions" in config["processing"]:
    rule sort_bed:
        input:
            config["processing"]["restrict-regions"]
        output:
            temp("{OUTDIR}/called/regions.sorted.bed")
        conda:
            "../envs/bedops.yaml"
        shell:
            "sort-bed {input} > {output}"

    rule compose_regions:
        input:
            f"{OUTDIR}/called/regions.sorted.bed"
        output:
            f"{OUTDIR}/called/{{contig}}.regions.bed"
        params:
            contig = lambda wc: get_contig_file_name(wc.contig)
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {params.contig} {input} > {output}"

rule call_variants:
    input:
        bam=lambda wc: get_sample_bams(wc.sample),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions=f"{OUTDIR}/called/{{contig}}.regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=f"{OUTDIR}/called/{{sample}}.{{contig}}.g.vcf.gz"
    log:
        f"{LOGDIR}/gatk/haplotypecaller/{{sample}}.{{contig}}.log"
    threads: get_resource("call_variants","threads")
    resources:
        mem_mb = get_resource("call_variants","mem"),
        walltime = get_resource("call_variants","walltime")
    params:
        extra=get_call_variants_params
    wrapper:
        "0.79.0/bio/gatk/haplotypecaller"


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
        mem_mb = get_resource("combine_calls","mem"),
        walltime = get_resource("combine_calls","walltime")
    wrapper:
        "0.79.0/bio/gatk/combinegvcfs"


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
        mem_mb = get_resource("genotype_variants","mem"),
        walltime = get_resource("genotype_variants","walltime")
    wrapper:
        "0.79.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(f"{OUTDIR}/genotyped/all.{{contig}}.vcf.gz", contig=get_contigs())
    output:
        vcf=f"{OUTDIR}/genotyped/all.vcf.gz"
    log:
        f"{LOGDIR}/picard/merge-genotyped.log"
    threads: get_resource("merge_variants","threads")
    resources:
        mem_mb = get_resource("merge_variants","mem"),
        walltime = get_resource("merge_variants","walltime")
    params:
        extra = ""
    wrapper:
        "0.79.0/bio/picard/mergevcfs"

rule merge_bams:
    input: lambda wc: get_sample_bams(wc.sample),
    output:
        f"{OUTDIR}/merged_bams/{{sample}}.bam"
    threads: get_resource("merge_bams","threads")
    resources:
        mem = get_resource("merge_bams","mem"),
        walltime = get_resource("merge_bams","walltime")
    wrapper:
        "0.79.0/bio/samtools/merge"

rule samtools_index_merged:
    input:
        f"{OUTDIR}/merged_bams/{{sample}}.bam"
    output:
        f"{OUTDIR}/merged_bams/{{sample}}.bai"
    threads: get_resource("samtools_index","threads")
    resources:
        mem = get_resource("samtools_index","mem"),
        walltime = get_resource("samtools_index","walltime")
    log:
        f"{LOGDIR}/samtools/index_merged/{{sample}}.log"
    wrapper:
        "0.79.0/bio/samtools/index"

rule mutect:
    input:
        bam=lambda wc: get_merged_bam(wc.sample)[0],
        bai=lambda wc: get_merged_bam(wc.sample)[1],
        cbam=lambda wc: get_merged_bam(samples.loc[(wc.sample),"control"])[0],
        cbai=lambda wc: get_merged_bam(samples.loc[(wc.sample),"control"])[1],
        ref=config["ref"]["genome"]
    output:
        f"{OUTDIR}/mutect/{{sample}}.vcf.gz"
    params:
        regions=lambda wc: get_mutect_params(wc.sample)[0],
        normal=lambda wc: get_mutect_params(wc.sample)[1]
    threads: get_resource("mutect","threads")
    resources:
        mem = get_resource("mutect","mem"),
        walltime = get_resource("mutect","walltime")
    conda: "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/mutect.{{sample}}.log"
    shell:"""
        gatk Mutect2 --callable-depth 1 --native-pair-hmm-threads 16 -R {input.ref} {params.regions} -I {input.bam} {params.normal} -O {output}
    """
