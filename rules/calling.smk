if "restrict_regions" in config["processing"]:
    rule sort_bed:
        input:
            in_file = config["processing"]["restrict_regions"]
        output:
            temp(f"{OUTDIR}/called/regions.sorted.bed")
        params:
            extra = ""
        resources:
            threads = get_resource("sort_bed","threads"),
            mem_mb = get_resource("sort_bed","mem_mb"),
            runtime = get_resource("sort_bed","runtime")
        benchmark:
            f"{LOGDIR}/benchmarks/sort_bed.txt"
        wrapper:
            "v3.5.0/bio/bedtools/sort"

    rule compose_regions:
        input:
            f"{OUTDIR}/called/regions.sorted.bed"
        output:
            f"{OUTDIR}/called/{{contig}}.regions.bed"
        params:
            contig = lambda wc: get_contig_file_name(wc.contig)
        resources:
            threads = get_resource("compose_regions","threads"),
            mem_mb = get_resource("compose_regions","mem_mb"),
            runtime = get_resource("compose_regions","runtime")
        benchmark:
            f"{LOGDIR}/benchmarks/{{contig}}.compose_regions.txt"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {params.contig} {input} > {output}"

rule call_variants:
    input:
        bam=lambda wc: get_sample_bams(wc.sample),
        ref=config["ref"]["genome"],
        known=config["ref"]["known_variants"],
        regions=f"{OUTDIR}/called/{{contig}}.regions.bed" if config["processing"].get("restrict_regions") else []
    output:
        gvcf=f"{OUTDIR}/called/{{sample}}.{{contig}}.g.vcf.gz"
    log:
        f"{LOGDIR}/gatk/haplotypecaller/{{sample}}.{{contig}}.log"
    threads: get_resource("call_variants","threads")
    resources:
        mem_mb = get_resource("call_variants","mem_mb"),
        runtime = get_resource("call_variants","runtime")
    params:
        extra=get_call_variants_params
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.{{contig}}.call_variants.txt"
    wrapper:
        "v3.5.0/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs= lambda wc: get_vcf_in_group(wc)
    output:
        gvcf=f"{OUTDIR}/called/{{group}}.{{contig}}.g.vcf.gz"
    log:
        f"{LOGDIR}/gatk/combinegvcfs.{{group}}.{{contig}}.log"
    params:
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("combine_calls","threads"))
    resources:
        mem_mb = get_resource("combine_calls","mem_mb"),
        runtime = get_resource("combine_calls","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.{{contig}}.combine_calls.txt"
    wrapper:
        "v3.5.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf=f"{OUTDIR}/called/{{group}}.{{contig}}.g.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/genotyped/{{group}}.{{contig}}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("genotype_variants","threads"))
    log:
        f"{LOGDIR}/gatk/genotypegvcfs.{{group}}.{{contig}}.log"
    resources:
        mem_mb = get_resource("genotype_variants","mem_mb"),
        runtime = get_resource("genotype_variants","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.{{contig}}.genotype_variants.txt"
    wrapper:
        "v3.5.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda wc: expand(f"{OUTDIR}/genotyped/{wc.group}.{{contig}}.vcf.gz", contig=get_contigs())
    output:
        vcf=f"{OUTDIR}/genotyped/{{group}}.vcf.gz"
    log:
        f"{LOGDIR}/picard/{{group}}.merge-genotyped.log"
    params:
        extra = "",
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("merge_variants","threads"))
    resources:
        mem_mb = get_resource("merge_variants","mem_mb"),
        runtime = get_resource("merge_variants","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.merge_variants.txt"
    wrapper:
        "v3.5.0/bio/picard/mergevcfs"

rule merge_bams:
    input: lambda wc: get_sample_bams(wc.sample),
    output:
        f"{OUTDIR}/merged_bams/{{sample}}.bam"
    threads: get_resource("merge_bams","threads")
    resources:
        mem_mb = get_resource("merge_bams","mem_mb"),
        runtime = get_resource("merge_bams","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.merge_bams.txt"
    wrapper:
        "v3.5.0/bio/samtools/merge"

rule samtools_index_merged:
    input:
        f"{OUTDIR}/merged_bams/{{sample}}.bam"
    output:
        f"{OUTDIR}/merged_bams/{{sample}}.bai"
    log:
        f"{LOGDIR}/samtools/index_merged/{{sample}}.log"
    threads: get_resource("samtools_index","threads")
    resources:
        mem_mb = get_resource("samtools_index","mem_mb"),
        runtime = get_resource("samtools_index","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.samtools_index_merged.txt"
    wrapper:
        "v3.5.0/bio/samtools/index"

rule mutect:
    input:
        map=lambda wc: get_merged_bam(wc.sample)[0],
        bai=lambda wc: get_merged_bam(wc.sample)[1],
        cbam=lambda wc: get_merged_bam(samples.loc[(wc.sample),"control"])[0],
        cbai=lambda wc: get_merged_bam(samples.loc[(wc.sample),"control"])[1],
        fasta=config["ref"]["genome"],
    output:
        vcf=f"{OUTDIR}/mutect/{{sample}}.vcf.gz",
        f1r2=f"{OUTDIR}/mutect/{{sample}}.f1r2.tar.gz"
    params:
        extra=lambda wc: get_mutect_params(wc.sample) + get_regions_param()
    log:
        f"{LOGDIR}/gatk/mutect.{{sample}}.log"
    threads: get_resource("mutect","threads")
    resources:
        mem_mb = get_resource("mutect","mem_mb"),
        runtime = get_resource("mutect","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.mutect.txt"
    wrapper:
        "v3.5.0/bio/gatk/mutect"
