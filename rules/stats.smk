rule index_snpeff_file:
    input:
        f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz"
    output:
        f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz.tbi"
    params:
        extra=""
    resources:
        mem_mb = get_resource("index_snpeff_file","mem"),
        runtime = get_resource("index_snpeff_file","walltime")
    log:
        f"{LOGDIR}/gatk/{{group}}.index_snpeff_file.log"
    wrapper:
        "v3.5.0/bio/bcftools/index"

rule vcf_to_tsv:
    input:
        vcf=f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz",
        tbi=f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz.tbi"
    output:
        tab=temp(f"{OUTDIR}/tables/{{group}}.calls.tsv")
    log:
        f"{LOGDIR}/gatk/variantstotable.{{group}}.snpeff.log"
    params:
        extra="-F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -GF GT -GF DP -GF AD",
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("vcf_to_tsv","threads"))
    resources:
        mem_mb = get_resource("vcf_to_tsv","mem"),
        runtime = get_resource("vcf_to_tsv","walltime")
    wrapper:
        "v3.5.0/bio/gatk/variantstotable"

rule compress_tsv:
    input:
        f"{OUTDIR}/tables/{{group}}.calls.tsv"
    output:
        report(f"{OUTDIR}/tables/{{group}}.calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    log:
        f"{LOGDIR}/bgzip/{{group}}.snpeff.log"
    params:
        extra=""
    threads: get_resource("vcf_to_tsv","threads")
    resources:
        mem_mb = get_resource("vcf_to_tsv","mem"),
        runtime = get_resource("vcf_to_tsv","walltime")
    wrapper:
        "v3.5.0/bio/bgzip"

rule plot_stats:
    input:
        f"{OUTDIR}/tables/{{group}}.calls.tsv.gz"
    output:
        depths=report(f"{OUTDIR}/plots/{{group}}.depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report(f"{OUTDIR}/plots/{{group}}.allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    threads: get_resource("plot_stats","threads")
    resources:
        mem_mb = get_resource("plots_stats","mem"),
        runtime = get_resource("plots_stats","walltime")
    script:
        "../scripts/plot-depths.py"
