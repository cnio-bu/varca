rule vcf_to_tsv:
    input:
        f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz"
    output:
        report(f"{OUTDIR}/tables/{{group}}.calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    threads: get_resource("vcf_to_tsv","threads")
    resources:
        mem = get_resource("vcf_to_tsv","mem"),
        walltime = get_resource("vcf_to_tsv","walltime")
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


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
        mem = get_resource("plots_stats","mem"),
        walltime = get_resource("plots_stats","walltime")
    script:
        "../scripts/plot-depths.py"
