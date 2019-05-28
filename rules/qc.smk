rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}.zip"
    threads: get_resource("fastqc","threads")
    resources:
        mem = get_resource("fastqc","mem")
    wrapper:
        "0.27.1/bio/fastqc"


rule samtools_stats:
    input:
        f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/qc/samtools-stats/{{sample}}-{{unit}}.txt"
    log:
        f"{LOGDIR}/samtools-stats/{{sample}}-{{unit}}.log"
    threads: get_resource("samtools_stats","threads")
    resources:
        mem = get_resource("samtools_stats","mem")
    wrapper:
        "0.27.1/bio/samtools/stats"


rule multiqc:
    input:
         expand(f"{OUTDIR}/qc/fastqc/{{u.sample}}-{{u.unit}}.zip", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/samtools-stats/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/dedup/{{u.sample}}-{{u.unit}}.metrics.txt", u=units.itertuples()),
         f"{OUTDIR}/snpeff/all.csv"
    output:
        report(f"{OUTDIR}/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        f"{LOGDIR}/multiqc.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem = get_resource("multiqc","mem")
    wrapper:
        "0.27.1/bio/multiqc"
