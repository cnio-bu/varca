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
        "https://bitbucket.org/tdido/snakemake-wrappers/raw/195d5bbcbbdab720554bc4046c5015e1c6e04114/bio/fastqc"


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
        "0.35.0/bio/samtools/stats"

if "restrict-regions" in config["processing"]:
    rule bed_to_interval:
        input:
            file=config["processing"]["restrict-regions"],
            SD=config["ref"]["genome"]
        output:
            temp("{OUTDIR}/regions.intervals")
        conda:
            "../envs/picard.yaml"
        shell:
            "picard BedToIntervalList I={input.file} O={output} SD={input.SD}"

    rule picard_collect_hs_metrics:
        input:
            bam=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam",
            reference=config["ref"]["genome"],
            bait_intervals=f"{OUTDIR}/regions.intervals",
            target_intervals=f"{OUTDIR}/regions.intervals"
        output:
            f"{OUTDIR}/qc/picard/{{sample}}-{{unit}}.txt"
        log:
            f"{LOGDIR}/picard_collect_hs_metrics/{{sample}}-{{unit}}.log"
        threads: 1
        resources:
            mem = 4000
        wrapper:
            "0.35.0/bio/picard/collecthsmetrics"

rule multiqc:
    input:
         expand(f"{OUTDIR}/qc/fastqc/{{u.sample}}-{{u.unit}}.zip", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/samtools-stats/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/dedup/{{u.sample}}-{{u.unit}}.metrics.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/picard/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()) if config["processing"].get("restrict-regions") else [],
         f"{OUTDIR}/snpeff/all.csv"
    output:
        report(f"{OUTDIR}/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        f"{LOGDIR}/multiqc.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem = get_resource("multiqc","mem")
    wrapper:
        "0.35.0/bio/multiqc"
