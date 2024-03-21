rule fastqc:
    input:
        lambda wc: get_fastq(wc)[wc.read]
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}-{{read}}_fastqc.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}-{{read}}_fastqc.zip"
    threads: get_resource("fastqc","threads")
    resources:
        mem_mb = get_resource("fastqc","mem"),
        runtime = get_resource("fastqc","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}-{{read}}.fastqc.txt"
    wrapper:
        "v3.5.0/bio/fastqc"

rule samtools_stats:
    input:
        f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/qc/samtools-stats/{{sample}}-{{unit}}.txt"
    log:
        f"{LOGDIR}/samtools-stats/{{sample}}-{{unit}}.log"
    threads: get_resource("samtools_stats","threads")
    resources:
        mem_mb = get_resource("samtools_stats","mem"),
        runtime = get_resource("samtools_stats","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.samtools_stats.txt"
    wrapper:
        "v3.5.0/bio/samtools/stats"

rule genome_dict:
    input:
        genome=config["ref"]["genome"],
        fai=f"{config['ref']['genome']}.fai"
    output:
        dict=re.sub("\.fa","",os.path.splitext(config["ref"]["genome"])[0]) + ".dict"
    log:
        f"{LOGDIR}/picard/genome_dict/CreateSequenceDictionary.log"
    threads: get_resource("genome_dict","threads")
    resources:
        mem_mb = get_resource("genome_dict","mem"),
        runtime = get_resource("genome_dict","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/genome_dict.txt"
    wrapper:
        "v3.5.0/bio/picard/createsequencedictionary"

if "restrict_regions" in config["processing"]:
    rule bed_to_interval:
        input:
            bed=config["processing"]["restrict_regions"],
            dict=re.sub("\.fa.*","",os.path.splitext(config["ref"]["genome"])[0]) + ".dict"
        output:
            temp(f"{OUTDIR}/regions.intervals")
        log:
            f"{LOGDIR}/picard/bed_to_interval/bedtointervals.log"
        params:
            extra = ""
        resources:
            mem_mb =  get_resource("bed_to_interval","mem"),
            runtime = get_resource("bed_to_interval","walltime")
        benchmark:
            f"{LOGDIR}/benchmarks/bed_to_interval.txt"
        wrapper:
            "v3.5.0/bio/picard/bedtointervallist"

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
        params:
            extra = ""
        threads: get_resource("picard_collect_hs_metrics","threads")
        resources:
            mem_mb = get_resource("picard_collect_hs_metrics","mem"),
            runtime = get_resource("picard_collect_hs_metrics","walltime")
        benchmark:
            f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.picard_collect_hs_metrics.txt"
        wrapper:
            "v3.5.0/bio/picard/collecthsmetrics"

rule multiqc:
    input:
         [expand(f"{OUTDIR}/qc/fastqc/{row.sample}-{row.unit}-{{r}}_fastqc.zip", r=["r1"]) for row in units.itertuples() if (str(getattr(row, 'fq2')) == "nan")],
         [expand(f"{OUTDIR}/qc/fastqc/{row.sample}-{row.unit}-{{r}}_fastqc.zip", r=["r1","r2"]) for row in units.itertuples() if (str(getattr(row, 'fq2')) != "nan")],
         expand(f"{OUTDIR}/qc/samtools-stats/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/dedup/{{u.sample}}-{{u.unit}}.metrics.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/picard/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()) if config["processing"].get("restrict_regions") else [],
         expand(f"{OUTDIR}/snpeff/{{u.group}}.csv", u=samples.itertuples())
    output:
        report(f"{OUTDIR}/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        f"{LOGDIR}/multiqc.log"
    threads: get_resource("multiqc","threads")
    resources:
        mem_mb = get_resource("multiqc","mem"),
        runtime = get_resource("multiqc","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/multiqc.txt"
    wrapper:
        "v3.5.0/bio/multiqc"
