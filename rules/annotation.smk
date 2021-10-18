rule snpeff_download:
    output:
        directory(f"resources/snpeff/{config['ref']['name']}")
    log:
        f"{LOGDIR}/snpeff/download/{config['ref']['name']}.log"
    params:
        reference=f"{config['ref']['name']}"
    resources:
        mem_mb = get_resource("snpeff","mem"),
        walltime = get_resource("snpeff","walltime")
    wrapper:
        "0.79.0/bio/snpeff/download"

rule snpeff:
    input:
        calls = f"{OUTDIR}/filtered/all.vcf.gz",
        db = f"resources/snpeff/{config['ref']['name']}"
    output:
        calls=report(f"{OUTDIR}/annotated/all.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats=f"{OUTDIR}/snpeff/all.csv"
    log:
        f"{LOGDIR}/snpeff/snpeff.log"
    threads: get_resource("snpeff","threads")
    resources:
        mem_mb = get_resource("snpeff","mem"),
        walltime = get_resource("snpeff","walltime")
    params:
        extra=""
    wrapper:
        "0.79.0/bio/snpeff/annotate"

rule vep_gatk:
    input:
        f"{OUTDIR}/filtered/all.vcf.gz"
    output:
        f"{OUTDIR}/annotated/all.vep.vcf.gz"
    params:
        get_vep_params()
    conda:"../envs/vep.yaml"
    threads: get_resource("vep","threads")
    resources:
        mem = get_resource("vep","mem"),
        walltime = get_resource("vep","walltime")
    log:
        f"{LOGDIR}/vep/gatk_vep.log"
    shell: """
        vep -i {input} -o {output} --vcf --compress_output gzip --force_overwrite {params} --fork {threads}
    """

rule vep_mutect:
    input:
        f"{OUTDIR}/mutect_filter/{{sample}}_passlable_filtered.vcf.gz"
    output:
        f"{OUTDIR}/annotated/{{sample}}_mutect.vep.vcf.gz"
    params:
        get_vep_params()
    conda:"../envs/vep.yaml"
    threads: get_resource("vep","threads")
    resources:
        mem = get_resource("vep","mem"),
        walltime = get_resource("vep","walltime")
    log:
        f"{LOGDIR}/vep/mutect_{{sample}}_vep.log"
    shell: """
        vep -i {input} -o {output} --vcf --compress_output gzip --force_overwrite {params} --fork {threads}
    """
