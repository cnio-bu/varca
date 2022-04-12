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
        calls = f"{OUTDIR}/filtered/{{group}}.vcf.gz",
        db = f"resources/snpeff/{config['ref']['name']}"
    output:
        calls=report(f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats=f"{OUTDIR}/snpeff/{{group}}.csv"
    log:
        f"{LOGDIR}/snpeff/{{group}}.snpeff.log"
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
        f"{OUTDIR}/filtered/{{group}}.vcf.gz"
    output:
        f"{OUTDIR}/annotated/{{group}}.vep.vcf.gz"
    params:
        get_vep_params()
    conda:"../envs/vep.yaml"
    threads: get_resource("vep","threads")
    resources:
        mem_mb = get_resource("vep","mem"),
        walltime = get_resource("vep","walltime")
    log:
        f"{LOGDIR}/vep/{{group}}.gatk_vep.log"
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
        mem_mb = get_resource("vep","mem"),
        walltime = get_resource("vep","walltime")
    log:
        f"{LOGDIR}/vep/mutect_{{sample}}_vep.log"
    shell: """
        vep -i {input} -o {output} --format 'vcf' --vcf --compress_output gzip --force_overwrite {params} --fork {threads}
    """
