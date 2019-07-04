def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf=f"{OUTDIR}/genotyped/all.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/all.{{vartype}}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        f"{LOGDIR}/gatk/selectvariants/{{vartype}}.log"
    threads: get_resource("select_calls","threads")
    resources:
        mem = get_resource("select_calls","mem")
    wrapper:
        "0.35.0/bio/gatk/selectvariants"


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf=f"{OUTDIR}/filtered/all.{{vartype}}.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/all.{{vartype}}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    threads: get_resource("hard_filter_calls","threads")
    resources:
        mem = get_resource("hard_filter_calls","mem")
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{vartype}}.log"
    wrapper:
        "0.35.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf=f"{OUTDIR}/filtered/all.{{vartype}}.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/all.{{vartype}}.recalibrated.vcf.gz")
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        f"{LOGDIR}/gatk/variantrecalibrator/{{vartype}}.log"
    threads: get_resource("recalibrate_calls","threads")
    resources:
        mem = get_resource("recalibrate_calls","mem")
    wrapper:
        "0.35.0/bio/gatk/variantrecalibrator"


rule merge_calls:
    input:
        vcf=expand(f"{OUTDIR}/filtered/all.{{vartype}}.{{filtertype}}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf=f"{OUTDIR}/filtered/all.vcf.gz"
    log:
        f"{LOGDIR}/picard/merge-filtered.log"
    threads: get_resource("merge_calls","threads")
    resources:
        mem = get_resource("merge_calls","mem")
    wrapper:
        "0.35.0/bio/picard/mergevcfs"

rule filter_mutect_calls:
    input:
        vcf=f"{OUTDIR}/mutect/{{sample}}.vcf.gz",
        ref=config["ref"]["genome"]
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlable.vcf.gz"
    conda: "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/mutect_filter.{{sample}}.log"
    shell:"""
        gatk FilterMutectCalls -R {input.ref} -V {input.vcf} -O {output}
    """
