def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf=f"{OUTDIR}/genotyped/{{group}}.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/{{group}}.{{vartype}}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        f"{LOGDIR}/gatk/selectvariants/{{group}}.{{vartype}}.log"
    threads: get_resource("select_calls","threads")
    resources:
        mem_mb = get_resource("select_calls","mem"),
        walltime = get_resource("select_calls","walltime")
    wrapper:
        "0.79.0/bio/gatk/selectvariants"


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/{{group}}.{{vartype}}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    threads: get_resource("hard_filter_calls","threads")
    resources:
        mem_mb = get_resource("hard_filter_calls","mem"),
        walltime = get_resource("hard_filter_calls","walltime")
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{group}}.{{vartype}}.log"
    wrapper:
        "0.79.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.vcf.gz",
        ref=config["ref"]["genome"],
        hapmap=config["params"]["gatk"]["VariantRecalibrator"]["hapmap"],
        omni=config["params"]["gatk"]["VariantRecalibrator"]["omni"],
        g1k=config["params"]["gatk"]["VariantRecalibrator"]["g1k"],
        dbsnp=config["params"]["gatk"]["VariantRecalibrator"]["dbsnp"],
        aux=config["params"]["gatk"]["VariantRecalibrator"]["aux"]
    output:
        vcf=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.recalibrated.vcf.gz",
        tranches=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.tranches"
    params:
        mode=lambda wc: "SNP"
                if wc.vartype == "snvs"
                else "INDEL",
        resources=config["params"]["gatk"]["VariantRecalibrator"]["parameters"],
        annotation=config["params"]["gatk"]["VariantRecalibrator"]["annotation"],
        extra=config["params"]["gatk"]["VariantRecalibrator"]["extra"]
    log:
        f"{LOGDIR}/gatk/variantrecalibrator/{{group}}.{{vartype}}.log"
    threads: get_resource("recalibrate_calls","threads")
    resources:
        mem_mb = get_resource("recalibrate_calls","mem"),
        walltime = get_resource("recalibrate_calls","walltime")
    wrapper:
        "0.79.0/bio/gatk/variantrecalibrator"

rule merge_calls:
    input:
        vcfs=lambda wc: expand(f"{OUTDIR}/filtered/{wc.group}.{{vartype}}.{{filtertype}}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf=f"{OUTDIR}/filtered/{{group}}.vcf.gz"
    log:
        f"{LOGDIR}/picard/{{group}}.merge-filtered.log"
    threads: get_resource("merge_calls","threads")
    resources:
        mem_mb = get_resource("merge_calls","mem"),
        walltime = get_resource("merge_calls","walltime")
    params:
        extra = ""
    wrapper:
        "0.79.0/bio/picard/mergevcfs"

rule learn_read_orientation_model:
    input:
        f1r2=f"{OUTDIR}/mutect/{{sample}}.f1r2.tar.gz"
    output:
        rom=f"{OUTDIR}/mutect/{{sample}}_read_orientation_model.tar.gz"
    conda: "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/read_orientation_model.{{sample}}.log"
    shell:"""
        gatk LearnReadOrientationModel -I {input.f1r2} -O {output.rom}
    """


rule filter_mutect_calls:
    input:
        vcf=f"{OUTDIR}/mutect/{{sample}}.vcf.gz",
        ref=config["ref"]["genome"],
        rom=f"{OUTDIR}/mutect/{{sample}}_read_orientation_model.tar.gz" if config["filtering"]["mutect"]["lrom"] else []
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlable.vcf.gz"
    params:
        rom="--ob-priors " if config["filtering"]["mutect"]["lrom"] else ""
    conda: "../envs/gatk.yaml"
    log:
        f"{LOGDIR}/gatk/mutect_filter.{{sample}}.log"
    shell:"""
        gatk FilterMutectCalls -R {input.ref} -V {input.vcf} {params.rom} {input.rom} -O {output}
    """

rule filter_mutect_2:
    input:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlable.vcf.gz",
        ref=config["ref"]["genome"]
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlable_filtered.vcf.gz"
    params:
        filters={"DPfilter": config["filtering"]["mutect"]["depth"]}
    threads: get_resource("hard_filter_calls","threads")
    resources:
        mem_mb = get_resource("hard_filter_calls","mem"),
        walltime = get_resource("hard_filter_calls","walltime")
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{sample}}_mutect.log"
    wrapper:
        "0.79.0/bio/gatk/variantfiltration"
