include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        f"{OUTDIR}/annotated/all.vcf.gz",
        f"{OUTDIR}/qc/multiqc.html",
        f"{OUTDIR}/plots/depths.svg",
        f"{OUTDIR}/plots/allele-freqs.svg"


##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
