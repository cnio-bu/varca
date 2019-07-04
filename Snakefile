import pandas as pd
from snakemake.utils import validate

singularity: "docker://continuumio/miniconda3:4.6.14"

report: "report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

samples = pd.read_csv(config["samples"],sep="\t").set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"],sep="\t", dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

# contigs in reference genome
contigs = pd.read_csv(config["ref"]["genome"] + ".fai", sep="\t",
                        header=None, usecols=[0], squeeze=True, dtype=str)

include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        f"{OUTDIR}/annotated/all.vcf.gz",
        ["{OUTDIR}/mutect/{sample}.vcf.gz".format(OUTDIR=OUTDIR,sample=s[1]) for s in samples.itertuples() if (s[2] != "-")],
        ["{OUTDIR}/mutect_filter/{sample}_passlable.vcf.gz".format(OUTDIR=OUTDIR,sample=s[1]) for s in samples.itertuples() if (s[2] != "-")],
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
