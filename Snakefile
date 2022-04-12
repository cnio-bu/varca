import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version
min_version("7.3.7")

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

include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        ["{OUTDIR}/annotated/{group}.snpeff.vcf.gz".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()],
        ["{OUTDIR}/annotated/{group}.vep.vcf.gz".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()],
        ["{OUTDIR}/annotated/{sample}_mutect.vep.vcf.gz".format(OUTDIR=OUTDIR,sample=getattr(row, 'sample')) for row in samples.itertuples() if (getattr(row, 'control') != "-")],
        f"{OUTDIR}/qc/multiqc.html",
        [f"{OUTDIR}/plots/{{group}}.depths.svg".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()],
        [f"{OUTDIR}/plots/{{group}}.allele-freqs.svg".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()]


##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
