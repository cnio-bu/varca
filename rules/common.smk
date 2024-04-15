##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),
    group="|".join([str(x) for x in samples["group"]])

##### Helper functions #####
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        if os.stat(config["contigs"]).st_size != 0:
            contigs = pd.read_csv(config["contigs"],sep="\t",header=None,usecols=[0],dtype=str).squeeze("columns")
        else:
            contigs = pd.read_table(fai, header=None, usecols=[0], dtype=str).squeeze("columns")
        return [re.sub(r"\*","___",x) for x in contigs]

def get_resource(rule,resource):
    try:
        return config["resources"][rule][resource]
    except KeyError:
        return config["resources"]["default"][resource]

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("{OUTDIR}/trimmed/{sample}-{unit}.{group}.fastq.gz",
                      OUTDIR=OUTDIR, group=[1, 2], **wildcards)
    # single end sample
    return f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.fastq.gz".format(**wildcards)

def get_trimmed_reads_qc(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(wildcards.sample, wildcards.unit):
        # paired-end sample
        return {"r1": "{OUTDIR}/trimmed/{sample}-{unit}.1.fastq.gz".format(OUTDIR=OUTDIR, sample=wildcards.sample, unit=wildcards.unit), "r2": "{OUTDIR}/trimmed/{sample}-{unit}.2.fastq.gz".format(OUTDIR=OUTDIR, sample=wildcards.sample, unit=wildcards.unit)}
    # single end sample
    return {"r1": "{OUTDIR}/trimmed/{sample}-{unit}.fastq.gz".format(OUTDIR=OUTDIR, sample=wildcards.sample, unit=wildcards.unit)}


def get_sample_bams(sample):
    """Get all aligned reads of given sample."""
    return expand(f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam",
                    sample=sample,
                    unit=units.loc[sample].unit)

def get_merged_bam(sample):
    """Merge aligned reads if there are multiple units."""
    bams = get_sample_bams(sample)
    if len(bams) > 1:
        return f"{OUTDIR}/merged_bams/{sample}.bam",f"{OUTDIR}/merged_bams/{sample}.bai"
    else:
        return bams[0],bams[0].replace(".bam",".bai")

def get_regions_param(regions=config["processing"].get("restrict_regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region_padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default=f"--intervals {wildcards.contig}") +
            config["params"]["gatk"]["HaplotypeCaller"])


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam"
    if config["processing"]["remove_duplicates"]:
        # case 2: remove duplicates
        f = f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"
    if bai:
        if config["processing"].get("restrict_regions"):
            # case 3: need an index because random access is required
            f = f.replace(".bam",".bai")
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f

def get_mutect_params(sample):
    if samples.loc[(sample),"control"] != sample:
        normal_call="-I "+get_merged_bam(samples.loc[(sample),"control"])[0]+" -normal "+samples.loc[(sample),"control"]+" "
    else:
        normal_call=" "
    return normal_call

def get_contig_file_name(contigfn):
    contig = re.sub("___","*",contigfn)
    return contig

def get_vcf_in_group(wc):
    vcfs = []
    for row in samples.itertuples():
        if str(getattr(row, 'group')) == str(wc.group):
            vcfs = vcfs + ["{OUTDIR}/called/{sample}.{contig}.g.vcf.gz".format(OUTDIR=OUTDIR, sample=getattr(row, 'sample'), contig=wc.contig)]
    return vcfs
