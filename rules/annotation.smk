rule snpeff:
    input:
        f"{OUTDIR}/filtered/all.vcf.gz",
    output:
        vcf=report(f"{OUTDIR}/annotated/all.vcf.gz", caption="../report/vcf.rst", category="Calls"),
        csvstats=f"{OUTDIR}/snpeff/all.csv"
    log:
        "logs/snpeff.log"
    threads: get_resource("snpeff","threads")
    resources:
        mem = get_resource("snpeff","mem")
    params:
        reference=config["ref"]["name"],
        extra="-Xmx{}m".format(get_resource("snpeff","mem"))
    wrapper:
        "0.35.0/bio/snpeff"
