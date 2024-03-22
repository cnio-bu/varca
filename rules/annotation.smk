rule snpeff_download:
    output:
        directory(f"resources/snpeff/{config['ref']['name']}")
    log:
        f"{LOGDIR}/snpeff/download/{config['ref']['name']}.log"
    params:
        reference=f"{config['ref']['name']}"
    resources:
        threads = get_resource("snpeff_download","threads"),
        mem_mb = get_resource("snpeff_download","mem"),
        runtime = get_resource("snpeff_download","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/snpeff_download.txt"
    wrapper:
        "v3.5.0/bio/snpeff/download"

rule snpeff:
    input:
        calls = f"{OUTDIR}/filtered/{{group}}.vcf.gz",
        db = f"resources/snpeff/{config['ref']['name']}"
    output:
        calls=report(f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz",
        caption="../report/vcf.rst", category="Calls"),
        csvstats=f"{OUTDIR}/snpeff/{{group}}.csv",
        stats="snpeff/{{group}}.html"
    log:
        f"{LOGDIR}/snpeff/{{group}}.snpeff.log"
    resources:
        mem_mb = get_resource("snpeff","mem"),
        runtime = get_resource("snpeff","walltime")
    params:
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("snpeff","threads")),
        extra=""
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.snpeff.txt"
    wrapper:
        "v3.5.0/bio/snpeff/annotate"

rule get_vep_cache:
    output:
        directory(f"{config['annotation']['vep']['cache_directory']}/cache")
    params:
        species=f"{config['annotation']['vep']['species']}",
        build=f"{config['annotation']['vep']['assembly']}",
        release=f"{config['annotation']['vep']['cache_version']}"
    log:
        f"{LOGDIR}/vep/get_vep_cache.log"
    resources:
        threads = get_resource("get_vep_cache","threads"),
        mem_mb = get_resource("get_vep_cache","mem"),
        runtime = get_resource("get_vep_cache","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/get_vep_cache.tsv"
    wrapper:
        "v3.5.0/bio/vep/cache"

rule download_vep_plugins:
    output:
        directory(f"{config['annotation']['vep']['cache_directory']}/plugins")
    params:
        release=f"{config['annotation']['vep']['cache_version']}"
    log:
        f"{LOGDIR}/vep/download_vep_plugins.log"
    resources:
        threads = get_resource("download_vep_plugins","threads"),
        mem_mb = get_resource("download_vep_plugins","mem"),
        runtime = get_resource("download_vep_plugins","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/download_vep_plugin.tsv"
    wrapper:
        "v3.5.0/bio/vep/plugins"

rule vep_gatk:
    input:
        calls=f"{OUTDIR}/filtered/{{group}}.vcf.gz",
        cache=f"{config['annotation']['vep']['cache_directory']}/cache",
        plugins=f"{config['annotation']['vep']['cache_directory']}/plugins"
    output:
        calls=f"{OUTDIR}/annotated/{{group}}.vep.vcf.gz",
        stats=f"{OUTDIR}/annotated/{{group}}.vep.vcf.gz_summary.html"
    params:
        plugins=f"{config['annotation']['vep']['plugins']}",
        extra=f"{config['annotation']['vep']['extra']}"
    log:
        f"{LOGDIR}/vep/{{group}}.gatk_vep.log"
    threads: get_resource("vep_gatk","threads")
    resources:
        mem_mb = get_resource("vep_gatk","mem"),
        runtime = get_resource("vep_gatk","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.vep_gatk.txt"
    wrapper:
        "v3.5.0/bio/vep/annotate"

rule vep_mutect:
    input:
        calls=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered.vcf.gz",
        cache=f"{config['annotation']['vep']['cache_directory']}/cache",
        plugins=f"{config['annotation']['vep']['cache_directory']}/plugins"
    output:
        calls=f"{OUTDIR}/annotated/{{sample}}_mutect.vep.vcf.gz",
        stats=f"{OUTDIR}/annotated/{{sample}}_mutect.vep.vcf.gz_summary.html"
    params:
        plugins=f"{config['annotation']['vep']['plugins']}",
        extra=f"{config['annotation']['vep']['extra']}"
    log:
        f"{LOGDIR}/vep/mutect_{{sample}}_vep.log"
    threads: get_resource("vep_mutect","threads")
    resources:
        mem_mb = get_resource("vep_mutect","mem"),
        runtime = get_resource("vep_mutect","walltime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.vep_mutect.txt"
    wrapper:
        "v3.5.0/bio/vep/annotate"
