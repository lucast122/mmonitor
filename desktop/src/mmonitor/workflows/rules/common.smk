# Common rules shared across pipelines

# ============ Quality Filtering ============
rule filtlong:
    """Filter reads by quality and length using Filtlong"""
    input:
        fastq=get_input_fastq
    output:
        filtered=OUTPUT_DIR / SAMPLE / "filtered" / "filtered.fastq.gz"
    params:
        min_length=config.get("filtlong", {}).get("min_length", 1000),
        min_mean_q=config.get("filtlong", {}).get("min_mean_q", 7.0),
        keep_percent=config.get("filtlong", {}).get("keep_percent", 90.0),
        target_bases=lambda wildcards: config.get("filtlong", {}).get("target_bases") or "",
        max_length=lambda wildcards: config.get("filtlong", {}).get("max_length") or ""
    threads: 4
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "filtlong.log"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        filtlong \
            --min_length {params.min_length} \
            --min_mean_q {params.min_mean_q} \
            --keep_percent {params.keep_percent} \
            $([ -n "{params.target_bases}" ] && echo "--target_bases {params.target_bases}") \
            $([ -n "{params.max_length}" ] && echo "--max_length {params.max_length}") \
            {input.fastq} 2> {log} | gzip > {output.filtered}
        """

rule skip_filtlong:
    """Skip filtering if disabled - just symlink/copy input"""
    input:
        fastq=get_input_fastq
    output:
        filtered=OUTPUT_DIR / SAMPLE / "filtered" / "unfiltered.fastq.gz"
    shell:
        """
        if [[ "{input.fastq}" == *.gz ]]; then
            ln -sf $(realpath {input.fastq}) {output.filtered}
        else
            gzip -c {input.fastq} > {output.filtered}
        fi
        """

def get_filtered_reads(wildcards=None):
    """Get filtered or unfiltered reads based on config"""
    if config.get("filtlong", {}).get("enabled", True):
        return OUTPUT_DIR / SAMPLE / "filtered" / "filtered.fastq.gz"
    else:
        return OUTPUT_DIR / SAMPLE / "filtered" / "unfiltered.fastq.gz"

# ============ Statistics ============
rule fastq_stats:
    """Calculate FASTQ statistics using seqkit"""
    input:
        fastq=get_filtered_reads
    output:
        stats=OUTPUT_DIR / SAMPLE / "stats" / "fastq_stats.tsv"
    threads: 4
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "fastq_stats.log"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        seqkit stats -T -a {input.fastq} > {output.stats} 2> {log}
        """

# ============ Result Upload ============
rule upload_results:
    """Upload results to MMonitor server"""
    input:
        results="{results_file}"
    output:
        upload_status="{results_file}.uploaded"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        sample=SAMPLE,
        project=config.get("project_name", ""),
        subproject=config.get("subproject_name", ""),
        date=config.get("sample_date", "")
    log:
        "{results_file}.upload.log"
    script:
        "../scripts/upload_results.py"
