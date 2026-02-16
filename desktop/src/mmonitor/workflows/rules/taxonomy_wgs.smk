# WGS Taxonomy Pipeline using Centrifuger
# Workflow: filtlong -> Centrifuger classify -> kreport -> parse results -> upload

# ============ Centrifuger Classification ============
rule centrifuger_classify:
    """Run Centrifuger for WGS taxonomic classification"""
    input:
        fastq=get_filtered_reads
    output:
        output=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger" / f"{SAMPLE}_output.tsv"
    params:
        db=get_db_path("centrifuger"),
        min_hitlen=config.get("centrifuger", {}).get("min_hitlen", 22),
        output_dir=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger"
    threads: get_threads("centrifuger")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "centrifuger_classify.log"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        mkdir -p {params.output_dir}

        centrifuger \
            -x {params.db} \
            -u {input.fastq} \
            -t {threads} \
            --min-hitlen {params.min_hitlen} \
            > {output.output} 2> {log}
        """

# ============ Centrifuger Kreport ============
rule centrifuger_kreport:
    """Generate Kraken-style report from Centrifuger output"""
    input:
        output=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger" / f"{SAMPLE}_output.tsv"
    output:
        kreport=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger" / f"{SAMPLE}_kreport.tsv"
    params:
        db=get_db_path("centrifuger")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "centrifuger_kreport.log"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        centrifuger-kreport \
            -x {params.db} \
            {input.output} > {output.kreport} 2> {log}
        """

# ============ Parse Centrifuger Results ============
rule parse_centrifuger_results:
    """Parse Centrifuger output into standardized JSON format"""
    input:
        kreport=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger" / f"{SAMPLE}_kreport.tsv",
        stats=OUTPUT_DIR / SAMPLE / "stats" / "fastq_stats.tsv"
    output:
        json=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger_parsed.json"
    params:
        sample=SAMPLE,
        project=config.get("project_name", ""),
        subproject=config.get("subproject_name", ""),
        date=config.get("sample_date", "")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "parse_centrifuger.log"
    script:
        "../scripts/parse_centrifuger_output.py"

# ============ Upload Centrifuger Results ============
rule upload_centrifuger_results:
    """Upload Centrifuger results to MMonitor server"""
    input:
        results=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger_parsed.json"
    output:
        status=OUTPUT_DIR / SAMPLE / "taxonomy_wgs" / "centrifuger_results.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "upload_centrifuger.log"
    script:
        "../scripts/upload_taxonomy_results.py"
