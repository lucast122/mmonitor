# 16S rRNA Taxonomy Pipeline using EMU
# Workflow: filtlong -> EMU abundance -> parse results -> upload
# All rules use {sample} wildcard for multi-sample support

# ============ EMU Abundance Estimation ============
rule emu_abundance:
    """Run EMU for 16S rRNA abundance estimation"""
    input:
        fastq=get_filtered_reads,
        db_ready=DB_BASE / "emu" / ".db_ready"
    output:
        abundance=OUTDIR + "/{sample}/taxonomy_16s/emu/{sample}_rel-abundance.tsv",
        alignments=OUTDIR + "/{sample}/taxonomy_16s/emu/{sample}_emu_alignments.sam"
    params:
        db=get_db_path("emu"),
        min_abundance=config.get("emu", {}).get("min_abundance", 0.0001),
        output_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/taxonomy_16s/emu",
        output_basename=lambda wildcards: wildcards.sample
    threads: get_threads("emu")
    log:
        OUTDIR + "/{sample}/logs/emu_abundance.log"
    conda:
        "../envs/taxonomy.yaml"
    shell:
        """
        mkdir -p {params.output_dir}

        emu abundance \
            --db {params.db} \
            --threads {threads} \
            --output-dir {params.output_dir} \
            --output-basename {params.output_basename} \
            --min-abundance {params.min_abundance} \
            --keep-files \
            --N 50 \
            {input.fastq} 2>&1 | tee {log}
        """

# ============ Parse EMU Results ============
rule parse_emu_results:
    """Parse EMU output into standardized JSON format"""
    input:
        abundance=OUTDIR + "/{sample}/taxonomy_16s/emu/{sample}_rel-abundance.tsv",
        stats=OUTDIR + "/{sample}/stats/fastq_stats.tsv"
    output:
        json=OUTDIR + "/{sample}/taxonomy_16s/emu_parsed.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        project=config.get("project_name", ""),
        subproject=config.get("subproject_name", ""),
        date=config.get("sample_date", ""),
        min_abundance=config.get("emu", {}).get("min_abundance", 0.0001),
        emu_db=get_db_path("emu")
    log:
        OUTDIR + "/{sample}/logs/parse_emu.log"
    script:
        "../scripts/parse_emu_output.py"

# ============ Upload EMU Results ============
rule upload_emu_results:
    """Upload EMU results to MMonitor server"""
    input:
        results=OUTDIR + "/{sample}/taxonomy_16s/emu_parsed.json"
    output:
        status=OUTDIR + "/{sample}/taxonomy_16s/emu_results.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTDIR + "/{sample}/logs/upload_emu.log"
    script:
        "../scripts/upload_taxonomy_results.py"
