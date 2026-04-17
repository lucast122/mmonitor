# WGS Taxonomy Pipeline using Centrifuger
# Workflow: filtlong -> Centrifuger classify -> kreport -> parse results -> upload
# All rules use {sample} wildcard for multi-sample support
#
# When multiple samples are configured, centrifuger is run once with --sample-sheet
# to avoid reloading the index for each sample.

# ============ Centrifuger Classification (batch) ============
def _all_filtered_reads_for_centrifuger(wildcards):
    """Get filtered reads for all samples (used by the batch rule)."""
    inputs = []
    for sample in SAMPLES:
        if config.get("filtlong", {}).get("enabled", True):
            inputs.append(OUTDIR + f"/{sample}/filtered/filtered.fastq.gz")
        else:
            inputs.append(OUTDIR + f"/{sample}/filtered/unfiltered.fastq.gz")
    return inputs


rule centrifuger_classify_batch:
    """Run Centrifuger for all samples in one invocation using --sample-sheet.

    This loads the index only once, which is critical for large databases
    like GTDB (~176 GB). The sample sheet tells centrifuger to write each
    sample's output to a separate file.
    """
    input:
        fastqs=_all_filtered_reads_for_centrifuger
    output:
        outputs=expand(
            OUTDIR + "/{sample}/taxonomy_wgs/centrifuger/{sample}_output.tsv",
            sample=SAMPLES.keys()
        )
    params:
        db=get_db_path("centrifuger"),
        min_hitlen=config.get("centrifuger", {}).get("min_hitlen", 22),
        sample_names=list(SAMPLES.keys()),
    threads: get_threads("centrifuger")
    log:
        OUTDIR + "/logs/centrifuger_classify_batch.log"
    conda:
        "../envs/taxonomy.yaml"
    run:
        import tempfile
        sample_names = params.sample_names
        fastqs = input.fastqs
        outputs = output.outputs

        # Ensure output directories exist
        for out_path in outputs:
            Path(out_path).parent.mkdir(parents=True, exist_ok=True)

        if len(sample_names) == 1:
            # Single sample: run directly without sample sheet
            shell(
                "centrifuger "
                "-x {params.db} "
                "-u {input.fastqs[0]} "
                "-t {threads} "
                "--min-hitlen {params.min_hitlen} "
                "> {output.outputs[0]} 2> {log}"
            )
        else:
            # Multiple samples: use --sample-sheet to load index only once
            # Format: read1 read2 barcode UMI output
            # For single-end nanopore: read1 . . . output
            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.tsv', delete=False, prefix='centrifuger_samples_'
            ) as f:
                sheet_path = f.name
                for fastq, out_path in zip(fastqs, outputs):
                    f.write(f"{fastq}\t.\t.\t.\t{out_path}\n")

            try:
                shell(
                    "centrifuger "
                    "-x {params.db} "
                    "--sample-sheet " + sheet_path + " "
                    "-t {threads} "
                    "--min-hitlen {params.min_hitlen} "
                    "2> {log}"
                )
            finally:
                import os
                os.unlink(sheet_path)


# ============ Centrifuger Kreport ============
rule centrifuger_kreport:
    """Generate Kraken-style report from Centrifuger output"""
    input:
        output=OUTDIR + "/{sample}/taxonomy_wgs/centrifuger/{sample}_output.tsv"
    output:
        kreport=OUTDIR + "/{sample}/taxonomy_wgs/centrifuger/{sample}_kreport.tsv"
    params:
        db=get_db_path("centrifuger")
    log:
        OUTDIR + "/{sample}/logs/centrifuger_kreport.log"
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
        kreport=OUTDIR + "/{sample}/taxonomy_wgs/centrifuger/{sample}_kreport.tsv",
        stats=OUTDIR + "/{sample}/stats/fastq_stats.tsv"
    output:
        json=OUTDIR + "/{sample}/taxonomy_wgs/centrifuger_parsed.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        project=config.get("project_name", ""),
        subproject=config.get("subproject_name", ""),
        date=config.get("sample_date", "")
    log:
        OUTDIR + "/{sample}/logs/parse_centrifuger.log"
    script:
        "../scripts/parse_centrifuger_output.py"

# ============ Upload Centrifuger Results ============
rule upload_centrifuger_results:
    """Upload Centrifuger results to MMonitor server"""
    input:
        results=OUTDIR + "/{sample}/taxonomy_wgs/centrifuger_parsed.json"
    output:
        status=OUTDIR + "/{sample}/taxonomy_wgs/centrifuger_results.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTDIR + "/{sample}/logs/upload_centrifuger.log"
    script:
        "../scripts/upload_taxonomy_results.py"
