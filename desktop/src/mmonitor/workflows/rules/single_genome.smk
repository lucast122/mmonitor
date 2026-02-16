# Single Genome Annotation Pipeline
# Workflow: Input genome -> Bakta annotation -> Create index -> Upload to server
#
# This pipeline is for annotating a single genome (not from binning)
# Useful for testing and for annotating reference genomes

import os
from pathlib import Path

# Get genome ID from config
GENOME_ID = config.get("genome_id", "input_genome")

# ============ Copy Input Genome ============
rule copy_input_genome:
    """Copy input genome FASTA to working directory"""
    input:
        fasta=config.get("input_fasta", [])
    output:
        fasta=OUTPUT_DIR / SAMPLE / "genome" / f"{GENOME_ID}.fasta"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        cp {input.fasta} {output.fasta}
        """

# ============ Bakta Annotation for Single Genome ============
rule bakta_single_genome:
    """Annotate a single genome using Bakta"""
    input:
        fasta=OUTPUT_DIR / SAMPLE / "genome" / f"{GENOME_ID}.fasta",
        db_ready=DB_BASE / "bakta" / ".db_ready"
    output:
        gff=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.gff3",
        faa=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.faa",
        tsv=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.tsv",
        fna=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna"
    params:
        db=get_db_path("bakta"),
        output_dir=str(OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID),
        min_contig=config.get("bakta", {}).get("min_contig_length", 1),
        prefix=GENOME_ID
    threads: get_threads("bakta")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / f"bakta_{GENOME_ID}.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.output_dir}

        bakta \
            --db {params.db} \
            --output {params.output_dir} \
            --prefix {params.prefix} \
            --threads {threads} \
            --min-contig-length {params.min_contig} \
            --force \
            {input.fasta} \
            2>&1 | tee {log}
        """

# ============ Create FASTA Index ============
rule create_fasta_index:
    """Create samtools faidx index for genome browser"""
    input:
        fasta=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna"
    output:
        fai=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input.fasta}
        """

# ============ Prepare Upload Package ============
rule prepare_single_genome_upload:
    """Prepare genome data for upload to MMonitor server"""
    input:
        fasta=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna",
        fai=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna.fai",
        gff=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.gff3",
        tsv=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.tsv"
    output:
        json=OUTPUT_DIR / SAMPLE / "upload" / f"{GENOME_ID}_upload.json"
    params:
        sample=SAMPLE,
        genome_id=GENOME_ID,
        taxonomy=config.get("taxonomy", {}),
        project=config.get("project_name", "")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / f"prepare_upload_{GENOME_ID}.log"
    script:
        "../scripts/prepare_single_genome_upload.py"

# ============ Upload to Server ============
rule upload_single_genome:
    """Upload genome and annotations to MMonitor server"""
    input:
        upload_json=OUTPUT_DIR / SAMPLE / "upload" / f"{GENOME_ID}_upload.json",
        fasta=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna",
        fai=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.fna.fai",
        gff=OUTPUT_DIR / SAMPLE / "annotation" / GENOME_ID / f"{GENOME_ID}.gff3"
    output:
        status=OUTPUT_DIR / SAMPLE / "upload" / f"{GENOME_ID}_complete.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTPUT_DIR / SAMPLE / "logs" / f"upload_{GENOME_ID}.log"
    script:
        "../scripts/upload_single_genome.py"

# ============ Complete Pipeline Target ============
rule annotate_genome_complete:
    """Complete single genome annotation and upload pipeline"""
    input:
        OUTPUT_DIR / SAMPLE / "upload" / f"{GENOME_ID}_complete.json"
    output:
        OUTPUT_DIR / SAMPLE / "annotation_complete.txt"
    shell:
        """
        echo "Single genome annotation complete for {GENOME_ID}" > {output}
        echo "Timestamp: $(date)" >> {output}
        """
