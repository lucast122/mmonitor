# Functional Annotation Pipeline
# Uses eggNOG-mapper for COG/KEGG/GO annotations
# Uses InterProScan for PFAM/TIGRFAM domain annotations

# ============ eggNOG-mapper ============
rule eggnog_mapper:
    """Run eggNOG-mapper for COG, KEGG, GO annotations"""
    input:
        proteins=OUTPUT_DIR / SAMPLE / "annotation" / "bakta" / "{bin_id}" / "{bin_id}.faa",
        db_ready=DB_BASE / "eggnog" / ".db_ready"
    output:
        annotations=OUTPUT_DIR / SAMPLE / "functional" / "eggnog" / "{bin_id}.emapper.annotations",
        orthologs=OUTPUT_DIR / SAMPLE / "functional" / "eggnog" / "{bin_id}.emapper.seed_orthologs"
    params:
        db=get_db_path("eggnog"),
        output_prefix=lambda wildcards: str(OUTPUT_DIR / SAMPLE / "functional" / "eggnog" / wildcards.bin_id),
        sensmode=config.get("eggnog", {}).get("sensmode", "default")
    threads: get_threads("eggnog")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "eggnog_{bin_id}.log"
    conda:
        "../envs/functional.yaml"
    shell:
        """
        emapper.py \
            -i {input.proteins} \
            --output {params.output_prefix} \
            --data_dir {params.db} \
            --cpu {threads} \
            --sensmode {params.sensmode} \
            --override \
            2>&1 | tee {log}
        """


# ============ InterProScan ============
rule interproscan:
    """Run InterProScan for detailed PFAM, TIGRFAM, CDD domain annotations"""
    input:
        proteins=OUTPUT_DIR / SAMPLE / "annotation" / "bakta" / "{bin_id}" / "{bin_id}.faa"
    output:
        tsv=OUTPUT_DIR / SAMPLE / "functional" / "interpro" / "{bin_id}.tsv",
        gff=OUTPUT_DIR / SAMPLE / "functional" / "interpro" / "{bin_id}.gff3"
    params:
        output_dir=OUTPUT_DIR / SAMPLE / "functional" / "interpro",
        applications=config.get("interproscan", {}).get("applications", "Pfam,TIGRFAM,CDD"),
        interproscan_path=config.get("interproscan", {}).get("path", "interproscan.sh")
    threads: get_threads("interproscan")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "interproscan_{bin_id}.log"
    shell:
        """
        mkdir -p {params.output_dir}

        {params.interproscan_path} \
            -i {input.proteins} \
            -d {params.output_dir} \
            -f TSV,GFF3 \
            -b {wildcards.bin_id} \
            -appl {params.applications} \
            --cpu {threads} \
            --goterms \
            --pathways \
            2>&1 | tee {log}
        """


# ============ Aggregate Functional Annotations ============
def get_all_eggnog_outputs(wildcards):
    """Get all eggNOG output files for HQ bins"""
    checkpoint_output = checkpoints.filter_hq_bins.get(**wildcards).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTPUT_DIR / SAMPLE / "functional" / "eggnog" / "{bin_id}.emapper.annotations",
        bin_id=bin_ids
    )

def get_all_interpro_outputs(wildcards):
    """Get all InterProScan output files for HQ bins"""
    checkpoint_output = checkpoints.filter_hq_bins.get(**wildcards).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTPUT_DIR / SAMPLE / "functional" / "interpro" / "{bin_id}.tsv",
        bin_id=bin_ids
    )

def get_all_bakta_gff(wildcards):
    """Get all Bakta GFF files for HQ bins"""
    checkpoint_output = checkpoints.filter_hq_bins.get(**wildcards).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTPUT_DIR / SAMPLE / "annotation" / "bakta" / "{bin_id}" / "{bin_id}.gff3",
        bin_id=bin_ids
    )


rule aggregate_functional:
    """Aggregate eggNOG + InterProScan results into summaries"""
    input:
        eggnog_files=get_all_eggnog_outputs,
        interpro_files=get_all_interpro_outputs,
        bakta_gff_files=get_all_bakta_gff,
        gtdbtk_summary=OUTPUT_DIR / SAMPLE / "annotation" / "gtdbtk" / "gtdbtk.bac120.summary.tsv",
        checkm_quality=OUTPUT_DIR / SAMPLE / "binning" / "checkm2" / "quality_report.tsv"
    output:
        summary=OUTPUT_DIR / SAMPLE / "functional" / "summary" / "functional_summary.json",
        annotations=OUTPUT_DIR / SAMPLE / "functional" / "summary" / "all_annotations.tsv"
    params:
        sample=SAMPLE,
        project=config.get("project_name", ""),
        eggnog_dir=OUTPUT_DIR / SAMPLE / "functional" / "eggnog",
        interpro_dir=OUTPUT_DIR / SAMPLE / "functional" / "interpro",
        bakta_dir=OUTPUT_DIR / SAMPLE / "annotation" / "bakta"
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "aggregate_functional.log"
    script:
        "../scripts/aggregate_functional.py"


# ============ Prepare MAG Upload ============
rule prepare_mag_upload:
    """Prepare lightweight data for server upload (no full sequences)"""
    input:
        functional_summary=OUTPUT_DIR / SAMPLE / "functional" / "summary" / "functional_summary.json",
        checkm_quality=OUTPUT_DIR / SAMPLE / "binning" / "checkm2" / "quality_report.tsv",
        gtdbtk_summary=OUTPUT_DIR / SAMPLE / "annotation" / "gtdbtk" / "gtdbtk.bac120.summary.tsv",
        hq_bins_dir=OUTPUT_DIR / SAMPLE / "binning" / "hq_bins"
    output:
        upload_json=OUTPUT_DIR / SAMPLE / "upload" / "mags_upload.json"
    params:
        sample=SAMPLE,
        project=config.get("project_name", ""),
        subproject=config.get("subproject", ""),
        date=config.get("date", ""),
        bakta_dir=OUTPUT_DIR / SAMPLE / "annotation" / "bakta",
        functional_dir=OUTPUT_DIR / SAMPLE / "functional" / "summary"
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "prepare_mag_upload.log"
    script:
        "../scripts/prepare_mag_upload.py"


# ============ Upload MAGs to Server ============
rule upload_mags:
    """Upload MAG data to MMonitor server"""
    input:
        upload_json=OUTPUT_DIR / SAMPLE / "upload" / "mags_upload.json"
    output:
        status=OUTPUT_DIR / SAMPLE / "upload" / "upload_complete.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "upload_mags.log"
    script:
        "../scripts/upload_mags.py"


# ============ Full Functional Pipeline ============
rule functional_analysis:
    """Run complete functional annotation pipeline with upload"""
    input:
        upload_status=OUTPUT_DIR / SAMPLE / "upload" / "upload_complete.json"
    output:
        complete=OUTPUT_DIR / SAMPLE / "functional" / "functional_analysis_complete.json"
    shell:
        """
        echo '{{"status": "complete", "pipeline": "functional_analysis"}}' > {output.complete}
        """
