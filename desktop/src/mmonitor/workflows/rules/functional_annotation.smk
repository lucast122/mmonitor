# Functional Annotation Pipeline
# Uses eggNOG-mapper for COG/KEGG/GO annotations
# Uses InterProScan for PFAM/TIGRFAM domain annotations
# All rules use {sample} wildcard for multi-sample support

# ============ eggNOG-mapper ============
rule eggnog_mapper:
    """Run eggNOG-mapper for COG, KEGG, GO annotations"""
    input:
        proteins=OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.faa",
        db_ready=DB_BASE / "eggnog" / ".db_ready"
    output:
        annotations=OUTDIR + "/{sample}/functional/eggnog/{bin_id}.emapper.annotations",
        orthologs=OUTDIR + "/{sample}/functional/eggnog/{bin_id}.emapper.seed_orthologs"
    params:
        db=get_db_path("eggnog"),
        output_prefix=lambda wildcards: OUTDIR + f"/{wildcards.sample}/functional/eggnog/{wildcards.bin_id}",
        sensmode=config.get("eggnog", {}).get("sensmode", "default")
    threads: get_threads("eggnog")
    log:
        OUTDIR + "/{sample}/logs/eggnog_{bin_id}.log"
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
        proteins=OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.faa"
    output:
        tsv=OUTDIR + "/{sample}/functional/interpro/{bin_id}.tsv",
        gff=OUTDIR + "/{sample}/functional/interpro/{bin_id}.gff3"
    params:
        output_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/functional/interpro",
        applications=config.get("interproscan", {}).get("applications", "Pfam,TIGRFAM,CDD"),
        interproscan_path=config.get("interproscan", {}).get("path", "interproscan.sh")
    threads: get_threads("interproscan")
    log:
        OUTDIR + "/{sample}/logs/interproscan_{bin_id}.log"
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
    checkpoint_output = checkpoints.filter_hq_bins.get(sample=wildcards.sample).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTDIR + "/{sample}/functional/eggnog/{bin_id}.emapper.annotations",
        sample=wildcards.sample, bin_id=bin_ids
    )

def get_all_interpro_outputs(wildcards):
    """Get all InterProScan output files for HQ bins"""
    checkpoint_output = checkpoints.filter_hq_bins.get(sample=wildcards.sample).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTDIR + "/{sample}/functional/interpro/{bin_id}.tsv",
        sample=wildcards.sample, bin_id=bin_ids
    )

def get_all_bakta_gff(wildcards):
    """Get all Bakta GFF files for HQ bins"""
    checkpoint_output = checkpoints.filter_hq_bins.get(sample=wildcards.sample).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.gff3",
        sample=wildcards.sample, bin_id=bin_ids
    )


rule aggregate_functional:
    """Aggregate eggNOG + InterProScan results into summaries"""
    input:
        eggnog_files=get_all_eggnog_outputs,
        interpro_files=get_all_interpro_outputs,
        bakta_gff_files=get_all_bakta_gff,
        gtdbtk_summary=OUTDIR + "/{sample}/annotation/gtdbtk/gtdbtk.bac120.summary.tsv",
        checkm_quality=OUTDIR + "/{sample}/binning/checkm2/quality_report.tsv"
    output:
        summary=OUTDIR + "/{sample}/functional/summary/functional_summary.json",
        annotations=OUTDIR + "/{sample}/functional/summary/all_annotations.tsv"
    params:
        sample=lambda wildcards: wildcards.sample,
        project=config.get("project_name", ""),
        eggnog_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/functional/eggnog",
        interpro_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/functional/interpro",
        bakta_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/annotation/bakta"
    log:
        OUTDIR + "/{sample}/logs/aggregate_functional.log"
    script:
        "../scripts/aggregate_functional.py"


# ============ Prepare MAG Upload ============
rule prepare_mag_upload:
    """Prepare lightweight data for server upload (no full sequences)"""
    input:
        functional_summary=OUTDIR + "/{sample}/functional/summary/functional_summary.json",
        checkm_quality=OUTDIR + "/{sample}/binning/checkm2/quality_report.tsv",
        gtdbtk_summary=OUTDIR + "/{sample}/annotation/gtdbtk/gtdbtk.bac120.summary.tsv",
        hq_bins_dir=OUTDIR + "/{sample}/binning/hq_bins"
    output:
        upload_json=OUTDIR + "/{sample}/upload/mags_upload.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        project=config.get("project_name", ""),
        subproject=config.get("subproject", ""),
        date=config.get("date", ""),
        bakta_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/annotation/bakta",
        functional_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/functional/summary"
    log:
        OUTDIR + "/{sample}/logs/prepare_mag_upload.log"
    script:
        "../scripts/prepare_mag_upload.py"


# ============ Upload MAGs to Server ============
rule upload_mags:
    """Upload MAG data to MMonitor server"""
    input:
        upload_json=OUTDIR + "/{sample}/upload/mags_upload.json"
    output:
        status=OUTDIR + "/{sample}/upload/upload_complete.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTDIR + "/{sample}/logs/upload_mags.log"
    script:
        "../scripts/upload_mags.py"


# ============ Full Functional Pipeline ============
rule functional_analysis:
    """Run complete functional annotation pipeline with upload"""
    input:
        upload_status=OUTDIR + "/{sample}/upload/upload_complete.json"
    output:
        complete=OUTDIR + "/{sample}/functional/functional_analysis_complete.json"
    shell:
        """
        echo '{{"status": "complete", "pipeline": "functional_analysis"}}' > {output.complete}
        """
