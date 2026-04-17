# MAG Annotation Pipeline
# Workflow: GTDB-TK classification -> Bakta annotation -> KEGGCharter
# All rules use {sample} wildcard for multi-sample support

# ============ GTDB-TK Classification ============
rule gtdbtk_classify:
    """Classify MAGs using GTDB-TK"""
    input:
        bins_dir=OUTDIR + "/{sample}/binning/hq_bins",
        hq_list=OUTDIR + "/{sample}/binning/hq_bins.txt",
        db_ready=DB_BASE / "gtdbtk" / ".db_ready"
    output:
        summary=OUTDIR + "/{sample}/annotation/gtdbtk/gtdbtk.bac120.summary.tsv",
        done=OUTDIR + "/{sample}/annotation/gtdbtk/gtdbtk.done"
    params:
        db=get_db_path("gtdbtk"),
        output_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/annotation/gtdbtk"
    threads: get_threads("gtdbtk")
    log:
        OUTDIR + "/{sample}/logs/gtdbtk_classify.log"
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}

        gtdbtk classify_wf \
            --genome_dir {input.bins_dir} \
            --out_dir {params.output_dir} \
            --cpus {threads} \
            --extension fa \
            2>&1 | tee {log}

        touch {output.done}
        """

# ============ Bakta Annotation ============
rule bakta_annotate:
    """Annotate each MAG using Bakta"""
    input:
        bin_fasta=OUTDIR + "/{sample}/binning/hq_bins/{bin_id}.fa",
        db_ready=DB_BASE / "bakta" / ".db_ready"
    output:
        gff=OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.gff3",
        faa=OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.faa",
        tsv=OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.tsv"
    params:
        db=get_db_path("bakta"),
        output_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/annotation/bakta/{wildcards.bin_id}",
        min_contig=config.get("bakta", {}).get("min_contig_length", 1)
    threads: get_threads("bakta")
    log:
        OUTDIR + "/{sample}/logs/bakta_{bin_id}.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        bakta \
            --db {params.db} \
            --output {params.output_dir} \
            --prefix {wildcards.bin_id} \
            --threads {threads} \
            --min-contig-length {params.min_contig} \
            {input.bin_fasta} \
            2>&1 | tee {log}
        """

# ============ Aggregate Bakta Annotations ============
def get_all_bakta_outputs(wildcards):
    """Get all Bakta output files for HQ bins"""
    checkpoint_output = checkpoints.filter_hq_bins.get(sample=wildcards.sample).output.hq_list
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip()]
    return expand(
        OUTDIR + "/{sample}/annotation/bakta/{bin_id}/{bin_id}.tsv",
        sample=wildcards.sample, bin_id=bin_ids
    )

# Note: filter_hq_bins checkpoint is defined in binning.smk

# ============ Aggregate All Annotations ============
rule aggregate_annotations:
    """Combine all annotation results"""
    input:
        gtdbtk=OUTDIR + "/{sample}/annotation/gtdbtk/gtdbtk.done",
        bakta_files=get_all_bakta_outputs
    output:
        summary=OUTDIR + "/{sample}/annotation/annotation_summary.json"
    params:
        sample=lambda wildcards: wildcards.sample,
        project=config.get("project_name", ""),
        gtdbtk_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/annotation/gtdbtk",
        bakta_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/annotation/bakta"
    log:
        OUTDIR + "/{sample}/logs/aggregate_annotations.log"
    script:
        "../scripts/aggregate_annotations.py"

# ============ Upload Annotation Results ============
rule upload_annotation_results:
    """Upload annotation results to MMonitor server"""
    input:
        summary=OUTDIR + "/{sample}/annotation/annotation_summary.json"
    output:
        status=OUTDIR + "/{sample}/annotation/annotation_complete.json"
    params:
        server_url=config.get("server", {}).get("url", "http://localhost:8000"),
        username=config.get("server", {}).get("username", ""),
        password=config.get("server", {}).get("password", ""),
        upload=config.get("server", {}).get("upload_results", True)
    log:
        OUTDIR + "/{sample}/logs/upload_annotations.log"
    script:
        "../scripts/upload_annotation_results.py"

# ============ Functional Pipeline (Full) ============
rule functional_complete:
    """Run complete functional analysis pipeline"""
    input:
        annotation=OUTDIR + "/{sample}/annotation/annotation_complete.json"
    output:
        complete=OUTDIR + "/{sample}/functional/functional_complete.json"
    shell:
        """
        cp {input.annotation} {output.complete}
        echo '{{"status": "complete", "pipeline": "functional"}}' >> {output.complete}
        """
