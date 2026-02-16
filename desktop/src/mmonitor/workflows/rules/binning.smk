# MAG Binning Pipeline
# Workflow: MetaBAT2 binning -> CheckM2 quality assessment

# ============ Map Reads to Assembly ============
rule map_reads_to_assembly:
    """Map reads back to polished assembly for coverage calculation"""
    input:
        reads=get_filtered_reads,
        assembly=OUTPUT_DIR / SAMPLE / "assembly" / "polished" / "consensus.fasta"
    output:
        bam=OUTPUT_DIR / SAMPLE / "binning" / "mapping" / "reads_to_assembly.sorted.bam",
        bai=OUTPUT_DIR / SAMPLE / "binning" / "mapping" / "reads_to_assembly.sorted.bam.bai"
    threads: get_threads("metabat2")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "map_reads_to_assembly.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.assembly} {input.reads} 2> {log} | \
            samtools sort -@ {threads} -o {output.bam} -

        samtools index {output.bam}
        """

# ============ Calculate Coverage Depth ============
rule jgi_summarize_bam:
    """Calculate contig depths for MetaBAT2"""
    input:
        bam=OUTPUT_DIR / SAMPLE / "binning" / "mapping" / "reads_to_assembly.sorted.bam"
    output:
        depth=OUTPUT_DIR / SAMPLE / "binning" / "coverage" / "depth.txt"
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "jgi_summarize_bam.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} 2>&1 | tee {log}
        """

# ============ MetaBAT2 Binning ============
rule metabat2_binning:
    """Bin contigs into MAGs using MetaBAT2"""
    input:
        assembly=OUTPUT_DIR / SAMPLE / "assembly" / "polished" / "consensus.fasta",
        depth=OUTPUT_DIR / SAMPLE / "binning" / "coverage" / "depth.txt"
    output:
        bins_dir=directory(OUTPUT_DIR / SAMPLE / "binning" / "bins"),
        done=OUTPUT_DIR / SAMPLE / "binning" / "metabat2.done"
    params:
        min_contig=config.get("metabat2", {}).get("min_contig", 2500),
        bin_prefix=OUTPUT_DIR / SAMPLE / "binning" / "bins" / "bin"
    threads: get_threads("metabat2")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "metabat2_binning.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        mkdir -p {output.bins_dir}

        metabat2 \
            -i {input.assembly} \
            -a {input.depth} \
            -o {params.bin_prefix} \
            --minContig {params.min_contig} \
            -t {threads} \
            2>&1 | tee {log}

        touch {output.done}
        """

# ============ CheckM2 Quality Assessment ============
rule checkm2_quality:
    """Assess MAG quality using CheckM2"""
    input:
        bins_dir=OUTPUT_DIR / SAMPLE / "binning" / "bins",
        done=OUTPUT_DIR / SAMPLE / "binning" / "metabat2.done",
        db_ready=DB_BASE / "checkm2" / ".db_ready"
    output:
        report=OUTPUT_DIR / SAMPLE / "binning" / "checkm2" / "quality_report.tsv"
    params:
        db=get_db_path("checkm2"),
        output_dir=OUTPUT_DIR / SAMPLE / "binning" / "checkm2"
    threads: get_threads("checkm2")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "checkm2_quality.log"
    conda:
        "../envs/checkm2.yaml"
    shell:
        """
        checkm2 predict \
            --input {input.bins_dir} \
            --output-directory {params.output_dir} \
            --threads {threads} \
            $([ -n "{params.db}" ] && echo "--database_path {params.db}") \
            -x fa \
            2>&1 | tee {log}

        # Rename output to expected name
        mv {params.output_dir}/quality_report.tsv {output.report} || true
        """

# ============ Filter High-Quality Bins ============
# NOTE: This MUST be a checkpoint (not a rule) because downstream rules in
# annotation.smk and functional_annotation.smk use
# checkpoints.filter_hq_bins.get() to dynamically resolve bin IDs.
checkpoint filter_hq_bins:
    """Filter bins by quality (completeness >= 50%, contamination <= 10%)"""
    input:
        bins_dir=OUTPUT_DIR / SAMPLE / "binning" / "bins",
        quality=OUTPUT_DIR / SAMPLE / "binning" / "checkm2" / "quality_report.tsv"
    output:
        hq_bins_dir=directory(OUTPUT_DIR / SAMPLE / "binning" / "hq_bins"),
        hq_list=OUTPUT_DIR / SAMPLE / "binning" / "hq_bins.txt"
    params:
        min_completeness=50.0,
        max_contamination=10.0
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "filter_hq_bins.log"
    script:
        "../scripts/filter_hq_bins.py"
