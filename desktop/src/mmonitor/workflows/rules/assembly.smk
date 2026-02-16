# Metagenomic Assembly Pipeline
# Workflow: Flye assembly -> Medaka polishing

# ============ Flye Assembly ============
rule flye_assembly:
    """Assemble reads using Flye (metagenomic mode)"""
    input:
        fastq=get_filtered_reads
    output:
        assembly=OUTPUT_DIR / SAMPLE / "assembly" / "flye" / "assembly.fasta",
        info=OUTPUT_DIR / SAMPLE / "assembly" / "flye" / "assembly_info.txt",
        graph=OUTPUT_DIR / SAMPLE / "assembly" / "flye" / "assembly_graph.gfa"
    params:
        mode=config.get("flye", {}).get("mode", "nano-raw"),
        min_overlap=config.get("flye", {}).get("min_overlap", 1000),
        meta=config.get("flye", {}).get("meta", True),
        output_dir=OUTPUT_DIR / SAMPLE / "assembly" / "flye"
    threads: get_threads("flye")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "flye_assembly.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        flye \
            --{params.mode} {input.fastq} \
            --out-dir {params.output_dir} \
            --threads {threads} \
            --min-overlap {params.min_overlap} \
            $([ "{params.meta}" = "True" ] && echo "--meta") \
            2>&1 | tee {log}
        """

# ============ Medaka Polishing ============
rule medaka_consensus:
    """Polish assembly using Medaka with auto-detect model (fallback to config)"""
    input:
        reads=get_filtered_reads,
        assembly=OUTPUT_DIR / SAMPLE / "assembly" / "flye" / "assembly.fasta"
    output:
        consensus=OUTPUT_DIR / SAMPLE / "assembly" / "polished" / "consensus.fasta",
        bam=OUTPUT_DIR / SAMPLE / "assembly" / "polished" / "calls_to_draft.bam"
    params:
        wrapper=str(Path(workflow.basedir) / "scripts" / "medaka_wrapper.sh"),
        model=config.get("medaka", {}).get("model", "r1041_e82_400bps_hac_v5.0.0"),
        output_dir=OUTPUT_DIR / SAMPLE / "assembly" / "polished"
    threads: get_threads("medaka")
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "medaka_consensus.log"
    conda:
        "../envs/medaka.yaml"
    shell:
        """
        bash {params.wrapper} \
            {input.reads} \
            {input.assembly} \
            {params.output_dir} \
            {threads} \
            {params.model} \
            {log}
        """

# ============ Assembly Statistics ============
rule assembly_stats:
    """Calculate assembly statistics using seqkit"""
    input:
        raw=OUTPUT_DIR / SAMPLE / "assembly" / "flye" / "assembly.fasta",
        polished=OUTPUT_DIR / SAMPLE / "assembly" / "polished" / "consensus.fasta"
    output:
        stats=OUTPUT_DIR / SAMPLE / "assembly" / "assembly_stats.tsv"
    log:
        OUTPUT_DIR / SAMPLE / "logs" / "assembly_stats.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        echo -e "stage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > {output.stats}
        seqkit stats -T {input.raw} | tail -n 1 | sed 's/^/raw\t/' >> {output.stats}
        seqkit stats -T {input.polished} | tail -n 1 | sed 's/^/polished\t/' >> {output.stats}
        """
