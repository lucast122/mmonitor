# Metagenomic Assembly Pipeline
# Workflow: Flye assembly -> Medaka polishing
# All rules use {sample} wildcard for multi-sample support

# ============ Flye Assembly ============
rule flye_assembly:
    """Assemble reads using Flye (metagenomic mode)"""
    input:
        fastq=get_filtered_reads
    output:
        assembly=OUTDIR + "/{sample}/assembly/flye/assembly.fasta",
        info=OUTDIR + "/{sample}/assembly/flye/assembly_info.txt",
        graph=OUTDIR + "/{sample}/assembly/flye/assembly_graph.gfa"
    params:
        mode=config.get("flye", {}).get("mode", "nano-raw"),
        min_overlap=config.get("flye", {}).get("min_overlap", 1000),
        meta=config.get("flye", {}).get("meta", True),
        output_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/assembly/flye"
    threads: get_threads("flye")
    log:
        OUTDIR + "/{sample}/logs/flye_assembly.log"
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
        assembly=OUTDIR + "/{sample}/assembly/flye/assembly.fasta"
    output:
        consensus=OUTDIR + "/{sample}/assembly/polished/consensus.fasta",
        bam=OUTDIR + "/{sample}/assembly/polished/calls_to_draft.bam"
    params:
        wrapper=str(Path(workflow.basedir) / "scripts" / "medaka_wrapper.sh"),
        model=config.get("medaka", {}).get("model", "r1041_e82_400bps_hac_v5.0.0"),
        output_dir=lambda wildcards: OUTDIR + f"/{wildcards.sample}/assembly/polished"
    threads: get_threads("medaka")
    log:
        OUTDIR + "/{sample}/logs/medaka_consensus.log"
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
        raw=OUTDIR + "/{sample}/assembly/flye/assembly.fasta",
        polished=OUTDIR + "/{sample}/assembly/polished/consensus.fasta"
    output:
        stats=OUTDIR + "/{sample}/assembly/assembly_stats.tsv"
    log:
        OUTDIR + "/{sample}/logs/assembly_stats.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        echo -e "stage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len" > {output.stats}
        seqkit stats -T {input.raw} | tail -n 1 | sed 's/^/raw\t/' >> {output.stats}
        seqkit stats -T {input.polished} | tail -n 1 | sed 's/^/polished\t/' >> {output.stats}
        """
