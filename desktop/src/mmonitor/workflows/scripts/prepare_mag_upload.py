#!/usr/bin/env python3
"""
Prepare MAG data for upload to MMonitor server.

Creates a lightweight JSON that includes:
- Metadata: name, taxonomy, sample_name
- Quality: completeness, contamination, strain_heterogeneity
- Assembly stats: genome_size, n_contigs, n50, gc_content
- Gene counts: gene_count, cds_count, trna_count, rrna_count
- Functional summary: COG/PFAM/KEGG category aggregates
- Annotations: gene_id, contig, start, stop, strand, product, cog, pfam, kegg (NO sequences)
- Contig summaries: contig_id, length, gc_content

Does NOT include:
- Full FASTA sequences
- Full protein sequences
- Raw GFF content
"""

import json
import os
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd


def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    if not sequence:
        return None
    sequence = sequence.upper()
    gc = sum(1 for base in sequence if base in 'GC')
    total = sum(1 for base in sequence if base in 'ATGC')
    return (gc / total * 100) if total > 0 else None


def calculate_n50(lengths):
    """Calculate N50 from list of contig lengths."""
    if not lengths:
        return None
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total / 2:
            return length
    return sorted_lengths[0]


def parse_fasta_stats(fasta_path):
    """Parse FASTA file to get contig stats without storing sequences."""
    contigs = []
    current_contig = None
    current_seq = []
    genome_size = 0

    if not os.path.exists(fasta_path):
        return [], 0, None, None

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous contig
                if current_contig:
                    seq = ''.join(current_seq)
                    length = len(seq)
                    gc = calculate_gc_content(seq)
                    contigs.append({
                        'contig_id': current_contig,
                        'length': length,
                        'gc_content': gc,
                    })
                    genome_size += length

                # Start new contig
                current_contig = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last contig
        if current_contig:
            seq = ''.join(current_seq)
            length = len(seq)
            gc = calculate_gc_content(seq)
            contigs.append({
                'contig_id': current_contig,
                'length': length,
                'gc_content': gc,
            })
            genome_size += length

    # Calculate overall stats
    lengths = [c['length'] for c in contigs]
    n50 = calculate_n50(lengths)
    overall_gc = sum(c['gc_content'] * c['length'] for c in contigs if c['gc_content']) / genome_size if genome_size > 0 else None

    return contigs, genome_size, n50, overall_gc


def parse_checkm_quality(checkm_path, bin_id):
    """Parse CheckM2 quality report for a specific bin."""
    if not os.path.exists(checkm_path):
        return {}

    try:
        df = pd.read_csv(checkm_path, sep='\t')
        row = df[df['Name'] == bin_id]
        if row.empty:
            # Try without extension
            row = df[df['Name'] == bin_id.replace('.fa', '')]
        if row.empty:
            return {}

        return {
            'completeness': float(row['Completeness'].values[0]),
            'contamination': float(row['Contamination'].values[0]),
        }
    except Exception as e:
        print(f"Warning: Error parsing CheckM quality: {e}", file=sys.stderr)
        return {}


def parse_gtdbtk_taxonomy(gtdbtk_path, bin_id):
    """Parse GTDB-TK classification for a specific bin."""
    if not os.path.exists(gtdbtk_path):
        return {}

    try:
        df = pd.read_csv(gtdbtk_path, sep='\t')
        row = df[df['user_genome'] == bin_id]
        if row.empty:
            row = df[df['user_genome'] == bin_id.replace('.fa', '')]
        if row.empty:
            return {}

        classification = row['classification'].values[0]
        # Parse GTDB taxonomy string
        tax_dict = {}
        for level in classification.split(';'):
            if '__' in level:
                rank, name = level.split('__', 1)
                rank_map = {'d': 'domain', 'p': 'phylum', 'c': 'class',
                           'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'}
                if rank in rank_map and name:
                    tax_dict[rank_map[rank]] = name

        return tax_dict
    except Exception as e:
        print(f"Warning: Error parsing GTDB-TK taxonomy: {e}", file=sys.stderr)
        return {}


def parse_bakta_summary(bakta_tsv_path):
    """Parse Bakta TSV to get gene counts."""
    counts = {
        'gene_count': 0,
        'cds_count': 0,
        'trna_count': 0,
        'rrna_count': 0,
    }

    if not os.path.exists(bakta_tsv_path):
        return counts

    try:
        with open(bakta_tsv_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                feature_type = parts[1].lower() if len(parts) > 1 else ''
                counts['gene_count'] += 1

                if feature_type == 'cds':
                    counts['cds_count'] += 1
                elif feature_type == 'trna':
                    counts['trna_count'] += 1
                elif feature_type == 'rrna':
                    counts['rrna_count'] += 1
    except Exception as e:
        print(f"Warning: Error parsing Bakta TSV: {e}", file=sys.stderr)

    return counts


def parse_annotations_tsv(annotations_path, bin_id):
    """Parse combined annotations TSV for a specific bin."""
    annotations = []

    if not os.path.exists(annotations_path):
        return annotations

    try:
        df = pd.read_csv(annotations_path, sep='\t')
        bin_df = df[df['bin_id'] == bin_id]

        for _, row in bin_df.iterrows():
            annotations.append({
                'gene_id': row.get('gene_id', ''),
                'contig': row.get('contig', ''),
                'start': int(row['start']) if pd.notna(row.get('start')) else None,
                'stop': int(row['stop']) if pd.notna(row.get('stop')) else None,
                'strand': row.get('strand', ''),
                'product': row.get('product', ''),
                'cog': row.get('cog', ''),
                'cog_category': row.get('cog_category', ''),
                'kegg': row.get('kegg', ''),
                'ec_number': row.get('ec_number', ''),
                'pfam': row.get('pfam', ''),
                'pfam_description': row.get('pfam_description', ''),
                'go_terms': row.get('go_terms', '').split('|') if row.get('go_terms') else [],
            })
    except Exception as e:
        print(f"Warning: Error parsing annotations TSV: {e}", file=sys.stderr)

    return annotations


def main():
    """Main function to prepare MAG upload data."""
    # Get inputs from snakemake
    functional_summary_path = snakemake.input.functional_summary
    checkm_path = snakemake.input.checkm_quality
    gtdbtk_path = snakemake.input.gtdbtk_summary
    hq_bins_dir = snakemake.input.hq_bins_dir

    output_path = snakemake.output.upload_json

    sample = snakemake.params.sample
    project = snakemake.params.project
    subproject = snakemake.params.subproject
    date = snakemake.params.date
    bakta_dir = snakemake.params.bakta_dir
    functional_dir = snakemake.params.functional_dir

    # Load functional summary
    with open(functional_summary_path, 'r') as f:
        functional_summary = json.load(f)

    # Find all HQ bins
    bin_files = list(Path(hq_bins_dir).glob('*.fa'))

    mags_data = []

    for bin_file in bin_files:
        bin_id = bin_file.stem

        print(f"Processing MAG: {bin_id}")

        # Get FASTA stats
        contigs, genome_size, n50, gc_content = parse_fasta_stats(str(bin_file))

        # Get CheckM quality
        quality = parse_checkm_quality(checkm_path, bin_id)

        # Get GTDB-TK taxonomy
        taxonomy = parse_gtdbtk_taxonomy(gtdbtk_path, bin_id)

        # Get Bakta gene counts
        bakta_tsv = Path(bakta_dir) / bin_id / f"{bin_id}.tsv"
        gene_counts = parse_bakta_summary(str(bakta_tsv))

        # Get functional summary for this MAG
        mag_functional = functional_summary.get('mags', {}).get(bin_id, {})

        # Get annotations (without sequences)
        annotations_tsv = Path(functional_dir) / "all_annotations.tsv"
        annotations = parse_annotations_tsv(str(annotations_tsv), bin_id)

        # Count genes per contig
        contig_gene_counts = defaultdict(int)
        for annot in annotations:
            contig = annot.get('contig')
            if contig:
                contig_gene_counts[contig] += 1

        # Add gene counts to contigs
        for contig in contigs:
            contig['gene_count'] = contig_gene_counts.get(contig['contig_id'], 0)

        # Build MAG data structure
        mag_data = {
            'name': bin_id,
            'sample_name': sample,
            'project': project,
            'subproject': subproject,
            'date': date,

            # Taxonomy
            'taxonomy': taxonomy,

            # Quality metrics
            'completeness': quality.get('completeness'),
            'contamination': quality.get('contamination'),
            'strain_heterogeneity': None,  # Not available from CheckM2

            # Assembly stats
            'genome_size': genome_size,
            'n_contigs': len(contigs),
            'n50': n50,
            'gc_content': gc_content,

            # Gene counts
            'gene_count': gene_counts['gene_count'],
            'cds_count': gene_counts['cds_count'],
            'trna_count': gene_counts['trna_count'],
            'rrna_count': gene_counts['rrna_count'],

            # Functional summary (aggregated)
            'functional_summary': {
                'cog_distribution': mag_functional.get('cog_distribution', {}),
                'kegg_kos': mag_functional.get('kegg_kos', {}),
                'pfam_domains': mag_functional.get('pfam_domains', {}),
                'genes_with_cog': mag_functional.get('genes_with_cog', 0),
                'genes_with_kegg': mag_functional.get('genes_with_kegg', 0),
                'genes_with_pfam': mag_functional.get('genes_with_pfam', 0),
                'genes_with_go': mag_functional.get('genes_with_go', 0),
            },

            # Contig summaries (no sequences)
            'contigs': contigs,

            # Annotations (no sequences)
            'annotations': annotations,
        }

        mags_data.append(mag_data)

    # Write output
    output = {
        'sample': sample,
        'project': project,
        'subproject': subproject,
        'date': date,
        'mags': mags_data,
        'total_mags': len(mags_data),
    }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"Prepared upload data for {len(mags_data)} MAGs")


if __name__ == '__main__':
    main()
