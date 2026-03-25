#!/usr/bin/env python3
"""
Prepare single genome data for upload to MMonitor server.

Parses Bakta annotation output and creates a JSON file with:
- Metadata: name, taxonomy, sample_name
- Assembly stats: genome_size, n_contigs, n50, gc_content
- Gene counts: gene_count, cds_count, trna_count, rrna_count
- Annotations: gene_id, contig, start, stop, strand, product
- Contig summaries: contig_id, length, gc_content
"""

import json
import os
import re
import sys
from pathlib import Path
from collections import defaultdict


def calculate_gc_content(sequence):
    """Calculate GC content of a sequence."""
    if not sequence:
        return None
    sequence = sequence.upper()
    gc = sum(1 for base in sequence if base in 'GC')
    total = sum(1 for base in sequence if base in 'ATGC')
    return round((gc / total * 100), 2) if total > 0 else None


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
        print(f"Warning: FASTA file not found: {fasta_path}", file=sys.stderr)
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
    overall_gc = round(sum(c['gc_content'] * c['length'] for c in contigs if c['gc_content']) / genome_size, 2) if genome_size > 0 else None

    return contigs, genome_size, n50, overall_gc


def parse_bakta_tsv(tsv_path):
    """Parse Bakta TSV to get gene counts and annotations."""
    counts = {
        'gene_count': 0,
        'cds_count': 0,
        'trna_count': 0,
        'rrna_count': 0,
    }
    annotations = []

    if not os.path.exists(tsv_path):
        print(f"Warning: Bakta TSV not found: {tsv_path}", file=sys.stderr)
        return counts, annotations

    try:
        with open(tsv_path, 'r') as f:
            header = None
            for line in f:
                if line.startswith('#'):
                    # Parse header
                    if 'Sequence Id' in line:
                        header = line[1:].strip().split('\t')
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue

                # Default column mapping for Bakta TSV
                # Columns: Sequence Id, Type, Start, Stop, Strand, Locus Tag, Gene, Product, DbXrefs
                seq_id = parts[0] if len(parts) > 0 else ''
                feature_type = parts[1].lower() if len(parts) > 1 else ''
                start = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 0
                stop = int(parts[3]) if len(parts) > 3 and parts[3].isdigit() else 0
                strand = parts[4] if len(parts) > 4 else '+'
                locus_tag = parts[5] if len(parts) > 5 else ''
                gene = parts[6] if len(parts) > 6 else ''
                product = parts[7] if len(parts) > 7 else ''
                dbxrefs = parts[8] if len(parts) > 8 else ''

                counts['gene_count'] += 1

                if feature_type == 'cds':
                    counts['cds_count'] += 1
                elif feature_type == 'trna':
                    counts['trna_count'] += 1
                elif feature_type in ('rrna', 'tmrna'):
                    counts['rrna_count'] += 1

                # Parse dbxrefs for COG, KEGG, etc.
                cog = ''
                kegg = ''
                pfam = ''
                ec_number = ''
                go_terms = []

                if dbxrefs:
                    for ref in dbxrefs.split(','):
                        ref = ref.strip()
                        if ref.startswith('COG:'):
                            cog = ref.replace('COG:', '')
                        elif ref.startswith('KEGG:'):
                            kegg = ref.replace('KEGG:', '')
                        elif ref.startswith('Pfam:'):
                            pfam = ref.replace('Pfam:', '')
                        elif ref.startswith('EC:'):
                            ec_number = ref.replace('EC:', '')
                        elif ref.startswith('GO:'):
                            go_terms.append(ref)

                annotations.append({
                    'gene_id': locus_tag or f"{seq_id}_{start}_{stop}",
                    'contig': seq_id,
                    'start': start,
                    'stop': stop,
                    'strand': strand,
                    'product': product or gene,
                    'cog': cog,
                    'kegg': kegg,
                    'pfam': pfam,
                    'ec_number': ec_number,
                    'go_terms': go_terms,
                })

    except Exception as e:
        print(f"Warning: Error parsing Bakta TSV: {e}", file=sys.stderr)

    return counts, annotations


def parse_gff_for_annotations(gff_path):
    """Parse GFF3 file for detailed annotations."""
    annotations = []

    if not os.path.exists(gff_path):
        return annotations

    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                seq_id = parts[0]
                feature_type = parts[2].lower()
                start = int(parts[3])
                stop = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # Only process CDS, tRNA, rRNA features
                if feature_type not in ('cds', 'trna', 'rrna', 'gene'):
                    continue

                # Parse attributes
                attrs = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attrs[key] = value

                gene_id = attrs.get('ID', attrs.get('locus_tag', f"{seq_id}_{start}_{stop}"))
                product = attrs.get('product', attrs.get('Name', ''))

                # Parse Dbxref for functional annotations
                cog = ''
                kegg = ''
                pfam = ''
                ec_number = ''
                go_terms = []

                if 'Dbxref' in attrs:
                    for ref in attrs['Dbxref'].split(','):
                        if ref.startswith('COG:'):
                            cog = ref.replace('COG:', '')
                        elif ref.startswith('KEGG:'):
                            kegg = ref.replace('KEGG:', '')
                        elif ref.startswith('Pfam:'):
                            pfam = ref.replace('Pfam:', '')

                if 'eC_number' in attrs:
                    ec_number = attrs['eC_number']

                annotations.append({
                    'gene_id': gene_id,
                    'contig': seq_id,
                    'start': start,
                    'stop': stop,
                    'strand': strand,
                    'product': product,
                    'cog': cog,
                    'kegg': kegg,
                    'pfam': pfam,
                    'ec_number': ec_number,
                    'go_terms': go_terms,
                })

    except Exception as e:
        print(f"Warning: Error parsing GFF: {e}", file=sys.stderr)

    return annotations


def build_functional_summary(annotations):
    """Build functional annotation summary from annotations."""
    cog_dist = defaultdict(int)
    kegg_kos = defaultdict(int)
    pfam_domains = defaultdict(int)

    genes_with_cog = 0
    genes_with_kegg = 0
    genes_with_pfam = 0
    genes_with_go = 0

    for ann in annotations:
        if ann.get('cog'):
            # COG categories are single letters
            for cat in ann['cog']:
                if cat.isalpha():
                    cog_dist[cat] += 1
            genes_with_cog += 1
        if ann.get('kegg'):
            kegg_kos[ann['kegg']] += 1
            genes_with_kegg += 1
        if ann.get('pfam'):
            pfam_domains[ann['pfam']] += 1
            genes_with_pfam += 1
        if ann.get('go_terms'):
            genes_with_go += 1

    return {
        'cog_distribution': dict(cog_dist),
        'kegg_kos': dict(sorted(kegg_kos.items(), key=lambda x: x[1], reverse=True)[:50]),
        'pfam_domains': dict(sorted(pfam_domains.items(), key=lambda x: x[1], reverse=True)[:50]),
        'genes_with_cog': genes_with_cog,
        'genes_with_kegg': genes_with_kegg,
        'genes_with_pfam': genes_with_pfam,
        'genes_with_go': genes_with_go,
    }


def main():
    """Main function to prepare genome upload data."""
    # Get inputs from snakemake
    fasta_path = str(snakemake.input.fasta)
    gff_path = str(snakemake.input.gff)
    tsv_path = str(snakemake.input.tsv)

    output_path = str(snakemake.output.json)

    sample = snakemake.params.sample
    genome_id = snakemake.params.genome_id
    taxonomy = snakemake.params.taxonomy
    project = snakemake.params.project

    print(f"Preparing upload data for genome: {genome_id}")

    # Parse FASTA stats
    contigs, genome_size, n50, gc_content = parse_fasta_stats(fasta_path)
    print(f"  FASTA: {len(contigs)} contigs, {genome_size} bp, N50={n50}")

    # Parse Bakta TSV for gene counts
    gene_counts, tsv_annotations = parse_bakta_tsv(tsv_path)
    print(f"  Bakta TSV: {gene_counts['gene_count']} genes")

    # Parse GFF for detailed annotations if TSV annotations are limited
    gff_annotations = parse_gff_for_annotations(gff_path)
    print(f"  GFF: {len(gff_annotations)} features")

    # Use whichever has more detailed annotations
    annotations = tsv_annotations if len(tsv_annotations) >= len(gff_annotations) else gff_annotations

    # Build functional summary
    functional_summary = build_functional_summary(annotations)

    # Count genes per contig
    contig_gene_counts = defaultdict(int)
    for annot in annotations:
        contig = annot.get('contig')
        if contig:
            contig_gene_counts[contig] += 1

    # Add gene counts to contigs
    for contig in contigs:
        contig['gene_count'] = contig_gene_counts.get(contig['contig_id'], 0)

    # Build genome data structure
    genome_data = {
        'name': genome_id,
        'sample_name': sample,
        'project': project,

        # Taxonomy (if provided)
        'taxonomy': taxonomy if taxonomy else {},

        # Quality metrics (reference genomes are assumed high quality)
        'completeness': 99.5,
        'contamination': 0.2,
        'strain_heterogeneity': None,

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

        # Functional summary
        'functional_summary': functional_summary,

        # Contig summaries (no sequences)
        'contigs': contigs,

        # Annotations (no sequences)
        'annotations': annotations,

        # File paths (will be set during upload)
        'fasta_path': fasta_path,
        'gff_path': gff_path,
    }

    # Write output
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(genome_data, f, indent=2)

    print(f"Prepared upload data: {output_path}")
    print(f"  Genome size: {genome_size:,} bp")
    print(f"  Contigs: {len(contigs)}")
    print(f"  Annotations: {len(annotations)}")


if __name__ == '__main__':
    main()
