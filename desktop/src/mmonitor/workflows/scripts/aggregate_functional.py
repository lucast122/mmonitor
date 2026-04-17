#!/usr/bin/env python3
"""
Aggregate functional annotations from eggNOG-mapper and InterProScan.

Creates:
1. A summary JSON with COG/PFAM/KEGG category distributions
2. A combined TSV with all gene annotations

Input files:
- eggNOG-mapper .emapper.annotations files
- InterProScan .tsv files
- Bakta .gff3 files (for gene locations)

Output:
- functional_summary.json: Aggregated statistics per MAG
- all_annotations.tsv: Combined annotations for database upload
"""

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# COG functional categories
COG_CATEGORIES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'X': 'Mobilome: prophages, transposons',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
}


def parse_eggnog_annotations(filepath):
    """Parse eggNOG-mapper annotation file."""
    annotations = {}

    if not os.path.exists(filepath):
        return annotations

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 21:
                    continue

                query = parts[0]
                annotations[query] = {
                    'seed_ortholog': parts[1] if len(parts) > 1 else '',
                    'evalue': float(parts[2]) if len(parts) > 2 and parts[2] else None,
                    'score': float(parts[3]) if len(parts) > 3 and parts[3] else None,
                    'cog': parts[4] if len(parts) > 4 else '',
                    'kegg_ko': parts[11] if len(parts) > 11 else '',
                    'kegg_pathway': parts[12] if len(parts) > 12 else '',
                    'kegg_module': parts[13] if len(parts) > 13 else '',
                    'kegg_reaction': parts[14] if len(parts) > 14 else '',
                    'kegg_rclass': parts[15] if len(parts) > 15 else '',
                    'kegg_brite': parts[16] if len(parts) > 16 else '',
                    'kegg_tc': parts[17] if len(parts) > 17 else '',
                    'cazy': parts[18] if len(parts) > 18 else '',
                    'pfam': parts[20] if len(parts) > 20 else '',
                    'go_terms': parts[9].split(',') if len(parts) > 9 and parts[9] else [],
                    'ec_number': parts[10] if len(parts) > 10 else '',
                    'description': parts[7] if len(parts) > 7 else '',
                }
    except Exception as e:
        print(f"Warning: Error parsing {filepath}: {e}", file=sys.stderr)

    return annotations


def parse_interpro_annotations(filepath):
    """Parse InterProScan TSV output file."""
    annotations = {}

    if not os.path.exists(filepath):
        return annotations

    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 11:
                    continue

                query = parts[0]
                if query not in annotations:
                    annotations[query] = {
                        'pfam': [],
                        'tigrfam': [],
                        'cdd': [],
                        'go_terms': [],
                        'pathways': [],
                    }

                analysis = parts[3]
                accession = parts[4]
                description = parts[5]

                if analysis == 'Pfam':
                    annotations[query]['pfam'].append({
                        'accession': accession,
                        'description': description,
                        'start': int(parts[6]) if parts[6] else None,
                        'end': int(parts[7]) if parts[7] else None,
                    })
                elif analysis == 'TIGRFAM':
                    annotations[query]['tigrfam'].append({
                        'accession': accession,
                        'description': description,
                    })
                elif analysis == 'CDD':
                    annotations[query]['cdd'].append({
                        'accession': accession,
                        'description': description,
                    })

                # GO terms are in column 14
                if len(parts) > 13 and parts[13]:
                    go_terms = parts[13].split('|')
                    annotations[query]['go_terms'].extend(go_terms)

                # Pathways in column 15
                if len(parts) > 14 and parts[14]:
                    pathways = parts[14].split('|')
                    annotations[query]['pathways'].extend(pathways)

    except Exception as e:
        print(f"Warning: Error parsing {filepath}: {e}", file=sys.stderr)

    return annotations


def parse_bakta_gff(filepath):
    """Parse Bakta GFF3 to get gene locations."""
    genes = {}

    if not os.path.exists(filepath):
        return genes

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('>'):
                    continue
                if '##FASTA' in line:
                    break

                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                if feature_type not in ['CDS', 'gene', 'tRNA', 'rRNA']:
                    continue

                attributes = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value

                gene_id = attributes.get('ID', '')
                if not gene_id:
                    continue

                genes[gene_id] = {
                    'contig': parts[0],
                    'type': feature_type,
                    'start': int(parts[3]),
                    'stop': int(parts[4]),
                    'strand': parts[6],
                    'product': attributes.get('product', attributes.get('Name', '')),
                }
    except Exception as e:
        print(f"Warning: Error parsing {filepath}: {e}", file=sys.stderr)

    return genes


def aggregate_cog_categories(annotations):
    """Count COG category occurrences."""
    cog_counts = defaultdict(int)

    for gene_id, annot in annotations.items():
        cogs = annot.get('cog', '')
        if cogs:
            for cog in cogs:
                if cog in COG_CATEGORIES:
                    cog_counts[cog] += 1

    return dict(cog_counts)


def aggregate_kegg_pathways(annotations):
    """Count KEGG pathway occurrences."""
    pathway_counts = defaultdict(int)

    for gene_id, annot in annotations.items():
        kegg_ko = annot.get('kegg_ko', '')
        if kegg_ko:
            for ko in kegg_ko.split(','):
                ko = ko.strip()
                if ko:
                    pathway_counts[ko] += 1

    return dict(pathway_counts)


def aggregate_pfam_domains(annotations):
    """Count PFAM domain occurrences."""
    pfam_counts = defaultdict(int)

    for gene_id, annot in annotations.items():
        # From eggNOG
        pfam_str = annot.get('pfam', '')
        if isinstance(pfam_str, str) and pfam_str:
            for pfam in pfam_str.split(','):
                pfam = pfam.strip()
                if pfam:
                    pfam_counts[pfam] += 1
        # From InterProScan
        elif isinstance(pfam_str, list):
            for pfam_entry in pfam_str:
                if isinstance(pfam_entry, dict):
                    pfam_counts[pfam_entry.get('accession', '')] += 1

    return dict(pfam_counts)


def main():
    """Main function to aggregate functional annotations."""
    # Get inputs from snakemake
    eggnog_files = snakemake.input.eggnog_files
    interpro_files = snakemake.input.interpro_files
    bakta_gff_files = snakemake.input.bakta_gff_files

    output_summary = snakemake.output.summary
    output_annotations = snakemake.output.annotations

    sample = snakemake.params.sample
    project = snakemake.params.project

    # Process each MAG
    all_mags = {}
    all_annotations_rows = []

    # Get bin IDs from file paths
    bin_ids = set()
    for f in eggnog_files:
        bin_id = Path(f).stem.replace('.emapper', '')
        bin_ids.add(bin_id)

    for bin_id in bin_ids:
        # Find files for this bin
        eggnog_file = None
        interpro_file = None
        bakta_gff = None

        for f in eggnog_files:
            if bin_id in str(f):
                eggnog_file = f
                break

        for f in interpro_files:
            if bin_id in str(f):
                interpro_file = f
                break

        for f in bakta_gff_files:
            if bin_id in str(f):
                bakta_gff = f
                break

        # Parse annotations
        eggnog_annots = parse_eggnog_annotations(eggnog_file) if eggnog_file else {}
        interpro_annots = parse_interpro_annotations(interpro_file) if interpro_file else {}
        bakta_genes = parse_bakta_gff(bakta_gff) if bakta_gff else {}

        # Merge annotations
        merged = {}
        for gene_id in set(list(eggnog_annots.keys()) + list(interpro_annots.keys())):
            merged[gene_id] = {
                **eggnog_annots.get(gene_id, {}),
                **interpro_annots.get(gene_id, {}),
            }

        # Calculate summaries
        cog_distribution = aggregate_cog_categories(merged)
        kegg_distribution = aggregate_kegg_pathways(merged)
        pfam_distribution = aggregate_pfam_domains(merged)

        all_mags[bin_id] = {
            'cog_distribution': cog_distribution,
            'cog_categories': {cat: COG_CATEGORIES.get(cat, 'Unknown') for cat in cog_distribution.keys()},
            'kegg_kos': kegg_distribution,
            'pfam_domains': pfam_distribution,
            'total_genes_annotated': len(merged),
            'genes_with_cog': sum(1 for g in merged.values() if g.get('cog')),
            'genes_with_kegg': sum(1 for g in merged.values() if g.get('kegg_ko')),
            'genes_with_pfam': sum(1 for g in merged.values() if g.get('pfam')),
            'genes_with_go': sum(1 for g in merged.values() if g.get('go_terms')),
        }

        # Create annotation rows for TSV
        for gene_id, annot in merged.items():
            gene_info = bakta_genes.get(gene_id, {})

            # Get primary PFAM
            pfam = annot.get('pfam', '')
            pfam_desc = ''
            if isinstance(pfam, list) and pfam:
                pfam = pfam[0].get('accession', '') if isinstance(pfam[0], dict) else ''
                pfam_desc = pfam[0].get('description', '') if isinstance(pfam[0], dict) else ''
            elif isinstance(pfam, str):
                pfam = pfam.split(',')[0].strip() if pfam else ''

            row = {
                'bin_id': bin_id,
                'gene_id': gene_id,
                'contig': gene_info.get('contig', ''),
                'start': gene_info.get('start', ''),
                'stop': gene_info.get('stop', ''),
                'strand': gene_info.get('strand', ''),
                'product': gene_info.get('product', annot.get('description', '')),
                'cog': annot.get('cog', ''),
                'cog_category': ','.join([COG_CATEGORIES.get(c, '') for c in annot.get('cog', '') if c in COG_CATEGORIES]),
                'kegg': annot.get('kegg_ko', ''),
                'ec_number': annot.get('ec_number', ''),
                'pfam': pfam,
                'pfam_description': pfam_desc,
                'go_terms': '|'.join(annot.get('go_terms', [])) if isinstance(annot.get('go_terms'), list) else annot.get('go_terms', ''),
            }
            all_annotations_rows.append(row)

    # Write summary JSON
    summary = {
        'sample': sample,
        'project': project,
        'mags': all_mags,
        'total_mags': len(all_mags),
    }

    os.makedirs(os.path.dirname(output_summary), exist_ok=True)
    with open(output_summary, 'w') as f:
        json.dump(summary, f, indent=2)

    # Write annotations TSV
    if all_annotations_rows:
        df = pd.DataFrame(all_annotations_rows)
        df.to_csv(output_annotations, sep='\t', index=False)
    else:
        # Create empty file with headers
        with open(output_annotations, 'w') as f:
            f.write('\t'.join(['bin_id', 'gene_id', 'contig', 'start', 'stop', 'strand',
                             'product', 'cog', 'cog_category', 'kegg', 'ec_number',
                             'pfam', 'pfam_description', 'go_terms']) + '\n')

    print(f"Aggregated functional annotations for {len(all_mags)} MAGs")


if __name__ == '__main__':
    main()
