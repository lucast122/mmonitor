#!/usr/bin/env python3
"""
Parse Centrifuger kreport output into standardized JSON format for MMonitor.
"""
import json
import pandas as pd
from pathlib import Path
from datetime import datetime


# Map centrifuger rank codes to taxonomy level names
RANK_MAP = {
    'D': 'superkingdom',
    'P': 'phylum',
    'C': 'class',
    'O': 'order',
    'F': 'family',
    'G': 'genus',
    'S': 'species',
}


def parse_kreport(kreport_file: Path, num_reads: int = 0) -> list:
    """Parse Kraken-style report into taxonomy records.

    The kreport has columns: percentage, clade_reads, taxon_reads, rank, tax_id, name.
    We extract species-level entries and build the taxonomy hierarchy by tracking
    the current lineage as we walk through the report.
    """
    columns = ['percentage', 'clade_reads', 'taxon_reads', 'rank', 'tax_id', 'name']
    df = pd.read_csv(kreport_file, sep='\t', names=columns, header=None)

    # Track current lineage as we walk through the hierarchical report
    current_lineage = {}
    rank_order = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    results = []
    for _, row in df.iterrows():
        rank_code = str(row['rank']).strip()
        name = str(row['name']).strip()
        percentage = float(row['percentage'])
        clade_reads = int(row['clade_reads'])
        tax_id = str(row['tax_id'])

        # Map rank code to taxonomy level
        rank_name = RANK_MAP.get(rank_code, '')

        if rank_name:
            # Update current lineage at this rank
            current_lineage[rank_name] = name
            # Clear lower ranks (they'll be set by subsequent rows)
            rank_idx = rank_order.index(rank_name) if rank_name in rank_order else -1
            for lower_rank in rank_order[rank_idx + 1:]:
                current_lineage.pop(lower_rank, None)

        # Only emit species-level entries with reads
        if rank_code == 'S' and percentage > 0:
            abundance = percentage / 100.0
            estimated_counts = round(abundance * num_reads) if num_reads > 0 else clade_reads

            result = {
                'tax_id': tax_id,
                'abundance': abundance,
                'estimated_counts': estimated_counts,
                'species': current_lineage.get('species', name),
                'genus': current_lineage.get('genus', ''),
                'family': current_lineage.get('family', ''),
                'order': current_lineage.get('order', ''),
                'class': current_lineage.get('class', ''),
                'phylum': current_lineage.get('phylum', ''),
                'superkingdom': current_lineage.get('superkingdom', ''),
            }
            results.append(result)

    return results


def parse_fastq_stats(stats_file: Path) -> dict:
    """Parse seqkit stats output."""
    df = pd.read_csv(stats_file, sep='\t')
    if len(df) == 0:
        return {}

    row = df.iloc[0]
    return {
        'num_reads': int(row.get('num_seqs', 0)),
        'total_bases': int(row.get('sum_len', 0)),
        'min_length': int(row.get('min_len', 0)),
        'avg_length': float(row.get('avg_len', 0)),
        'max_length': int(row.get('max_len', 0))
    }


def main():
    # Snakemake provides these variables
    kreport_file = Path(snakemake.input.kreport)
    stats_file = Path(snakemake.input.stats)
    output_file = Path(snakemake.output.json)

    # Parse stats first to get read count for estimated_counts
    stats = parse_fastq_stats(stats_file)
    num_reads = stats.get('num_reads', 0)

    # Parse kreport into taxonomy records
    taxonomy_results = parse_kreport(kreport_file, num_reads=num_reads)

    # Build output structure (same format as EMU parsed output)
    output = {
        'pipeline': 'taxonomy_wgs',
        'tool': 'centrifuger',
        'version': '1.0.12',
        'timestamp': datetime.now().isoformat(),
        'sample': {
            'name': snakemake.params.sample,
            'project': snakemake.params.project,
            'subproject': snakemake.params.subproject,
            'date': snakemake.params.date
        },
        'statistics': stats,
        'taxonomy': taxonomy_results
    }

    # Write output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)

if __name__ == '__main__':
    main()
