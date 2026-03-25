#!/usr/bin/env python3
"""
Parse EMU abundance output into standardized JSON format for MMonitor.
"""
import json
import pandas as pd
from pathlib import Path
from datetime import datetime


def _safe_str(val):
    """Convert value to string, replacing NaN/None with empty string."""
    if pd.isna(val):
        return ''
    return str(val)


def load_emu_taxonomy(db_path: str) -> dict:
    """Load EMU taxonomy.tsv to map tax_id -> taxonomy info."""
    taxonomy_file = Path(db_path) / "taxonomy.tsv"
    if not taxonomy_file.exists():
        return {}

    df = pd.read_csv(taxonomy_file, sep='\t', dtype={'tax_id': str})
    tax_map = {}
    for _, row in df.iterrows():
        tid = _safe_str(row.get('tax_id', ''))
        if tid:
            tax_map[tid] = {
                'species': _safe_str(row.get('species', '')),
                'genus': _safe_str(row.get('genus', '')),
                'family': _safe_str(row.get('family', '')),
                'order': _safe_str(row.get('order', '')),
                'class': _safe_str(row.get('class', '')),
                'phylum': _safe_str(row.get('phylum', '')),
                'superkingdom': _safe_str(row.get('superkingdom', '')),
            }
    return tax_map


def parse_emu_abundance(abundance_file: Path, emu_db_path: str = '', num_reads: int = 0) -> list:
    """Parse EMU rel-abundance.tsv file."""
    df = pd.read_csv(abundance_file, sep='\t')

    # Load taxonomy lookup from EMU database
    tax_map = load_emu_taxonomy(emu_db_path) if emu_db_path else {}

    results = []
    for _, row in df.iterrows():
        abundance = row.get('abundance', 0)
        if pd.isna(abundance):
            abundance = 0.0

        tax_id = _safe_str(row.get('tax_id', ''))

        # Try to get taxonomy from the row itself first
        species = _safe_str(row.get('species', ''))
        genus = _safe_str(row.get('genus', ''))
        family = _safe_str(row.get('family', ''))
        order = _safe_str(row.get('order', ''))
        tax_class = _safe_str(row.get('class', ''))
        phylum = _safe_str(row.get('phylum', ''))
        superkingdom = _safe_str(row.get('superkingdom', ''))

        # If taxonomy fields are empty, look up from database
        if not species and tax_id in tax_map:
            info = tax_map[tax_id]
            species = info.get('species', '')
            genus = info.get('genus', '')
            family = info.get('family', '')
            order = info.get('order', '')
            tax_class = info.get('class', '')
            phylum = info.get('phylum', '')
            superkingdom = info.get('superkingdom', '')

        # Estimated counts = relative abundance * total reads
        estimated_counts = round(float(abundance) * num_reads) if num_reads > 0 else 0

        result = {
            'tax_id': tax_id,
            'abundance': float(abundance),
            'estimated_counts': estimated_counts,
            'species': species,
            'genus': genus,
            'family': family,
            'order': order,
            'class': tax_class,
            'phylum': phylum,
            'superkingdom': superkingdom,
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
    abundance_file = Path(snakemake.input.abundance)
    stats_file = Path(snakemake.input.stats)
    output_file = Path(snakemake.output.json)

    # Get EMU database path for taxonomy lookup
    emu_db_path = snakemake.params.get('emu_db', '')

    # Parse inputs — get stats first so we can compute estimated counts
    stats = parse_fastq_stats(stats_file)
    num_reads = stats.get('num_reads', 0)
    taxonomy_results = parse_emu_abundance(abundance_file, emu_db_path, num_reads=num_reads)

    # Build output structure
    output = {
        'pipeline': 'taxonomy_16s',
        'tool': 'emu',
        'version': '3.4.5',
        'timestamp': datetime.now().isoformat(),
        'sample': {
            'name': snakemake.params.sample,
            'project': snakemake.params.project,
            'subproject': snakemake.params.subproject,
            'date': snakemake.params.date
        },
        'parameters': {
            'min_abundance': snakemake.params.min_abundance
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
