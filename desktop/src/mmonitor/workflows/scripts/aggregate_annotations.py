#!/usr/bin/env python3
"""
Aggregate all annotation results from GTDB-TK and Bakta.
"""
import json
import pandas as pd
from pathlib import Path
from datetime import datetime


def parse_gtdbtk_results(gtdbtk_dir: Path) -> dict:
    """Parse GTDB-TK classification results."""
    results = {
        'bacteria': [],
        'archaea': []
    }

    # Parse bacterial classifications
    bac_summary = gtdbtk_dir / "gtdbtk.bac120.summary.tsv"
    if bac_summary.exists():
        df = pd.read_csv(bac_summary, sep='\t')
        for _, row in df.iterrows():
            results['bacteria'].append({
                'bin_id': row['user_genome'],
                'classification': row.get('classification', ''),
                'fastani_reference': row.get('fastani_reference', ''),
                'fastani_ani': row.get('fastani_ani', None),
                'fastani_af': row.get('fastani_af', None),
                'closest_placement_reference': row.get('closest_placement_reference', ''),
                'closest_placement_taxonomy': row.get('closest_placement_taxonomy', ''),
                'red_value': row.get('red_value', None),
                'note': row.get('note', '')
            })

    # Parse archaeal classifications
    arc_summary = gtdbtk_dir / "gtdbtk.ar53.summary.tsv"
    if arc_summary.exists():
        df = pd.read_csv(arc_summary, sep='\t')
        for _, row in df.iterrows():
            results['archaea'].append({
                'bin_id': row['user_genome'],
                'classification': row.get('classification', ''),
                'fastani_reference': row.get('fastani_reference', ''),
                'fastani_ani': row.get('fastani_ani', None),
                'fastani_af': row.get('fastani_af', None),
                'closest_placement_reference': row.get('closest_placement_reference', ''),
                'closest_placement_taxonomy': row.get('closest_placement_taxonomy', ''),
                'red_value': row.get('red_value', None),
                'note': row.get('note', '')
            })

    return results


def parse_bakta_results(bakta_dir: Path) -> dict:
    """Parse Bakta annotation results for all bins."""
    results = {}

    for bin_dir in bakta_dir.iterdir():
        if not bin_dir.is_dir():
            continue

        bin_id = bin_dir.name
        tsv_file = bin_dir / f"{bin_id}.tsv"

        if not tsv_file.exists():
            continue

        # Parse Bakta TSV output
        df = pd.read_csv(tsv_file, sep='\t', comment='#')

        # Count features by type
        feature_counts = df['Type'].value_counts().to_dict() if 'Type' in df.columns else {}

        # Get gene annotations
        genes = []
        for _, row in df.iterrows():
            gene_info = {
                'locus_tag': row.get('Locus Tag', ''),
                'type': row.get('Type', ''),
                'start': int(row.get('Start', 0)),
                'stop': int(row.get('Stop', 0)),
                'strand': row.get('Strand', ''),
                'gene': row.get('Gene', ''),
                'product': row.get('Product', ''),
                'db_xrefs': row.get('DbXrefs', '')
            }
            genes.append(gene_info)

        results[bin_id] = {
            'total_features': len(genes),
            'feature_counts': feature_counts,
            'genes': genes
        }

    return results


def aggregate_results(gtdbtk_results: dict, bakta_results: dict,
                     sample: str, project: str) -> dict:
    """Combine GTDB-TK and Bakta results into unified format."""

    # Build unified MAG list
    mags = []

    # Combine bacterial results
    for bac in gtdbtk_results.get('bacteria', []):
        bin_id = bac['bin_id']
        mag_info = {
            'bin_id': bin_id,
            'domain': 'Bacteria',
            'taxonomy': bac.get('classification', ''),
            'gtdbtk': {
                'classification': bac.get('classification', ''),
                'ani': bac.get('fastani_ani'),
                'af': bac.get('fastani_af'),
                'reference': bac.get('fastani_reference', '') or bac.get('closest_placement_reference', ''),
                'red_value': bac.get('red_value'),
                'note': bac.get('note', '')
            }
        }

        # Add Bakta annotation if available
        if bin_id in bakta_results:
            mag_info['annotation'] = {
                'total_features': bakta_results[bin_id]['total_features'],
                'feature_counts': bakta_results[bin_id]['feature_counts']
            }

        mags.append(mag_info)

    # Combine archaeal results
    for arc in gtdbtk_results.get('archaea', []):
        bin_id = arc['bin_id']
        mag_info = {
            'bin_id': bin_id,
            'domain': 'Archaea',
            'taxonomy': arc.get('classification', ''),
            'gtdbtk': {
                'classification': arc.get('classification', ''),
                'ani': arc.get('fastani_ani'),
                'af': arc.get('fastani_af'),
                'reference': arc.get('fastani_reference', '') or arc.get('closest_placement_reference', ''),
                'red_value': arc.get('red_value'),
                'note': arc.get('note', '')
            }
        }

        # Add Bakta annotation if available
        if bin_id in bakta_results:
            mag_info['annotation'] = {
                'total_features': bakta_results[bin_id]['total_features'],
                'feature_counts': bakta_results[bin_id]['feature_counts']
            }

        mags.append(mag_info)

    # Summary statistics
    total_cds = sum(
        bakta_results.get(m['bin_id'], {}).get('feature_counts', {}).get('CDS', 0)
        for m in mags
    )
    total_trna = sum(
        bakta_results.get(m['bin_id'], {}).get('feature_counts', {}).get('tRNA', 0)
        for m in mags
    )
    total_rrna = sum(
        bakta_results.get(m['bin_id'], {}).get('feature_counts', {}).get('rRNA', 0)
        for m in mags
    )

    return {
        'pipeline': 'annotation',
        'sample_name': sample,
        'project_name': project,
        'timestamp': datetime.now().isoformat(),
        'summary': {
            'total_mags': len(mags),
            'bacteria_count': len(gtdbtk_results.get('bacteria', [])),
            'archaea_count': len(gtdbtk_results.get('archaea', [])),
            'total_cds': total_cds,
            'total_trna': total_trna,
            'total_rrna': total_rrna
        },
        'mags': mags
    }


def main():
    # Snakemake provides these variables
    gtdbtk_dir = Path(snakemake.params.gtdbtk_dir)
    bakta_dir = Path(snakemake.params.bakta_dir)
    output_file = Path(snakemake.output.summary)
    sample = snakemake.params.sample
    project = snakemake.params.project

    print(f"Aggregating annotations for sample: {sample}")

    # Parse results
    gtdbtk_results = parse_gtdbtk_results(gtdbtk_dir)
    bakta_results = parse_bakta_results(bakta_dir)

    # Aggregate
    combined = aggregate_results(gtdbtk_results, bakta_results, sample, project)

    # Write output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(combined, f, indent=2)

    print(f"Aggregated {combined['summary']['total_mags']} MAGs")
    print(f"  - Bacteria: {combined['summary']['bacteria_count']}")
    print(f"  - Archaea: {combined['summary']['archaea_count']}")
    print(f"  - Total CDS: {combined['summary']['total_cds']}")
    print(f"Output written to: {output_file}")


if __name__ == '__main__':
    main()
