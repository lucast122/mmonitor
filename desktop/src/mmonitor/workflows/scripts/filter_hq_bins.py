#!/usr/bin/env python3
"""
Filter MAG bins by quality metrics (completeness and contamination).
"""
import shutil
import pandas as pd
from pathlib import Path

def main():
    # Snakemake provides these variables
    bins_dir = Path(snakemake.input.bins_dir)
    quality_file = Path(snakemake.input.quality)
    output_dir = Path(snakemake.output.hq_bins_dir)
    output_list = Path(snakemake.output.hq_list)

    min_completeness = snakemake.params.min_completeness
    max_contamination = snakemake.params.max_contamination

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read CheckM2 quality report
    df = pd.read_csv(quality_file, sep='\t')

    # Filter by quality thresholds
    hq_bins = df[
        (df['Completeness'] >= min_completeness) &
        (df['Contamination'] <= max_contamination)
    ]

    # Copy high-quality bins to output directory
    hq_bin_ids = []
    for _, row in hq_bins.iterrows():
        bin_name = row['Name']
        hq_bin_ids.append(bin_name)

        # Find and copy the bin file
        for ext in ['.fa', '.fasta', '.fna']:
            src = bins_dir / f"{bin_name}{ext}"
            if src.exists():
                dst = output_dir / f"{bin_name}.fa"
                shutil.copy(src, dst)
                break

    # Write list of HQ bin IDs
    with open(output_list, 'w') as f:
        for bin_id in hq_bin_ids:
            f.write(f"{bin_id}\n")

    print(f"Filtered {len(hq_bin_ids)} high-quality bins from {len(df)} total bins")
    print(f"Criteria: completeness >= {min_completeness}%, contamination <= {max_contamination}%")

if __name__ == '__main__':
    main()
