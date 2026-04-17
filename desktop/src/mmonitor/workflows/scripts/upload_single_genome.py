#!/usr/bin/env python3
"""
Upload single genome data to MMonitor server.

Uploads:
- Genome metadata (name, taxonomy, quality, etc.)
- FASTA file to media directory
- FASTA index (.fai) to media directory
- GFF file to media directory
- Contig summaries
- Gene annotations
"""

import json
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

import requests


def authenticate(server_url, username, password):
    """Authenticate with the MMonitor server and get access token."""
    auth_url = f"{server_url}/api/v1/auth/login/"

    try:
        response = requests.post(auth_url, json={
            'username': username,
            'password': password
        })
        response.raise_for_status()
        data = response.json()
        return data.get('access') or data.get('tokens', {}).get('access')
    except requests.exceptions.RequestException as e:
        print(f"Authentication failed: {e}", file=sys.stderr)
        return None


def copy_file_to_media(source_path, media_root, user_id, mag_id, filename):
    """Copy a file to the media directory and return the relative path."""
    # Create directory structure: media/mags/{user_id}/{mag_id}/
    dest_dir = Path(media_root) / 'mags' / str(user_id) / str(mag_id)
    dest_dir.mkdir(parents=True, exist_ok=True)

    dest_path = dest_dir / filename
    shutil.copy2(source_path, dest_path)

    # Return path relative to media root
    return f"mags/{user_id}/{mag_id}/{filename}"


def upload_genome(server_url, token, genome_data, fasta_path, fai_path, gff_path, media_root):
    """Upload a genome to the server with file paths."""
    mags_url = f"{server_url}/api/v1/mags/"
    headers = {'Authorization': f'Bearer {token}'}

    try:
        # First, get user ID
        user_url = f"{server_url}/api/v1/auth/user/"
        user_response = requests.get(user_url, headers=headers)
        user_response.raise_for_status()
        user_id = user_response.json().get('id')

        # Build MAG payload
        mag_payload = {
            'name': genome_data['name'],
            'sample_name': genome_data['sample_name'],
            'taxonomy': genome_data.get('taxonomy', {}),
            'completeness': genome_data.get('completeness'),
            'contamination': genome_data.get('contamination'),
            'strain_heterogeneity': genome_data.get('strain_heterogeneity'),
            'genome_size': genome_data.get('genome_size'),
            'n_contigs': genome_data.get('n_contigs'),
            'n50': genome_data.get('n50'),
            'gc_content': genome_data.get('gc_content'),
            'gene_count': genome_data.get('gene_count'),
            'cds_count': genome_data.get('cds_count'),
            'trna_count': genome_data.get('trna_count'),
            'rrna_count': genome_data.get('rrna_count'),
            'functional_summary': genome_data.get('functional_summary', {}),
        }

        # Check if MAG already exists
        check_url = f"{mags_url}?name={genome_data['name']}"
        check_response = requests.get(check_url, headers=headers)

        if check_response.status_code == 200:
            existing = check_response.json()
            if existing.get('results') and len(existing['results']) > 0:
                # Update existing MAG
                mag_id = existing['results'][0]['id']
                response = requests.patch(
                    f"{mags_url}{mag_id}/",
                    json=mag_payload,
                    headers=headers
                )
                print(f"  Updating existing MAG (ID: {mag_id})")
            else:
                # Create new MAG
                response = requests.post(mags_url, json=mag_payload, headers=headers)
                print(f"  Creating new MAG")
        else:
            # Create new MAG
            response = requests.post(mags_url, json=mag_payload, headers=headers)

        response.raise_for_status()
        mag_response = response.json()
        mag_id = mag_response.get('id')

        print(f"  MAG created/updated: {genome_data['name']} (ID: {mag_id})")

        # Copy files to media directory
        print(f"  Copying files to media directory...")

        genome_name = genome_data['name']

        # Copy FASTA
        if os.path.exists(fasta_path):
            fasta_rel_path = copy_file_to_media(
                fasta_path, media_root, user_id, mag_id, f"{genome_name}.fasta"
            )
            print(f"    FASTA: {fasta_rel_path}")
        else:
            fasta_rel_path = None
            print(f"    Warning: FASTA not found: {fasta_path}")

        # Copy FAI
        if os.path.exists(fai_path):
            fai_rel_path = copy_file_to_media(
                fai_path, media_root, user_id, mag_id, f"{genome_name}.fasta.fai"
            )
            print(f"    FAI: {fai_rel_path}")
        else:
            fai_rel_path = None
            print(f"    Warning: FAI not found: {fai_path}")

        # Copy GFF
        if os.path.exists(gff_path):
            gff_rel_path = copy_file_to_media(
                gff_path, media_root, user_id, mag_id, f"{genome_name}.gff3"
            )
            print(f"    GFF: {gff_rel_path}")
        else:
            gff_rel_path = None
            print(f"    Warning: GFF not found: {gff_path}")

        # Update MAG with file paths
        file_paths_payload = {
            'fasta_path': fasta_rel_path,
            'fai_path': fai_rel_path,
            'gff_path': gff_rel_path,
        }

        path_response = requests.patch(
            f"{mags_url}{mag_id}/",
            json=file_paths_payload,
            headers=headers
        )

        if path_response.status_code == 200:
            print(f"  File paths updated in database")
        else:
            print(f"  Warning: Failed to update file paths: {path_response.text}")

        # Upload contigs
        if genome_data.get('contigs'):
            contigs_url = f"{mags_url}{mag_id}/contigs/bulk_create/"
            contigs_payload = {
                'contigs': [
                    {
                        'contig_id': c['contig_id'],
                        'length': c['length'],
                        'gc_content': c.get('gc_content'),
                        'gene_count': c.get('gene_count', 0),
                    }
                    for c in genome_data['contigs']
                ]
            }

            try:
                contig_response = requests.post(
                    contigs_url,
                    json=contigs_payload,
                    headers=headers
                )
                if contig_response.status_code in [200, 201]:
                    print(f"  Uploaded {len(genome_data['contigs'])} contigs")
                else:
                    print(f"  Warning: Failed to upload contigs: {contig_response.text}", file=sys.stderr)
            except Exception as e:
                print(f"  Warning: Contig upload error: {e}", file=sys.stderr)

        # Upload annotations in batches
        if genome_data.get('annotations'):
            annotations = genome_data['annotations']
            batch_size = 500
            annotations_url = f"{mags_url}{mag_id}/annotations/bulk_create/"

            total_uploaded = 0
            for i in range(0, len(annotations), batch_size):
                batch = annotations[i:i + batch_size]
                annotations_payload = {
                    'annotations': [
                        {
                            'gene_id': a.get('gene_id', ''),
                            'contig': a.get('contig', ''),
                            'start': a.get('start', 0),
                            'stop': a.get('stop', 0),
                            'strand': a.get('strand', '+'),
                            'product': a.get('product', ''),
                            'cog': a.get('cog', ''),
                            'cog_category': a.get('cog_category', ''),
                            'kegg': a.get('kegg', ''),
                            'ec_number': a.get('ec_number', ''),
                            'pfam': a.get('pfam', ''),
                            'pfam_description': a.get('pfam_description', ''),
                            'go_terms': a.get('go_terms', []),
                        }
                        for a in batch
                    ]
                }

                try:
                    annot_response = requests.post(
                        annotations_url,
                        json=annotations_payload,
                        headers=headers
                    )
                    if annot_response.status_code in [200, 201]:
                        total_uploaded += len(batch)
                    else:
                        print(f"  Warning: Failed to upload annotation batch: {annot_response.text}", file=sys.stderr)
                except Exception as e:
                    print(f"  Warning: Annotation upload error: {e}", file=sys.stderr)

            print(f"  Uploaded {total_uploaded} annotations")

        return True, mag_id

    except requests.exceptions.RequestException as e:
        print(f"Failed to upload genome {genome_data['name']}: {e}", file=sys.stderr)
        return False, None


def main():
    """Main function to upload genome to server."""
    # Get inputs from snakemake
    upload_json_path = str(snakemake.input.upload_json)
    fasta_path = str(snakemake.input.fasta)
    fai_path = str(snakemake.input.fai)
    gff_path = str(snakemake.input.gff)

    output_path = str(snakemake.output.status)

    server_url = snakemake.params.server_url
    username = snakemake.params.username
    password = snakemake.params.password
    do_upload = snakemake.params.upload

    # Get media root from server URL (assume it's local)
    # This should be configured properly in production
    media_root = os.environ.get('MMONITOR_MEDIA_ROOT', '/home/minion-computer/mmonitor_dev/MMonitor/server/media')

    # Load upload data
    with open(upload_json_path, 'r') as f:
        genome_data = json.load(f)

    result = {
        'status': 'skipped',
        'mag_id': None,
        'timestamp': datetime.now().isoformat(),
        'genome_name': genome_data['name'],
    }

    if not do_upload:
        print("Upload disabled in config, skipping...")
        result['status'] = 'disabled'
    elif not username or not password:
        print("No credentials provided, skipping upload...")
        result['status'] = 'no_credentials'
    else:
        # Authenticate
        print(f"Authenticating with {server_url}...")
        token = authenticate(server_url, username, password)

        if not token:
            result['status'] = 'auth_failed'
        else:
            print(f"Authenticated successfully")

            # Upload genome
            print(f"Uploading genome: {genome_data['name']}")
            success, mag_id = upload_genome(
                server_url, token, genome_data,
                fasta_path, fai_path, gff_path,
                media_root
            )

            if success:
                result['status'] = 'success'
                result['mag_id'] = mag_id
            else:
                result['status'] = 'failed'

    print(f"\nUpload complete: {result['status']}")

    # Write status
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)


if __name__ == '__main__':
    main()
