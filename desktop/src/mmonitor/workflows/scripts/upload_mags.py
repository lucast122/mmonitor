#!/usr/bin/env python3
"""
Upload MAG data to MMonitor server.

Uploads the lightweight MAG data prepared by prepare_mag_upload.py
to the MMonitor server API endpoints.
"""

import json
import os
import sys
from datetime import datetime

import requests


def authenticate(server_url, username, password):
    """Authenticate with the MMonitor server and get access token."""
    if not server_url.startswith('http'):
        server_url = f"https://{server_url}"
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


def upload_mag(server_url, token, mag_data):
    """Upload a single MAG to the server."""
    if not server_url.startswith('http'):
        server_url = f"https://{server_url}"
    mags_url = f"{server_url}/api/v1/mags/"
    headers = {'Authorization': f'Bearer {token}'}

    try:
        # First, create or update the MAG
        mag_payload = {
            'name': mag_data['name'],
            'sample_name': mag_data['sample_name'],
            'taxonomy': mag_data.get('taxonomy', {}),
            'completeness': mag_data.get('completeness'),
            'contamination': mag_data.get('contamination'),
            'strain_heterogeneity': mag_data.get('strain_heterogeneity'),
            'genome_size': mag_data.get('genome_size'),
            'n_contigs': mag_data.get('n_contigs'),
            'n50': mag_data.get('n50'),
            'gc_content': mag_data.get('gc_content'),
            'gene_count': mag_data.get('gene_count'),
            'cds_count': mag_data.get('cds_count'),
            'trna_count': mag_data.get('trna_count'),
            'rrna_count': mag_data.get('rrna_count'),
            'functional_summary': mag_data.get('functional_summary', {}),
        }

        # Check if MAG already exists
        check_url = f"{mags_url}?name={mag_data['name']}"
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
            else:
                # Create new MAG
                response = requests.post(mags_url, json=mag_payload, headers=headers)
        else:
            # Create new MAG
            response = requests.post(mags_url, json=mag_payload, headers=headers)

        response.raise_for_status()
        mag_response = response.json()
        mag_id = mag_response.get('id')

        print(f"  Uploaded MAG: {mag_data['name']} (ID: {mag_id})")

        # Upload contigs
        if mag_data.get('contigs'):
            contigs_url = f"{mags_url}{mag_id}/contigs/bulk_create/"
            contigs_payload = {
                'contigs': [
                    {
                        'contig_id': c['contig_id'],
                        'length': c['length'],
                        'gc_content': c.get('gc_content'),
                        'gene_count': c.get('gene_count', 0),
                    }
                    for c in mag_data['contigs']
                ]
            }

            try:
                contig_response = requests.post(
                    contigs_url,
                    json=contigs_payload,
                    headers=headers
                )
                if contig_response.status_code in [200, 201]:
                    print(f"    Uploaded {len(mag_data['contigs'])} contigs")
                else:
                    print(f"    Warning: Failed to upload contigs: {contig_response.text}", file=sys.stderr)
            except Exception as e:
                print(f"    Warning: Contig upload error: {e}", file=sys.stderr)

        # Upload annotations in batches
        if mag_data.get('annotations'):
            annotations = mag_data['annotations']
            batch_size = 500
            annotations_url = f"{mags_url}{mag_id}/annotations/bulk_create/"

            for i in range(0, len(annotations), batch_size):
                batch = annotations[i:i + batch_size]
                annotations_payload = {
                    'annotations': [
                        {
                            'gene_id': a['gene_id'],
                            'contig': a['contig'],
                            'start': a['start'],
                            'stop': a['stop'],
                            'strand': a['strand'],
                            'product': a['product'],
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
                        print(f"    Uploaded annotations batch {i // batch_size + 1}")
                    else:
                        print(f"    Warning: Failed to upload annotations batch: {annot_response.text}", file=sys.stderr)
                except Exception as e:
                    print(f"    Warning: Annotation upload error: {e}", file=sys.stderr)

            print(f"    Uploaded {len(annotations)} annotations total")

        return True

    except requests.exceptions.RequestException as e:
        print(f"Failed to upload MAG {mag_data['name']}: {e}", file=sys.stderr)
        return False


def main():
    """Main function to upload MAGs to server."""
    # Get inputs from snakemake
    upload_json_path = snakemake.input.upload_json
    output_path = snakemake.output.status

    server_url = snakemake.params.server_url
    username = snakemake.params.username
    password = snakemake.params.password
    do_upload = snakemake.params.upload

    # Load upload data
    with open(upload_json_path, 'r') as f:
        upload_data = json.load(f)

    result = {
        'status': 'skipped',
        'uploaded': 0,
        'failed': 0,
        'timestamp': datetime.now().isoformat(),
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

            # Upload each MAG
            mags = upload_data.get('mags', [])
            print(f"Uploading {len(mags)} MAGs...")

            for mag_data in mags:
                success = upload_mag(server_url, token, mag_data)
                if success:
                    result['uploaded'] += 1
                else:
                    result['failed'] += 1

            if result['failed'] == 0:
                result['status'] = 'success'
            elif result['uploaded'] > 0:
                result['status'] = 'partial'
            else:
                result['status'] = 'failed'

    print(f"\nUpload complete: {result['uploaded']} succeeded, {result['failed']} failed")

    # Write status
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)


if __name__ == '__main__':
    main()
