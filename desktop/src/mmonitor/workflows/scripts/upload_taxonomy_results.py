#!/usr/bin/env python3
"""
Upload taxonomy results to MMonitor server.
"""
import json
import requests
from pathlib import Path
from datetime import datetime

def upload_to_server(results: dict, server_url: str, username: str, password: str) -> dict:
    """Upload results to MMonitor Django server."""
    # Ensure URL has scheme
    if not server_url.startswith('http'):
        server_url = f"https://{server_url}"

    # Login to get token
    login_url = f"{server_url}/api/v1/auth/login/"
    login_response = requests.post(login_url, json={
        'username': username,
        'password': password
    })

    if login_response.status_code != 200:
        raise Exception(f"Login failed: {login_response.text}")

    login_data = login_response.json()
    # Token may be at top level or nested under 'tokens'
    token = login_data.get('access') or login_data.get('tokens', {}).get('access')
    if not token:
        raise Exception(f"No access token in login response: {login_data}")
    headers = {'Authorization': f'Bearer {token}'}

    # Upload to the taxonomy upload endpoint
    endpoint = f"{server_url}/api/v1/taxonomy/upload/"
    response = requests.post(endpoint, json=results, headers=headers)

    if response.status_code not in [200, 201]:
        raise Exception(f"Upload failed: {response.text}")

    return response.json()

def main():
    # Snakemake provides these variables
    results_file = Path(snakemake.input.results)
    output_file = Path(snakemake.output.status)

    server_url = snakemake.params.server_url
    username = snakemake.params.username
    password = snakemake.params.password
    do_upload = snakemake.params.upload

    # Load results
    with open(results_file) as f:
        results = json.load(f)

    # Upload if enabled
    status = {
        'timestamp': datetime.now().isoformat(),
        'input_file': str(results_file),
        'uploaded': False,
        'server_url': server_url
    }

    if do_upload and username and password:
        try:
            response = upload_to_server(results, server_url, username, password)
            status['uploaded'] = True
            status['response'] = response
            print(f"Successfully uploaded results to {server_url}")
        except Exception as e:
            status['error'] = str(e)
            print(f"Warning: Upload failed - {e}")
    else:
        print("Skipping upload (disabled or no credentials)")

    # Write status file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(status, f, indent=2)

if __name__ == '__main__':
    main()
