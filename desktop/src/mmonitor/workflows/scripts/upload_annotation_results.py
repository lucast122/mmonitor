#!/usr/bin/env python3
"""
Upload annotation results to MMonitor server.
"""
import json
import requests
from pathlib import Path
from datetime import datetime


def upload_to_server(results: dict, server_url: str, username: str, password: str) -> dict:
    """Upload annotation results to MMonitor Django server."""
    # Login to get token
    login_url = f"{server_url}/api/auth/login/"
    login_response = requests.post(login_url, json={
        'username': username,
        'password': password
    })

    if login_response.status_code != 200:
        raise Exception(f"Login failed: {login_response.text}")

    login_data = login_response.json()
    token = login_data.get('access') or login_data.get('tokens', {}).get('access')
    if not token:
        raise Exception(f"No access token in login response: {login_data}")
    headers = {'Authorization': f'Bearer {token}'}

    # Upload annotation results
    endpoint = f"{server_url}/api/annotation/upload/"

    response = requests.post(endpoint, json=results, headers=headers)

    if response.status_code not in [200, 201]:
        raise Exception(f"Upload failed: {response.text}")

    return response.json()


def main():
    # Snakemake provides these variables
    summary_file = Path(snakemake.input.summary)
    output_file = Path(snakemake.output.status)

    server_url = snakemake.params.server_url
    username = snakemake.params.username
    password = snakemake.params.password
    do_upload = snakemake.params.upload

    # Load results
    with open(summary_file) as f:
        results = json.load(f)

    # Status tracking
    status = {
        'pipeline': 'annotation',
        'timestamp': datetime.now().isoformat(),
        'input_file': str(summary_file),
        'uploaded': False,
        'server_url': server_url,
        'summary': results.get('summary', {})
    }

    # Upload if enabled
    if do_upload and username and password:
        try:
            response = upload_to_server(results, server_url, username, password)
            status['uploaded'] = True
            status['response'] = response
            print(f"Successfully uploaded annotation results to {server_url}")
        except Exception as e:
            status['error'] = str(e)
            print(f"Warning: Upload failed - {e}")
    else:
        print("Skipping upload (disabled or no credentials)")

    # Write status file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(status, f, indent=2)

    print(f"Status written to: {output_file}")


if __name__ == '__main__':
    main()
