import os
import gzip
import shutil
import argparse
import concurrent.futures
import requests
from urllib.parse import urljoin
from Bio import SeqIO

def download_file(url, local_filename):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

def process_genome(taxid, ftp_path, output_dir):
    filename = os.path.basename(ftp_path)
    local_path = os.path.join(output_dir, filename)
    
    # Download the file
    download_file(ftp_path, local_path)
    
    # Process the file
    unzipped_path = local_path.replace('.gz', '')
    with gzip.open(local_path, 'rt') as f_in, open(unzipped_path, 'w') as f_out:
        sequences = SeqIO.parse(f_in, 'fasta')
        for record in sequences:
            record.id = f"kraken:taxid|{taxid}|{record.id}"
            SeqIO.write(record, f_out, 'fasta')
    
    # Clean up
    os.remove(local_path)
    
    return unzipped_path

def download_and_process_genomes(assembly_summary, output_dir, num_threads):
    os.makedirs(output_dir, exist_ok=True)
    
    with open(assembly_summary, 'r') as f:
        lines = f.readlines()[2:]  # Skip header lines
    
    genome_data = []
    for line in lines:
        fields = line.strip().split('\t')
        taxid = fields[5]
        ftp_path = fields[19]
        if ftp_path.startswith('ftp://'):
            ftp_path = 'https://' + ftp_path[6:]
        genome_file = os.path.basename(ftp_path) + '_genomic.fna.gz'
        full_url = urljoin(ftp_path + '/', genome_file)
        genome_data.append((taxid, full_url))
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        future_to_genome = {executor.submit(process_genome, taxid, url, output_dir): (taxid, url) for taxid, url in genome_data}
        for future in concurrent.futures.as_completed(future_to_genome):
            taxid, url = future_to_genome[future]
            try:
                processed_file = future.result()
                print(f"Processed: {processed_file}")
            except Exception as exc:
                print(f"Error processing genome {taxid}: {exc}")

def main():
    parser = argparse.ArgumentParser(description="Download and process genomes for Centrifuge database")
    parser.add_argument("assembly_summary", help="Path to the assembly summary file")
    parser.add_argument("output_dir", help="Directory to store processed genomes")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use for downloading and processing")
    args = parser.parse_args()

    download_and_process_genomes(args.assembly_summary, args.output_dir, args.threads)

if __name__ == "__main__":
    main()