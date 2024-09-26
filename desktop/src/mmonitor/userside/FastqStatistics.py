import gzip
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from Bio import SeqIO
import time
import json

def process_fastq_file(file_path):
    sequences_temp = []
    qualities_temp = []
    local_gc_counts = []
    local_lengths = []
    local_q20_counts = []
    local_q30_counts = []

    is_gzipped = file_path.endswith(".gz")
    with (gzip.open(file_path, 'rt') if is_gzipped else open(file_path, 'r')) as f:
        read_count = 0  # Initialize read count
        for record in SeqIO.parse(f, "fastq"):
            if read_count >= 20:  # Adjust the threshold as needed
                break  # Stop reading if the threshold is exceeded
            seq = str(record.seq)
            sequences_temp.append(seq)  # Optional
            quality_scores = np.array(record.letter_annotations["phred_quality"], dtype=int)
            qualities_temp.append(quality_scores)
            local_gc_counts.append(seq.count('G') + seq.count('C'))
            local_lengths.append(len(seq))
            read_count += 1  # Increment read count

    return {
        'sequences': sequences_temp,
        'qualities': qualities_temp,
        'gc_counts': local_gc_counts,
        'lengths': local_lengths,
        'q20_counts': local_q20_counts,
        'q30_counts': local_q30_counts
    }


class FastqStatistics:

    def __init__(self, file_paths, multi=True, num_threads=64):
        self.file_paths = file_paths if isinstance(file_paths, list) else [file_paths]
        # Initialize lists and numpy arrays for aggregated data
        self.sequences = []
        self.qualities = []
        self.lengths = np.array([], dtype=int)
        self.gc_counts = np.array([], dtype=int)
        self.q20_counts = np.array([], dtype=int)
        self.q30_counts = np.array([], dtype=int)
        self.total_bases = 0
        self.load_files(multi, num_threads)

    def load_files(self, multi, num_threads):
        if multi:
            self.load_files_parallel(num_threads)
        else:
            for file_path in self.file_paths:
                result = process_fastq_file(file_path)
                self.aggregate_results(result)

    def load_files_parallel(self, num_threads):
        start_time = time.time()
        with ProcessPoolExecutor(max_workers=num_threads) as executor:
            results = list(executor.map(process_fastq_file, self.file_paths))
        for result in results:
            self.aggregate_results(result)
        end_time = time.time()
        print(f"Loaded files in {end_time - start_time} seconds")

    def aggregate_results(self, result):
        # Aggregate results here, similar to the load_file method adjustments
        self.sequences.extend(result['sequences'])
        self.qualities.extend(result['qualities'])
        self.lengths = np.concatenate([self.lengths, np.array(result['lengths'], dtype=int)])
        self.gc_counts = np.concatenate([self.gc_counts, np.array(result['gc_counts'], dtype=int)])
        self.q20_counts = np.concatenate([self.q20_counts, np.array(result['q20_counts'], dtype=int)])
        self.q30_counts = np.concatenate([self.q30_counts, np.array(result['q30_counts'], dtype=int)])
        self.total_bases += np.sum(result['lengths'])



    def number_of_reads(self):
        return len(self.sequences)

    def total_bases_sequenced(self):
        return self.total_bases

    def q20_q30_scores(self):
        q20_percentage = (self.q20_counts.sum() / self.total_bases) * 100
        q30_percentage = (self.q30_counts.sum() / self.total_bases) * 100
        return q20_percentage, q30_percentage

    def gc_content(self):
        if self.total_bases > 0:
            return (self.gc_counts.sum() / self.total_bases) * 100
        else:
            return 0  # or some other appropriate value indicating no GC content could be calculated

    def read_lengths_statistics(self):
        if len(self.lengths) == 0:
            return {}
        return {
            'min_length': np.min(self.lengths),
            'max_length': np.max(self.lengths),
            'mean_length': np.mean(self.lengths),
            'median_length': np.median(self.lengths)
        }

    def quality_statistics(self):
        all_quality_scores = np.concatenate(self.qualities)
        return {
            'min_quality': np.min(all_quality_scores),
            'max_quality': np.max(all_quality_scores),
            'mean_quality': np.mean(all_quality_scores),
        }

    def qualities_vs_lengths(self):
        avg_qualities = [np.mean(quals) for quals in self.qualities]
        return {
            'read_lengths': self.lengths.tolist(),
            'avg_qualities': avg_qualities
        }

    def gc_content_per_sequence(self):
        gc_contents = (self.gc_counts.sum() / self.lengths) * 100
        return gc_contents.tolist()
