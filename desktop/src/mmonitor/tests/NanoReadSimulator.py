import random
import gzip
import os
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class NanoPoreReadSimulator:
    def __init__(self, read_length, error_rate, num_reads):
        self.read_length = read_length
        self.error_rate = error_rate
        self.num_reads = num_reads

    def simulate_reads(self):
        reads = []
        for i in range(self.num_reads):
            seq = ''.join(random.choice('ATCG') for _ in range(self.read_length))
            seq = self.introduce_errors(seq)
            quality = [random.randint(0, 40) for _ in range(len(seq))]
            record = SeqRecord(Seq(seq), id=f"read_{i}", description="", letter_annotations={"phred_quality": quality})
            reads.append(record)
        return reads

    def introduce_errors(self, seq):
        seq_list = list(seq)
        i = 0
        while i < len(seq_list):
            if random.random() < self.error_rate:
                error_type = random.choice(['substitution', 'insertion', 'deletion'])
                if error_type == 'substitution':
                    seq_list[i] = random.choice('ATCG')
                elif error_type == 'insertion':
                    seq_list.insert(i, random.choice('ATCG'))
                    i += 1
                elif error_type == 'deletion' and len(seq_list) > 1:
                    seq_list.pop(i)
                    i -= 1
            i += 1
        return ''.join(seq_list)

    def save_reads_to_file(self, filename):
        reads = self.simulate_reads()
        with gzip.open(filename, 'wt') as handle:
            SeqIO.write(reads, handle, 'fastq')

    def generate_files(self, output_dir, interval=30):
        file_count = 0
        while True:
            filename = os.path.join(output_dir, f"simulated_reads_{file_count}.fastq.gz")
            self.save_reads_to_file(filename)
            file_count += 1
            time.sleep(interval)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Simulate NanoPore reads and save to a specified directory.")
    parser.add_argument("output_dir", help="Directory to save simulated read files")
    parser.add_argument("--read_length", type=int, default=1000, help="Length of simulated reads")
    parser.add_argument("--error_rate", type=float, default=0.01, help="Error rate for simulated reads")
    parser.add_argument("--num_reads", type=int, default=100, help="Number of reads per file")
    parser.add_argument("--interval", type=int, default=30, help="Interval between file generation in seconds")

    args = parser.parse_args()

    simulator = NanoPoreReadSimulator(args.read_length, args.error_rate, args.num_reads)
    
    print(f"Generating simulated reads in {args.output_dir}. Press Ctrl+C to stop.")
    try:
        simulator.generate_files(args.output_dir, args.interval)
    except KeyboardInterrupt:
        print("\nStopped generating reads.")