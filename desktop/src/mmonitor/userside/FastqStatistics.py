import gzip
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from Bio import SeqIO
import time
import json
import traceback

def process_fastq_file(file_path):
    try:
        sequences_temp = []
        qualities_temp = []
        gc_counts = []
        lengths = []
        
        is_gzipped = file_path.endswith(".gz")
        with (gzip.open(file_path, 'rt') if is_gzipped else open(file_path, 'r')) as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                qual = record.letter_annotations["phred_quality"]
                
                sequences_temp.append(seq)
                qualities_temp.append(qual)
                gc_counts.append(seq.count('G') + seq.count('C'))
                lengths.append(len(seq))
        
        return {
            'sequences': sequences_temp,
            'qualities': qualities_temp,
            'gc_counts': gc_counts,
            'lengths': lengths
        }
        
    except Exception as e:
        print(f"Error processing FASTQ file: {str(e)}")
        traceback.print_exc()
        return None


class FastqStatistics:

    def __init__(self, file_paths, multi=True, num_threads=4):  
        self.file_paths = file_paths if isinstance(file_paths, list) else [file_paths]
        # Initialize lists and numpy arrays for aggregated data
        self.sequences = []
        self.qualities = []
        self.lengths = np.array([], dtype=int)
        self.gc_counts = np.array([], dtype=int)
        self.q20_counts = np.array([], dtype=int)
        self.q30_counts = np.array([], dtype=int)
        self.total_bases = 0
        
        # Process files
        if multi and len(self.file_paths) > 1:
            with ProcessPoolExecutor(max_workers=num_threads) as executor:
                results = list(executor.map(process_fastq_file, self.file_paths))
        else:
            results = [process_fastq_file(fp) for fp in self.file_paths]
            
        # Aggregate results
        for result in results:
            self.sequences.extend(result['sequences'])
            self.qualities.extend(result['qualities'])
            self.lengths = np.append(self.lengths, result['lengths'])
            self.gc_counts = np.append(self.gc_counts, result['gc_counts'])
            for quality_scores in result['qualities']:
                self.q20_counts = np.append(self.q20_counts, np.sum(np.array(quality_scores) >= 20))
                self.q30_counts = np.append(self.q30_counts, np.sum(np.array(quality_scores) >= 30))
            
        # Calculate total bases
        self.total_bases = np.sum(self.lengths)


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
            return 0  

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
        all_quality_scores = np.concatenate([np.array(quals) for quals in self.qualities])
        return {
            'min_quality': np.min(all_quality_scores),
            'max_quality': np.max(all_quality_scores),
            'mean_quality': np.mean(all_quality_scores),
        }

    def qualities_vs_lengths(self):
        avg_qualities = [np.mean(np.array(quals)) for quals in self.qualities]
        return {
            'read_lengths': self.lengths.tolist(),
            'avg_qualities': avg_qualities
        }

    def gc_content_per_sequence(self):
        gc_contents = (self.gc_counts / self.lengths) * 100
        return gc_contents.tolist()

    def _convert_to_json_serializable(self, obj):
        """Convert numpy types to Python native types for JSON serialization."""
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8, np.int16, np.int32, np.int64,
                          np.uint8, np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: self._convert_to_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._convert_to_json_serializable(x) for x in obj]
        return obj

    def calculate_plot_data(self):
        """Pre-calculate all data needed for plots to avoid recalculation."""
        plot_data = {
            'version': 1,  # Increment this when plot data format changes
            'read_length_dist': {
                'data': self.read_lengths_statistics(),
                'histogram': np.histogram(self.lengths, bins=50),
            },
            'quality_dist': {
                'mean': float(np.mean([np.mean(np.array(q)) for q in self.qualities])),
                'histogram': np.histogram([np.mean(np.array(q)) for q in self.qualities], bins=40, range=(0, 40)),
            },
            'gc_content': {
                'mean': float(self.gc_content()),
                'histogram': np.histogram([(gc/l)*100 for gc, l in zip(self.gc_counts, self.lengths)], 
                                       bins=50, range=(0, 100)),
            },
            'quality_summary': {
                'q20_score': float((self.q20_counts.sum() / self.total_bases) * 100),
                'q30_score': float((self.q30_counts.sum() / self.total_bases) * 100),
            }
        }
        
        # Convert numpy arrays to lists for JSON serialization
        for key in plot_data:
            if isinstance(plot_data[key], dict) and 'histogram' in plot_data[key]:
                hist_counts, hist_bins = plot_data[key]['histogram']
                plot_data[key]['histogram'] = {
                    'counts': hist_counts.tolist(),
                    'bins': hist_bins.tolist()
                }
        
        # Convert numpy types to Python native types before JSON serialization
        for key, value in plot_data.items():
            plot_data[key] = self._convert_to_json_serializable(value)
        
        return plot_data

    def process_fastq_file(self, file_path):
        """Process a FASTQ file and calculate statistics."""
        try:
            sequences_temp = []
            qualities_temp = []
            gc_counts = []
            lengths = []
            
            is_gzipped = file_path.endswith(".gz")
            with (gzip.open(file_path, 'rt') if is_gzipped else open(file_path, 'r')) as f:
                for record in SeqIO.parse(f, "fastq"):
                    seq = str(record.seq)
                    qual = record.letter_annotations["phred_quality"]
                    
                    sequences_temp.append(seq)
                    qualities_temp.append(qual)
                    gc_counts.append(seq.count('G') + seq.count('C'))
                    lengths.append(len(seq))
            
            self.sequences = np.array(sequences_temp)
            self.qualities = np.array(qualities_temp)
            self.gc_counts = np.array(gc_counts)
            self.lengths = np.array(lengths)
            
            # Calculate plot data
            plot_data = self.calculate_plot_data()
            
            # Convert numpy types to native Python types
            plot_data = self._convert_to_json_serializable(plot_data)
            
            return plot_data
            
        except Exception as e:
            print(f"Error processing FASTQ file: {str(e)}")
            traceback.print_exc()
            return None
