import gzip
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from Bio import SeqIO
import time
import json
import traceback
import os

def process_fastq_file(file_path):
    """Process a single FASTQ file and return statistics or None if errors occur."""
    try:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            return None
            
        sequences_temp = []
        qualities_temp = []
        gc_counts = []
        lengths = []
        
        is_gzipped = file_path.endswith(".gz")
        
        try:
            with (gzip.open(file_path, 'rt') if is_gzipped else open(file_path, 'r')) as f:
                try:
                    for record in SeqIO.parse(f, "fastq"):
                        try:
                            seq = str(record.seq)
                            qual = record.letter_annotations["phred_quality"]
                            
                            # Validate sequence and quality have the same length
                            if len(seq) != len(qual):
                                print(f"Warning: Sequence and quality length mismatch in {file_path}, record {record.id} - skipping this record")
                                continue
                                
                            sequences_temp.append(seq)
                            qualities_temp.append(qual)
                            gc_counts.append(seq.count('G') + seq.count('C'))
                            lengths.append(len(seq))
                        except Exception as e:
                            print(f"Error processing record in {file_path}: {str(e)}")
                            continue
                except Exception as parse_error:
                    print(f"Error parsing FASTQ file {file_path}: {str(parse_error)}")
                    # If we have partial data, return it instead of None
                    if sequences_temp:
                        print(f"Returning partial data ({len(sequences_temp)} sequences) from {file_path}")
                    else:
                        return None
        except Exception as file_error:
            print(f"Error opening file {file_path}: {str(file_error)}")
            return None
        
        # Only return data if we successfully processed at least one sequence
        if not sequences_temp:
            print(f"No valid sequences found in {file_path}")
            return None
            
        return {
            'sequences': sequences_temp,
            'qualities': qualities_temp,
            'gc_counts': gc_counts,
            'lengths': lengths
        }
        
    except Exception as e:
        print(f"Unexpected error processing FASTQ file {file_path}: {str(e)}")
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
        
        # Filter out None results (failed files)
        valid_results = [r for r in results if r is not None]
        
        if not valid_results:
            print("Warning: All FASTQ files failed to process. Statistics will be empty.")
            return
            
        # Aggregate results from valid files
        for result in valid_results:
            try:
                self.sequences.extend(result['sequences'])
                self.qualities.extend(result['qualities'])
                self.lengths = np.append(self.lengths, result['lengths'])
                self.gc_counts = np.append(self.gc_counts, result['gc_counts'])
                
                for quality_scores in result['qualities']:
                    self.q20_counts = np.append(self.q20_counts, np.sum(np.array(quality_scores) >= 20))
                    self.q30_counts = np.append(self.q30_counts, np.sum(np.array(quality_scores) >= 30))
            except (KeyError, TypeError) as e:
                print(f"Error aggregating result: {str(e)}")
                continue
            
        # Calculate total bases
        self.total_bases = np.sum(self.lengths) if len(self.lengths) > 0 else 0


    def number_of_reads(self):
        return len(self.sequences)

    def total_bases_sequenced(self):
        return self.total_bases

    def q20_q30_scores(self):
        if self.total_bases > 0:
            q20_percentage = (self.q20_counts.sum() / self.total_bases) * 100
            q30_percentage = (self.q30_counts.sum() / self.total_bases) * 100
            return q20_percentage, q30_percentage
        return 0, 0

    def gc_content(self):
        if self.total_bases > 0:
            return (self.gc_counts.sum() / self.total_bases) * 100
        return 0  

    def read_lengths_statistics(self):
        if len(self.lengths) == 0:
            return {
                'min_length': 0,
                'max_length': 0,
                'mean_length': 0,
                'median_length': 0
            }
        return {
            'min_length': np.min(self.lengths),
            'max_length': np.max(self.lengths),
            'mean_length': np.mean(self.lengths),
            'median_length': np.median(self.lengths)
        }

    def quality_statistics(self):
        if not self.qualities:
            return {
                'min_quality': 0,
                'max_quality': 0,
                'mean_quality': 0,
            }
        try:
            all_quality_scores = np.concatenate([np.array(quals) for quals in self.qualities])
            return {
                'min_quality': np.min(all_quality_scores) if len(all_quality_scores) > 0 else 0,
                'max_quality': np.max(all_quality_scores) if len(all_quality_scores) > 0 else 0,
                'mean_quality': np.mean(all_quality_scores) if len(all_quality_scores) > 0 else 0,
            }
        except Exception as e:
            print(f"Error calculating quality statistics: {str(e)}")
            return {
                'min_quality': 0,
                'max_quality': 0,
                'mean_quality': 0,
            }

    def qualities_vs_lengths(self):
        if not self.qualities or not self.lengths.size:
            return {
                'read_lengths': [],
                'avg_qualities': []
            }
        try:
            avg_qualities = [np.mean(np.array(quals)) for quals in self.qualities]
            return {
                'read_lengths': self.lengths.tolist(),
                'avg_qualities': avg_qualities
            }
        except Exception as e:
            print(f"Error calculating qualities vs lengths: {str(e)}")
            return {
                'read_lengths': [],
                'avg_qualities': []
            }

    def gc_content_per_sequence(self):
        if not self.lengths.size or not self.gc_counts.size:
            return []
        try:
            gc_contents = (self.gc_counts / self.lengths) * 100
            return gc_contents.tolist()
        except Exception as e:
            print(f"Error calculating GC content per sequence: {str(e)}")
            return []

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
        try:
            if not self.lengths.size or not self.qualities:
                # Return empty data structure if no valid sequences
                return {
                    'version': 1,
                    'read_length_dist': {'data': self.read_lengths_statistics(), 'histogram': {'counts': [], 'bins': []}},
                    'quality_dist': {'mean': 0, 'histogram': {'counts': [], 'bins': []}},
                    'gc_content': {'mean': 0, 'histogram': {'counts': [], 'bins': []}},
                    'quality_summary': {'q20_score': 0, 'q30_score': 0}
                }
            
            # Generate histogram data safely
            length_hist = np.histogram(self.lengths, bins=min(50, len(self.lengths)) if len(self.lengths) > 0 else 1)
            
            # Safely calculate quality histogram
            quality_means = []
            for q in self.qualities:
                try:
                    quality_means.append(np.mean(np.array(q)))
                except Exception:
                    continue
            
            quality_hist = np.histogram(quality_means, bins=min(40, len(quality_means)) if len(quality_means) > 0 else 1, 
                                     range=(0, 40)) if quality_means else ([],[])
            
            # Safely calculate GC content histogram
            gc_percentages = []
            for gc, l in zip(self.gc_counts, self.lengths):
                try:
                    gc_percentages.append((gc/l)*100)
                except Exception:
                    continue
            
            gc_hist = np.histogram(gc_percentages, bins=min(50, len(gc_percentages)) if len(gc_percentages) > 0 else 1, 
                                range=(0, 100)) if gc_percentages else ([],[])
                            
            plot_data = {
                'version': 1,  # Increment this when plot data format changes
                'read_length_dist': {
                    'data': self.read_lengths_statistics(),
                    'histogram': length_hist,
                },
                'quality_dist': {
                    'mean': float(np.mean(quality_means)) if quality_means else 0,
                    'histogram': quality_hist,
                },
                'gc_content': {
                    'mean': float(self.gc_content()),
                    'histogram': gc_hist,
                },
                'quality_summary': {
                    'q20_score': float((self.q20_counts.sum() / self.total_bases) * 100) if self.total_bases > 0 else 0,
                    'q30_score': float((self.q30_counts.sum() / self.total_bases) * 100) if self.total_bases > 0 else 0,
                }
            }
            
            # Convert numpy arrays to lists for JSON serialization
            for key in plot_data:
                if isinstance(plot_data[key], dict) and 'histogram' in plot_data[key]:
                    hist_counts, hist_bins = plot_data[key]['histogram']
                    plot_data[key]['histogram'] = {
                        'counts': hist_counts.tolist() if hasattr(hist_counts, 'tolist') else [],
                        'bins': hist_bins.tolist() if hasattr(hist_bins, 'tolist') else []
                    }
            
            # Convert numpy types to Python native types before JSON serialization
            for key, value in plot_data.items():
                plot_data[key] = self._convert_to_json_serializable(value)
            
            return plot_data
            
        except Exception as e:
            print(f"Error generating plot data: {str(e)}")
            traceback.print_exc()
            # Return an empty but valid data structure
            return {
                'version': 1,
                'read_length_dist': {'data': {}, 'histogram': {'counts': [], 'bins': []}},
                'quality_dist': {'mean': 0, 'histogram': {'counts': [], 'bins': []}},
                'gc_content': {'mean': 0, 'histogram': {'counts': [], 'bins': []}},
                'quality_summary': {'q20_score': 0, 'q30_score': 0}
            }

    def process_fastq_file(self, file_path):
        """Process a FASTQ file and calculate statistics."""
        return process_fastq_file(file_path)
