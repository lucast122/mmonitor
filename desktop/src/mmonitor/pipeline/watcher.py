#!/usr/bin/env python3
"""
Folder watcher for real-time sequencing analysis.

Monitors a directory for new FASTQ files, groups them by sample (barcode),
concatenates when enough files accumulate, and triggers Snakemake analysis.

Designed to be used by both the CLI (mmonitor watch) and the GUI.
"""
import os
import re
import time
import gzip
import shutil
import logging
import threading
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import Callable, Optional

from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

logger = logging.getLogger(__name__)

FASTQ_EXTENSIONS = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}


def _is_fastq(path: str) -> bool:
    """Check if a file is a FASTQ file (not a concatenated one we created)."""
    p = Path(path)
    if '_concatenated' in p.name:
        return False
    return ''.join(p.suffixes[-2:]) in FASTQ_EXTENSIONS or p.suffix in FASTQ_EXTENSIONS


def infer_sample_name(filepath: str) -> str:
    """Infer sample name from file path.

    Looks for barcode/sample directory names in the path.
    Falls back to parent directory name.
    """
    parts = Path(filepath).parts
    for part in reversed(parts):
        if re.match(r'barcode\d+', part, re.IGNORECASE):
            return part
        if re.match(r'sample[_-]?\d+', part, re.IGNORECASE):
            return part
    # Fallback: use parent directory name
    return Path(filepath).parent.name


def preload_index(db_path: str):
    """Preload a database index into OS page cache.

    Reads all .cfr files into /dev/null so Linux keeps them in RAM.
    Subsequent centrifuger invocations will find the data already cached.
    """
    p = Path(db_path)
    if p.is_dir():
        cfr_files = list(p.rglob('*.cfr'))
    else:
        # Prefix: look for prefix.*.cfr
        cfr_files = list(Path(p).parent.glob(f'{p.name}.*.cfr'))

    if not cfr_files:
        logger.info("No .cfr index files found to preload")
        return

    total_size = sum(f.stat().st_size for f in cfr_files)
    logger.info(f"Preloading {len(cfr_files)} index files ({total_size / 1e9:.1f} GB) into page cache...")

    for f in cfr_files:
        with open(f, 'rb') as fh:
            while fh.read(1024 * 1024 * 64):  # 64 MB chunks
                pass

    logger.info("Index preloaded into page cache")


def concatenate_fastq(files: list, output_path: str):
    """Concatenate FASTQ files (handles both .gz and plain)."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(output, 'wt') if str(output).endswith('.gz') else open(output, 'w') as out_fh:
        for fpath in sorted(files):
            open_fn = gzip.open if fpath.endswith('.gz') else open
            with open_fn(fpath, 'rt') as in_fh:
                shutil.copyfileobj(in_fh, out_fh)

    logger.info(f"Concatenated {len(files)} files -> {output}")


class SampleTracker:
    """Tracks files per sample and decides when to trigger analysis."""

    def __init__(self, min_files: int = 5):
        self.min_files = min_files
        self._lock = threading.Lock()
        # sample_name -> {'files': set(), 'analyzed_count': int, 'last_analysis': datetime}
        self._samples = defaultdict(lambda: {
            'files': set(),
            'analyzed_count': 0,
            'last_analysis': None,
        })

    def add_file(self, filepath: str, sample_name: str = None):
        """Register a new file. Returns sample_name."""
        if sample_name is None:
            sample_name = infer_sample_name(filepath)
        with self._lock:
            self._samples[sample_name]['files'].add(filepath)
        return sample_name

    def get_pending_samples(self) -> dict:
        """Get samples that have enough new files for analysis.

        Returns dict of {sample_name: [file_paths]} for samples ready to analyze.
        """
        ready = {}
        with self._lock:
            for name, info in self._samples.items():
                total = len(info['files'])
                new_since_last = total - info['analyzed_count']
                if info['analyzed_count'] == 0:
                    # First analysis: need min_files total
                    if total >= self.min_files:
                        ready[name] = sorted(info['files'])
                else:
                    # Subsequent: need min_files new files
                    if new_since_last >= self.min_files:
                        ready[name] = sorted(info['files'])
        return ready

    def mark_analyzed(self, sample_name: str, file_count: int):
        """Mark a sample as analyzed with the number of files that were included.

        Args:
            sample_name: The sample that was analyzed.
            file_count: Number of files that were actually included in this
                        analysis run (NOT the current total — files may have
                        arrived during the run).
        """
        with self._lock:
            info = self._samples[sample_name]
            info['analyzed_count'] = file_count
            info['last_analysis'] = datetime.now()

    @property
    def samples(self) -> dict:
        """Get a snapshot of all tracked samples."""
        with self._lock:
            return {
                name: {
                    'total_files': len(info['files']),
                    'analyzed_count': info['analyzed_count'],
                    'new_files': len(info['files']) - info['analyzed_count'],
                    'last_analysis': info['last_analysis'],
                }
                for name, info in self._samples.items()
            }


class FastqWatchHandler(FileSystemEventHandler):
    """Watchdog handler that tracks new FASTQ files."""

    def __init__(self, tracker: SampleTracker,
                 on_new_file: Optional[Callable] = None):
        super().__init__()
        self.tracker = tracker
        self.on_new_file = on_new_file

    def on_created(self, event):
        if event.is_directory:
            return
        if _is_fastq(event.src_path):
            sample = self.tracker.add_file(event.src_path)
            logger.debug(f"New file: {event.src_path} -> sample {sample}")
            if self.on_new_file:
                self.on_new_file(event.src_path, sample)


class SequencingWatcher:
    """High-level folder watcher that triggers Snakemake analysis.

    Monitors a folder for FASTQ files, groups by sample, and periodically
    triggers analysis via a callback.

    Args:
        watch_dir: Directory to watch (recursively)
        analysis_callback: Called with (samples_dict, concat_dir) when analysis should run.
                          samples_dict = {sample_name: {"fastq": [concat_path]}}
        min_files: Minimum new files before triggering analysis per sample
        interval: Seconds between analysis checks
        concat_dir: Directory to write concatenated files
        on_status: Optional callback(message) for status updates
    """

    def __init__(self, watch_dir: str, analysis_callback: Callable,
                 min_files: int = 5, interval: int = 300,
                 concat_dir: str = None,
                 on_status: Optional[Callable[[str], None]] = None):
        self.watch_dir = Path(watch_dir)
        self.analysis_callback = analysis_callback
        self.interval = interval
        self.concat_dir = Path(concat_dir) if concat_dir else self.watch_dir / '_mmonitor_concat'
        self.on_status = on_status or (lambda msg: None)

        self.tracker = SampleTracker(min_files=min_files)
        self.handler = FastqWatchHandler(
            self.tracker,
            on_new_file=lambda path, sample: self.on_status(
                f"New file: {Path(path).name} -> {sample}"
            )
        )
        self.observer = None
        self._stop_event = threading.Event()
        self._analysis_thread = None

    def scan_existing(self):
        """Scan for existing FASTQ files in the watch directory."""
        count = 0
        for root, dirs, files in os.walk(self.watch_dir):
            for fname in files:
                fpath = os.path.join(root, fname)
                if _is_fastq(fpath):
                    self.tracker.add_file(fpath)
                    count += 1
        if count:
            self.on_status(f"Found {count} existing FASTQ files")
            status = self.tracker.samples
            for name, info in status.items():
                self.on_status(f"  {name}: {info['total_files']} files")

    def start(self):
        """Start watching and periodic analysis."""
        self.scan_existing()

        self.observer = Observer()
        self.observer.schedule(self.handler, str(self.watch_dir), recursive=True)
        self.observer.start()
        self.on_status(f"Watching {self.watch_dir} (interval: {self.interval}s)")

        self._stop_event.clear()
        self._analysis_thread = threading.Thread(target=self._analysis_loop, daemon=True)
        self._analysis_thread.start()

    def stop(self):
        """Stop watching."""
        self._stop_event.set()
        if self.observer:
            self.observer.stop()
            self.observer.join()
            self.observer = None
        if self._analysis_thread:
            self._analysis_thread.join(timeout=5)
        self.on_status("Watcher stopped")

    def _analysis_loop(self):
        """Periodically check for pending samples and trigger analysis.

        Checks immediately on start (for existing files), then every interval.
        """
        while not self._stop_event.is_set():
            self._check_and_analyze()
            self._stop_event.wait(self.interval)
            if self._stop_event.is_set():
                break

    def check_now(self):
        """Manually trigger an analysis check (useful for GUI)."""
        self._check_and_analyze()

    def _check_and_analyze(self):
        """Check for pending samples, concatenate, and trigger analysis."""
        pending = self.tracker.get_pending_samples()
        if not pending:
            self.on_status("No samples ready for analysis")
            return

        self.on_status(f"Triggering analysis for {len(pending)} sample(s): {', '.join(pending.keys())}")

        # Concatenate files per sample
        samples_dict = {}
        for sample_name, files in pending.items():
            concat_path = str(self.concat_dir / sample_name / f"{sample_name}_concatenated.fastq.gz")
            self.on_status(f"Concatenating {len(files)} files for {sample_name}...")
            concatenate_fastq(files, concat_path)
            samples_dict[sample_name] = {'fastq': [concat_path]}

        # Record how many files were included in this analysis
        files_in_analysis = {name: len(files) for name, files in pending.items()}

        # Trigger analysis callback
        try:
            self.analysis_callback(samples_dict, str(self.concat_dir))
            # Mark as analyzed with the count of files that were actually included
            for sample_name, count in files_in_analysis.items():
                self.tracker.mark_analyzed(sample_name, count)
                self.on_status(f"Analysis complete for {sample_name} ({count} files)")
        except Exception as e:
            self.on_status(f"Analysis failed: {e}")
            logger.error(f"Analysis callback failed: {e}", exc_info=True)
