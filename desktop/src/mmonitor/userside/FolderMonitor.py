import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import threading
import re
import os
from collections import defaultdict
import json

class SequencerFileHandler(FileSystemEventHandler):
    def __init__(self, folder_monitor):
        super().__init__()
        self.folder_monitor = folder_monitor

    def on_created(self, event):
        if not event.is_directory and event.src_path.endswith(('.fastq', '.fq', '.fastq.gz')):
            print(f"New file detected: {event.src_path}")
            self.folder_monitor.add_new_file(event.src_path, is_new=True)

class FolderMonitor:
    def __init__(self, folders_to_watch, gui_reference):
        self.folders_to_watch = folders_to_watch
        self.gui_reference = gui_reference
        self.observer = Observer()
        self.files_by_sample = defaultdict(list)
        self.analysis_threshold = 10  # Number of files to trigger analysis
        self.lock = threading.Lock()
        self.known_files = self.load_known_files()

    def start(self):
        event_handler = SequencerFileHandler(self)
        for folder in self.folders_to_watch:
            self.observer.schedule(event_handler, folder, recursive=True)
            self.scan_existing_files(folder)
        self.observer.start()

    def stop(self):
        self.observer.stop()
        self.observer.join()
        self.save_known_files()

    def scan_existing_files(self, folder):
        for root, _, files in os.walk(folder):
            for file in files:
                if file.endswith(('.fastq', '.fq', '.fastq.gz')):
                    file_path = os.path.join(root, file)
                    self.add_new_file(file_path, is_new=file_path not in self.known_files)

    def add_new_file(self, file_path, is_new=True):
        sample_name, barcode = self.extract_sample_info(file_path)
        with self.lock:
            self.files_by_sample[sample_name].append((file_path, is_new))
            self.gui_reference.add_new_file(file_path, sample_name, is_new)
            if is_new:
                self.known_files.add(file_path)
            if len([f for f, new in self.files_by_sample[sample_name] if new]) >= self.analysis_threshold:
                self.trigger_analysis(sample_name)

    def trigger_analysis(self, sample_name):
        print(f"Triggering analysis for sample: {sample_name}")
        # Call the appropriate analysis method in the GUI
        self.gui_reference.run_analysis_for_sample(sample_name)

    @staticmethod
    def extract_sample_info(file_path):
        path_parts = file_path.split(os.sep)
        sample_name = None
        barcode = None

        if 'fastq_pass' in path_parts:
            idx = path_parts.index('fastq_pass')
            if idx + 1 < len(path_parts):
                next_part = path_parts[idx + 1]
                if next_part.startswith('barcode'):
                    barcode = next_part
                    sample_name = f"Sample_{barcode}"
                else:
                    sample_name = next_part
        else:
            # Try to extract sample name from parent directory
            sample_name = os.path.basename(os.path.dirname(file_path))

        # If still no sample name, use a default
        if not sample_name:
            sample_name = f"Unknown_Sample_{hash(file_path)}"

        return sample_name, barcode

    def get_files_for_sample(self, sample_name):
        return self.files_by_sample.get(sample_name, [])

    def load_known_files(self):
        try:
            with open('known_files.json', 'r') as f:
                return set(json.load(f))
        except FileNotFoundError:
            return set()

    def save_known_files(self):
        with open('known_files.json', 'w') as f:
            json.dump(list(self.known_files), f)

