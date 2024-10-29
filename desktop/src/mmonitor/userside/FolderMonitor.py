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
        if not event.is_directory and event.src_path.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')) and "_concatenated" not in event.src_path:
            print(f"SequencerFileHandler: New file detected: {event.src_path}")
            self.folder_monitor.add_new_file(event.src_path)

class FolderMonitor:
    def __init__(self, folders, watcher_window):
        self.folders = folders
        self.watcher_window = watcher_window
        self.observer = Observer()
        self.handler = SequencerFileHandler(self)
        print(f"FolderMonitor: Initialized with folders: {self.folders}")

    def start(self):
        for folder in self.folders:
            # Changed recursive to True to monitor subdirectories
            self.observer.schedule(self.handler, folder, recursive=True)
            print(f"FolderMonitor: Scheduled observer for folder: {folder} (including subdirectories)")
        self.observer.start()
        print("FolderMonitor: Observer started")

    def stop(self):
        print("FolderMonitor: Stopping observer")
        self.observer.stop()
        self.observer.join()
        print("FolderMonitor: Observer stopped")

    def add_new_file(self, file_path):
        print(f"FolderMonitor: Adding new file: {file_path}")
        self.watcher_window.add_new_file(file_path)

class FileHandler(FileSystemEventHandler):
    def __init__(self, folder_monitor):
        self.folder_monitor = folder_monitor

    def on_created(self, event):
        if not event.is_directory and (event.src_path.endswith('.fastq') or event.src_path.endswith('.fastq.gz')):
            self.folder_monitor.add_new_file(event.src_path)
