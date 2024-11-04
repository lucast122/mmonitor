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
    def __init__(self, folders, callback):
        self.folders = folders
        self.callback = callback
        self.observer = None
        self.active_watches = set()  # Track active watches
        print(f"FolderMonitor: Initialized with folders: {folders}")

    def start(self):
        """Start monitoring folders"""
        if self.observer is not None:
            print("FolderMonitor: Observer already running, stopping first")
            self.stop()
        
        self.observer = Observer()
        
        for folder in self.folders:
            if folder not in self.active_watches:
                try:
                    self.observer.schedule(
                        FileSystemEventHandler(self.callback),
                        folder,
                        recursive=True
                    )
                    self.active_watches.add(folder)
                    print(f"FolderMonitor: Scheduled observer for folder: {folder} (including subdirectories)")
                except Exception as e:
                    print(f"FolderMonitor: Error scheduling watch for {folder}: {e}")
        
        try:
            self.observer.start()
            print("FolderMonitor: Observer started")
        except Exception as e:
            print(f"FolderMonitor: Error starting observer: {e}")

    def stop(self):
        """Stop monitoring folders"""
        if self.observer is not None:
            print("FolderMonitor: Stopping observer")
            self.observer.stop()
            self.observer.join()
            self.observer = None
            self.active_watches.clear()
            print("FolderMonitor: Observer stopped")

    def __del__(self):
        """Ensure observer is stopped when object is destroyed"""
        self.stop()

    def add_new_file(self, file_path):
        print(f"FolderMonitor: Adding new file: {file_path}")
        self.callback.add_new_file(file_path)

class FileHandler(FileSystemEventHandler):
    def __init__(self, folder_monitor):
        self.folder_monitor = folder_monitor

    def on_created(self, event):
        if not event.is_directory and (event.src_path.endswith('.fastq') or event.src_path.endswith('.fastq.gz')):
            self.folder_monitor.add_new_file(event.src_path)
