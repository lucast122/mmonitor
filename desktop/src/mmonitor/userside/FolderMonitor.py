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
        if not event.is_directory and event.src_path.endswith(('.fastq', '.fq', '.fastq.gz')) and "_concatenated" not in event.src_path:
            print(f"New file detected: {event.src_path}")
            self.folder_monitor.add_new_file(event.src_path)

class FolderMonitor:
     def __init__(self, folders, watcher_window):
        self.folders = folders
        self.watcher_window = watcher_window
        self.observer = Observer()
        self.running = False

     def start(self):
        self.running = True
        for folder in self.folders:
            event_handler = SequencerFileHandler(self)
            self.observer.schedule(event_handler, folder, recursive=True)
        self.observer_thread = threading.Thread(target=self.observer.start, daemon=True)
        self.observer_thread.start()

     def stop(self):
        self.running = False
        self.observer.stop()
        self.observer_thread.join()

     def add_new_file(self, file_path):
        if self.running:
            self.watcher_window.add_new_file(file_path)