import unittest
import tempfile
import os
import time
import tkinter as tk
from mmonitor.tests.NanoReadSimulator import NanoPoreReadSimulator
from mmonitor.userside.FolderWatcherWindow import FolderWatcherWindow

class MockGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.db_path = "mock_db_path"

class TestFolderWatcher(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.mock_gui = MockGUI()
        self.folder_watcher = FolderWatcherWindow(self.mock_gui)
        self.folder_watcher.folder_entry.insert(0, self.temp_dir)
        self.simulator = NanoPoreReadSimulator(read_length=100, error_rate=0.01, num_reads=10)

    def tearDown(self):
        self.folder_watcher.stop_watching()
        for file in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, file))
        os.rmdir(self.temp_dir)
        self.mock_gui.destroy()

    def test_folder_watcher_with_simulated_reads(self):
        self.folder_watcher.start_watching()

        for i in range(3):
            filename = os.path.join(self.temp_dir, f"simulated_reads_{i}.fastq.gz")
            self.simulator.save_reads_to_file(filename)
            self.folder_watcher.add_new_file(filename)
            time.sleep(0.1)
            self.mock_gui.update()

        self.assertEqual(len(self.folder_watcher.samples), 3)

        for sample in self.folder_watcher.samples.values():
            self.assertTrue(sample.get('analyzed', False))

    def test_auto_analyze_feature(self):
        self.folder_watcher.auto_analyze_var.set(True)
        self.folder_watcher.start_watching()

        filename = os.path.join(self.temp_dir, "test_auto_analyze.fastq.gz")
        self.simulator.save_reads_to_file(filename)
        self.folder_watcher.add_new_file(filename)
        time.sleep(0.1)
        self.mock_gui.update()

        sample_name = self.folder_watcher.suggest_sample_name(os.path.basename(filename))
        self.assertTrue(self.folder_watcher.samples[sample_name].get('analyzed', False))

if __name__ == '__main__':
    unittest.main()