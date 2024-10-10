import unittest
import os
import tempfile
import shutil
import sys
from unittest.mock import Mock, patch
import requests
import json

# Add the path to your project's root directory and the src directory
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, project_root)
sys.path.insert(0, src_path)

# Mock the ROOT constant with a string value
mock_root = project_root
sys.modules['build_mmonitor_pyinstaller'] = Mock(ROOT=mock_root)

# Mock the entire customtkinter module
sys.modules['customtkinter'] = Mock()

from mmonitor.userside.FolderWatcherWindow import FolderWatcherWindow
from mmonitor.userside.EmuRunner import EmuRunner
from mmonitor.userside.MMonitorCMD import MMonitorCMD

class TestEmuPipeline(unittest.TestCase):
    @patch('mmonitor.userside.FolderWatcherWindow.FolderWatcherWindow')
    def setUp(self, MockFolderWatcherWindow):
        self.test_folder = "/Users/timo/paper_mmonitor/watch/op0319_seq_on_0408/20240408_1735_MN35270_FAV63441_ce845c06/fastq_pass/barcode45"
        self.temp_output_dir = tempfile.mkdtemp()
        
        # Mock GUI reference
        self.mock_gui = Mock()
        self.mock_gui.db_path = self.temp_output_dir
        
        # Use the mocked FolderWatcherWindow
        self.folder_watcher = MockFolderWatcherWindow.return_value
        
        # Set up necessary attributes on the mocked FolderWatcherWindow
        self.folder_watcher.project_entry = Mock(get=lambda: 'test_project')
        self.folder_watcher.subproject_entry = Mock(get=lambda: 'test_subproject')
        self.folder_watcher.selected_date = "2024-04-08"
        self.folder_watcher.analysis_type = Mock(get=lambda: 'taxonomy-16s')
        self.folder_watcher.samples = {}

        # Create a real EmuRunner instance
        with patch('mmonitor.userside.EmuRunner.ROOT', new=mock_root):
            self.emu_runner = EmuRunner()

        # Set up Django server URL
        self.django_server_url = "https://www.mmonitor.org"  # Replace with your actual Django server URL

    def tearDown(self):
        shutil.rmtree(self.temp_output_dir)

    @patch('mmonitor.userside.FolderWatcherWindow.ROOT', new=mock_root)
    @patch('mmonitor.userside.EmuRunner.ROOT', new=mock_root)
    @patch('mmonitor.userside.MMonitorCMD.ROOT', new=mock_root)
    @patch('subprocess.run')
    @patch('requests.post')
    def test_emu_pipeline(self, mock_post, mock_subprocess_run):
        # Mock subprocess.run to simulate Emu running
        mock_subprocess_run.return_value = Mock(returncode=0, stdout="Emu analysis completed successfully")

        # Mock requests.post to simulate uploading to Django server
        mock_post.return_value = Mock(status_code=200, json=lambda: {"message": "Data uploaded successfully"})

        # Add the test folder to the samples
        sample_name = "barcode45"
        self.folder_watcher.samples[sample_name] = [os.path.join(self.test_folder, f) for f in os.listdir(self.test_folder) if f.endswith(('.fastq', '.fastq.gz'))]

        # Ensure there are sample files
        self.assertTrue(self.folder_watcher.samples[sample_name], "No sample files found")

        # Debug print
        print(f"EmuRunner custom_db_path: {self.emu_runner.custom_db_path}")

        # If custom_db_path is None, set it to a default value
        if self.emu_runner.custom_db_path is None:
            self.emu_runner.custom_db_path = os.path.join(mock_root, "src", "resources", "emu_db")
        
        # Verify that the custom_db_path exists and contains the necessary files
        self.assertTrue(os.path.exists(self.emu_runner.custom_db_path), f"Custom DB path does not exist: {self.emu_runner.custom_db_path}")
        self.assertTrue(os.path.exists(os.path.join(self.emu_runner.custom_db_path, "taxonomy.tsv")), "taxonomy.tsv not found in custom DB path")
        self.assertTrue(os.path.exists(os.path.join(self.emu_runner.custom_db_path, "species_taxid.fasta")), "species_taxid.fasta not found in custom DB path")

        # Run the Emu analysis
        input_file = self.folder_watcher.samples[sample_name][0]  # Use the first file for simplicity
        success = self.emu_runner.run_emu(
            input_file=input_file,
            output_dir=self.temp_output_dir,
            db_dir=self.emu_runner.custom_db_path,
            threads=12
        )

        # Assert that the Emu analysis was successful
        self.assertTrue(success, "Emu analysis should have been successful")

        # Assert that subprocess.run was called with the correct arguments
        mock_subprocess_run.assert_called_once()
        args, kwargs = mock_subprocess_run.call_args
        cmd = args[0]
        
        # Debug print
        print("Full Emu command:", ' '.join(cmd))

        # Check for essential components in the command
        self.assertIn('emu.py', cmd[1], "Emu script not found in command")
        self.assertIn('abundance', cmd[2], "Abundance command not found")
        self.assertIn(self.emu_runner.custom_db_path, cmd, f"Custom DB path {self.emu_runner.custom_db_path} not found in command")
        self.assertIn(self.temp_output_dir, cmd, "Output directory not found in command")
        self.assertIn(input_file, cmd, "Input file not found in command")

        # Check for specific parameters
        self.assertIn('--type', cmd, "--type parameter missing")
        self.assertIn('map-ont', cmd, "map-ont value missing for --type")
        self.assertIn('--min-abundance', cmd, "--min-abundance parameter missing")
        self.assertIn('--keep-files', cmd, "--keep-files parameter missing")
        self.assertIn('--N', cmd, "--N parameter missing")
        self.assertIn('--K', cmd, "--K parameter missing")

        # Create mock output files
        expected_files = [
            "emu_abundance.tsv",
            "emu_species_abundance.tsv",
            "emu_genus_abundance.tsv",
            "emu_family_abundance.tsv",
        ]
        for file in expected_files:
            file_path = os.path.join(self.temp_output_dir, file)
            with open(file_path, 'w') as f:
                f.write("Mock Emu output")
            self.assertTrue(os.path.exists(file_path), f"Expected output file {file} not found")

        # Simulate data processing
        processed_data = self.process_emu_results(self.temp_output_dir)

        # Upload results to Django server
        upload_url = f"{self.django_server_url}/api/upload_emu_results/"
        response = requests.post(upload_url, json=processed_data)

        # Assert that the upload was successful
        self.assertEqual(response.status_code, 200, "Failed to upload results to Django server")
        self.assertEqual(response.json()["message"], "Data uploaded successfully")

    def process_emu_results(self, output_dir):
        # Implement this method to process Emu output files and prepare data for upload
        processed_data = {
            "sample_name": "test_sample",
            "project": "test_project",
            "subproject": "test_subproject",
            "date": "2024-04-08",
            "results": {}
        }

        for file in ["emu_abundance.tsv", "emu_species_abundance.tsv", "emu_genus_abundance.tsv", "emu_family_abundance.tsv"]:
            file_path = os.path.join(output_dir, file)
            with open(file_path, 'r') as f:
                # Read and process the file content
                # This is a simplified example; you'll need to adapt it to your actual file format
                processed_data["results"][file] = f.read()

        return processed_data

if __name__ == '__main__':
    unittest.main()