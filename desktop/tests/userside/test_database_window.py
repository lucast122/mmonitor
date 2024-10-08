import unittest
from unittest.mock import patch, MagicMock, call, mock_open
import os
import sys
import tempfile
import shutil
from pathlib import Path
import multiprocessing
import subprocess

# Add the project root to the Python path
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

from mmonitor.userside.DatabaseWindow import DatabaseWindow

class TestDatabaseWindow(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.db_window = DatabaseWindow(None)
        self.db_window.base_dir = self.temp_dir
        self.db_window.emu_db_path = os.path.join(self.temp_dir, "emu_db")
        self.db_window.centrifuge_db_path = os.path.join(project_root, "src", "resources", "centrifuge_db")
        self.test_input = """
>NR_181927.1 Kaistella faecalis strain F4 16S ribosomal RNA, complete sequence kraken:taxid|123456
ACAATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGGGAGGCCTAACACATGCAAG
>NR_181928.1 Another species 16S ribosomal RNA, partial sequence taxid=789012
ACAATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGGGAGGCCTAACACATGCAAG
>NR_181929.1 Yet another species 16S ribosomal RNA
ACAATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGGGAGGCCTAACACATGCAAG
"""
        self.mock_db_window = MagicMock()
        self.mock_db_window.log_progress = MagicMock()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch('multiprocessing.Pool')
    @patch('builtins.open', new_callable=mock_open, read_data="Dummy content for testing")
    def test_create_seq2taxid_map(self, mock_file, mock_pool):
        mock_pool.return_value.__enter__.return_value.map.return_value = [
            (["NR_181927.1\t123456\n"], 0),
            (["NR_181928.1\t789012\n"], 0),
            ([], 1)
        ]

        mock_db_window = MagicMock()
        skipped, total = DatabaseWindow.create_seq2taxid_map(mock_db_window, "input.fna", "output.map")

        self.assertEqual(skipped, 1)
        self.assertEqual(total, 3)

        mock_file.assert_any_call("input.fna", "r")
        mock_file.assert_any_call("output.map", "w")

    @patch('subprocess.run')
    @patch('subprocess.Popen')
    @patch('multiprocessing.Pool')
    @patch('os.path.getsize')
    @patch('builtins.open', new_callable=mock_open, read_data="Dummy content for testing")
    @patch('os.path.exists')
    @patch('tempfile.TemporaryDirectory')
    def test_build_emu_db(self, mock_temp_dir, mock_exists, mock_file, mock_getsize, mock_pool, mock_popen, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.base_dir = "/path/to/base"
        mock_db_window.emu_db_path = "/path/to/emu_db"
        mock_db_window.create_seq2taxid_map.return_value = (0, 100)  # 0 skipped, 100 processed
        mock_getsize.return_value = 1000  # Non-zero file size
        mock_exists.return_value = True  # Simulate that all files exist
        mock_temp_dir.return_value.__enter__.return_value = "/tmp/dir"

        mock_process = MagicMock()
        mock_process.stdout.readline.side_effect = ["Output line 1", "Output line 2", ""]
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        DatabaseWindow._build_emu_db(mock_db_window, ["bacteria", "archaea"])

        # Check if necessary methods were called
        mock_db_window.log_progress.assert_called()
        mock_run.assert_called()
        mock_popen.assert_called()

        # Check if the emu_db_path was updated
        self.assertIn("custom_emu_db", mock_db_window.emu_db_path)

    @patch('subprocess.run')
    @patch('os.path.exists')
    @patch('os.path.getsize')
    @patch('builtins.open', new_callable=mock_open, read_data="Dummy content for testing")
    @patch('tempfile.TemporaryDirectory')
    def test_build_centrifuge_db(self, mock_temp_dir, mock_file, mock_getsize, mock_exists, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.centrifuge_db_path = os.path.join(project_root, "src", "resources", "centrifuge_db")
        mock_exists.return_value = True  # Simulate that all files exist
        mock_getsize.return_value = 1000  # Non-zero file size
        mock_temp_dir.return_value.__enter__.return_value = "/tmp/dir"

        DatabaseWindow._build_centrifuge_db(mock_db_window, ["bacteria", "archaea"])

        # Check if necessary methods were called
        mock_run.assert_called()

        # Check if the progress window was closed
        mock_db_window.close_progress_window.assert_called()

    @patch('subprocess.run')
    @patch('tempfile.TemporaryDirectory')
    @patch('builtins.open', new_callable=mock_open, read_data="Test content")
    @patch('os.path.exists')
    def test_download_files(self, mock_exists, mock_file, mock_temp_dir, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.base_dir = "/path/to/base"
        mock_db_window.emu_db_path = "/path/to/emu_db"
        mock_temp_dir.return_value.__enter__.return_value = "/tmp/dir"
        mock_exists.return_value = True  # Simulate that all files exist

        DatabaseWindow._build_emu_db(mock_db_window, ["bacteria", "archaea"])

        # Check if wget was called to download the necessary files
        mock_run.assert_any_call(["wget", "-O", "/tmp/dir/taxdump.tar.gz", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"], 
                                 capture_output=True, text=True, check=True)
        mock_run.assert_any_call(["wget", "-O", "/tmp/dir/bacteria_16S.fna.gz", "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"], 
                                 capture_output=True, text=True, check=True)

        # Check if the files were extracted
        mock_run.assert_any_call(["tar", "-xzvf", "/tmp/dir/taxdump.tar.gz", "-C", unittest.mock.ANY], 
                                 capture_output=True, text=True, check=True)
        mock_run.assert_any_call(["gunzip", "-f", "/tmp/dir/bacteria_16S.fna.gz"], 
                                 capture_output=True, text=True, check=True)

    @patch('os.path.getsize')
    @patch('builtins.open', new_callable=mock_open, read_data="")
    def test_empty_seq2taxid_map(self, mock_file, mock_getsize):
        mock_db_window = MagicMock()
        mock_getsize.return_value = 0  # Simulate empty file

        with self.assertRaises(ValueError) as context:
            DatabaseWindow._build_emu_db(mock_db_window, ["bacteria"])

        self.assertIn("seq2taxid map is empty", str(context.exception))

    @patch('subprocess.run')
    @patch('subprocess.Popen')
    @patch('multiprocessing.cpu_count')
    @patch('builtins.open', new_callable=mock_open, read_data="Test content")
    @patch('os.path.exists')
    def test_emu_build_command(self, mock_exists, mock_file, mock_cpu_count, mock_popen, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.base_dir = "/path/to/base"
        mock_db_window.emu_db_path = "src/resources/emu_db"
        mock_db_window.create_seq2taxid_map.return_value = (0, 100)  # 0 skipped, 100 processed
        mock_cpu_count.return_value = 4
        mock_exists.return_value = True  # Simulate that all files exist

        mock_process = MagicMock()
        mock_process.stdout.readline.side_effect = ["Output line 1", "Output line 2", ""]
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        DatabaseWindow._build_emu_db(mock_db_window, ['bacteria'])

        # Check if the emu build-database command was called with the correct parameters
        expected_command = [
            sys.executable, unittest.mock.ANY, "build-database",
            "--sequences", unittest.mock.ANY,
            "--seq2tax", unittest.mock.ANY,
            "--ncbi-taxonomy", unittest.mock.ANY,
            unittest.mock.ANY
        ]
        mock_popen.assert_called_with(expected_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    def test_create_domain_selection_dialog(self):
        # This test simulates user input for domain selection
        with patch('mmonitor.userside.DatabaseWindow.ctk.CTkToplevel') as mock_toplevel:
            mock_dialog = MagicMock()
            mock_toplevel.return_value = mock_dialog
            
            def side_effect(dialog):
                self.db_window.selected_domains = ["bacteria", "archaea"]
                dialog.destroy()
            
            mock_dialog.wait_window.side_effect = side_effect
            
            result = self.db_window.create_domain_selection_dialog()
            
            self.assertEqual(result, ["bacteria", "archaea"])

    @patch('mmonitor.userside.DatabaseWindow.urllib.request.urlopen')
    @patch('mmonitor.userside.DatabaseWindow.tarfile.open')
    @patch('builtins.open', new_callable=mock_open, read_data="1\t| Bacteria\t|\t| scientific name\t|\n")
    def test_download_and_process_taxdump(self, mock_file, mock_tarfile, mock_urlopen):
        mock_urlopen.return_value.__enter__.return_value.read.return_value = b"mock data"
        mock_tarfile.return_value.__enter__.return_value = MagicMock()
        
        self.db_window.download_and_process_taxdump()
        
        self.assertIn("Bacteria", self.db_window.species_taxid_map)
        self.assertEqual(self.db_window.species_taxid_map["Bacteria"], "1")

if __name__ == '__main__':
    unittest.main()