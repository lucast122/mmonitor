import os
import tempfile
import unittest
from unittest import mock
from unittest.mock import patch, MagicMock, call
import multiprocessing
import subprocess
from mmonitor.userside.DatabaseWindow import DatabaseWindow, process_chunk

class TestDatabaseWindow(unittest.TestCase):
    def setUp(self):
        self.test_input = """
>NR_181927.1 Kaistella faecalis strain F4 16S ribosomal RNA, complete sequence kraken:taxid|123456
ACAATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGGGAGGCCTAACACATGCAAG
CCGAGCGGTAGAAGATCTTCGGATCTTTGAGAGCGGCGTACGGGTGCGTAACACGTGTGCA
>NR_181928.1 Another species 16S ribosomal RNA, partial sequence taxid=789012
ACAATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGGGAGGCCTAACACATGCAAG
CCGAGCGGTAGAAGATCTTCGGATCTTTGAGAGCGGCGTACGGGTGCGTAACACGTGTGCA
>NR_181929.1 Yet another species 16S ribosomal RNA
ACAATGGAGAGTTTGATCCTGGCTCAGGATGAACGCTAGCGGGAGGCCTAACACATGCAAG
CCGAGCGGTAGAAGATCTTCGGATCTTTGAGAGCGGCGTACGGGTGCGTAACACGTGTGCA
"""
        self.mock_db_window = MagicMock()
        self.mock_db_window.log_progress = MagicMock()

    def test_process_chunk(self):
        chunk = self.test_input.strip().split('\n')
        result = process_chunk(chunk)
        self.assertEqual(len(result), 2)
        self.assertIn("NR_181927.1\t123456\n", result)
        self.assertIn("NR_181928.1\t789012\n", result)

    @patch('multiprocessing.Pool')
    def test_create_seq2taxid_map(self, mock_pool):
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_input:
            temp_input.write(self.test_input)
            temp_input.flush()
            input_filename = temp_input.name

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_output:
            output_filename = temp_output.name

        mock_db_window = MagicMock()
        DatabaseWindow.create_seq2taxid_map(mock_db_window, input_filename, output_filename)

        with open(output_filename, 'r') as f:
            content = f.read()
            self.assertIn("NR_181927.1\t123456", content)
            self.assertIn("NR_181928.1\t789012", content)

        # Clean up temporary files
        os.unlink(input_filename)
        os.unlink(output_filename)

    @patch('subprocess.run')
    def test_build_centrifuge_db(self, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.centrifuge_db_path = "/path/to/centrifuge_db"

        DatabaseWindow._build_centrifuge_db(mock_db_window, ['bacteria'], 'ncbi', '/path/to/index')

        mock_run.assert_called()

    @patch('subprocess.run')
    @patch('subprocess.Popen')
    @patch('multiprocessing.cpu_count')
    def test_build_emu_db(self, mock_cpu_count, mock_popen, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.base_dir = "/path/to/base"
        mock_db_window.emu_db_path = "/path/to/emu_db"
        mock_db_window.create_seq2taxid_map.return_value = None
        mock_cpu_count.return_value = 4

        DatabaseWindow._build_emu_db(mock_db_window, ['bacteria'])

        mock_run.assert_called()

    @patch('subprocess.run')
    @patch('tempfile.TemporaryDirectory')
    def test_download_files(self, mock_temp_dir, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.base_dir = "/path/to/base"
        mock_db_window.emu_db_path = "/path/to/emu_db"
        mock_temp_dir.return_value.__enter__.return_value = "/tmp/dir"

        DatabaseWindow._build_emu_db(mock_db_window, ['bacteria'])

        mock_run.assert_any_call(["wget", "-O", "/tmp/dir/taxdump.tar.gz", "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"], 
                                 capture_output=True, text=True, check=True)

    @patch('os.path.getsize')
    @patch('builtins.open', new_callable=mock.mock_open, read_data="")
    def test_empty_seq2taxid_map(self, mock_open, mock_getsize):
        mock_db_window = MagicMock()
        mock_getsize.return_value = 0  # Simulate empty file

        with self.assertRaises(ValueError):
            DatabaseWindow._build_emu_db(mock_db_window, ['bacteria'])

    @patch('subprocess.run')
    @patch('subprocess.Popen')
    @patch('multiprocessing.cpu_count')
    def test_emu_build_command(self, mock_cpu_count, mock_popen, mock_run):
        mock_db_window = MagicMock()
        mock_db_window.base_dir = "/path/to/base"
        mock_db_window.emu_db_path = "src/resources/emu_db"
        mock_db_window.create_seq2taxid_map.return_value = None
        mock_cpu_count.return_value = 4

        mock_process = MagicMock()
        mock_process.stdout.readline.side_effect = ["Output line 1", "Output line 2", ""]
        mock_process.poll.return_value = 0
        mock_popen.return_value = mock_process

        DatabaseWindow._build_emu_db(mock_db_window, ['bacteria'])

        expected_command = [
            "emu", "build-database", "emu_custom_db",
            "--sequences", mock.ANY,
            "--seq2tax", mock.ANY,
            "--taxonomy", mock.ANY,
            "--output", mock.ANY,
            "--threads", "12"
        ]
        mock_popen.assert_called_with(expected_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

if __name__ == '__main__':
    unittest.main()