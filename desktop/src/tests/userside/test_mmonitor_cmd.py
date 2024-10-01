import gzip
import sys
import os
import logging
import argparse

# Add the project root directory to the Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, project_root)

# Add the src directory to the Python path
src_dir = os.path.abspath(os.path.join(project_root, 'src'))
sys.path.insert(0, src_dir)

# Ensure that the 'mmonitor' module can be imported
try:
    import mmonitor
except ImportError:
    raise ImportError("Failed to import 'mmonitor'. Please check the Python path and module structure.")

import unittest
import shutil
import subprocess
from datetime import date
from build_mmonitor_pyinstaller import ROOT
from src.mmonitor.userside.MMonitorCMD import MMonitorCMD

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

class TestMMonitorCMD(unittest.TestCase):
    def setUp(self):
        logger.debug("Setting up TestMMonitorCMD")
        self.test_data_dir = os.path.join(ROOT, "src", "tests", "test_data")
        self.output_dir = os.path.join(ROOT, "src", "resources", "pipeline_out")
        self.db_config = os.path.join(ROOT, "src", "resources", "db_config.json")
        self.centrifuge_db = os.path.join(ROOT, "src", "resources", "centrifuge_db", "p_compressed")
        
        # Create a small test FASTQ file
        os.makedirs(self.test_data_dir, exist_ok=True)
        self.test_fastq = os.path.join(self.test_data_dir, "test_sample.fastq")
        with open(self.test_fastq, 'w') as f:
            f.write("@read1\nACGT\n+\n!!!!\n@read2\nTGCA\n+\n!!!!\n")
        
        # Create a dummy config file if it doesn't exist
        if not os.path.exists(self.db_config):
            with open(self.db_config, 'w') as f:
                f.write('{"dummy": "config"}')

    def tearDown(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        if os.path.exists(self.test_fastq):
            os.remove(self.test_fastq)

    def test_emu_analysis_cmd(self):
        logger.debug("Starting test_emu_analysis_cmd")
        cmd = [
            sys.executable,
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", "taxonomy-16s",
            "-c", self.db_config,
            "-i", self.test_fastq,
            "-s", "test_sample",
            "-d", date.today().strftime("%Y-%m-%d"),
            "-p", "test_project",
            "-u", "test_subproject",
            "--overwrite"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, f"Command failed with error: {result.stderr}")
        self.assertIn("Emu analysis completed successfully", result.stdout)
        logger.debug("Finished test_emu_analysis_cmd")

    def test_centrifuge_analysis_cmd(self):
        logger.debug("Starting test_centrifuge_analysis_cmd")
        sample_name = "test_centrifuge_sample"
        input_dir = os.path.join(self.test_data_dir, sample_name)
        os.makedirs(input_dir, exist_ok=True)
        
        # Create a dummy fastq file in the input directory
        dummy_file = os.path.join(input_dir, f"{sample_name}.fastq.gz")
        with gzip.open(dummy_file, 'wt') as f:
            f.write("@read1\nACGT\n+\n!!!!\n@read2\nTGCA\n+\n!!!!\n")
        
        cmd = [
            sys.executable,  # Use the current Python interpreter
            "-m", "mmonitor.userside.MMonitorCMD",
            "-a", "taxonomy-wgs",
            "-c", self.db_config,
            "-i", input_dir,
            "-s", sample_name,
            "-d", date.today().strftime("%Y-%m-%d"),
            "-p", "test_project",
            "-u", "test_subproject",
            "--overwrite",
            "--centrifuge-db", self.centrifuge_db
        ]
        
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{project_root}:{src_dir}:{env.get('PYTHONPATH', '')}"
        
        
        result = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if "command not found" in result.stderr:
            self.skipTest("Centrifuge is not installed")
        self.assertEqual(result.returncode, 0, f"Command failed with error: {result.stderr}")
        self.assertIn("Centrifuge analysis completed successfully", result.stdout)
        logger.debug("Finished test_centrifuge_analysis_cmd")

        # Clean up the created directory
        shutil.rmtree(input_dir)

    # Add more tests for other analysis types as needed

class TestFunctionalPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        logger.debug("Setting up TestFunctionalPipeline")
        cls.test_data_path = os.path.join(ROOT, "src", "tests", "test_data", "test_functional_sample.fastq.gz")
        cls.output_dir = os.path.join(ROOT, "src", "resources", "pipeline_out")
        cls.config_path = os.path.join(ROOT, "src", "resources", "db_config.json")
        
        # Create a dummy config file if it doesn't exist
        if not os.path.exists(cls.config_path):
            with open(cls.config_path, 'w') as f:
                f.write('{"dummy": "config"}')
        
        # Create a dummy test file if it doesn't exist
        if not os.path.exists(cls.test_data_path):
            os.makedirs(os.path.dirname(cls.test_data_path), exist_ok=True)
            with gzip.open(cls.test_data_path, 'wt') as f:
                f.write("@read1\nACGT\n+\n!!!!\n@read2\nTGCA\n+\n!!!!\n")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.output_dir, ignore_errors=True)

    def test_functional_pipeline(self):
        logger.debug("Starting test_functional_pipeline")
        cmd = MMonitorCMD()
        args = cmd.parse_arguments([
            "-a", "functional",
            "-c", self.config_path,
            "-i", self.test_data_path,
            "-s", "test_sample",
            "-p", "test_project",
            "-u", "test_subproject",
            "-d", "2023-05-01",
            "-t", "12",
            "--overwrite"
        ])
        cmd.initialize_from_args(args)
        
        try:
            cmd.run()
        except subprocess.CalledProcessError as e:
            if "minimap2 is not installed" in str(e):
                self.skipTest("minimap2 is not installed")
            else:
                self.fail(f"Functional pipeline raised an unexpected exception: {e}")
        logger.debug("Finished test_functional_pipeline")

        # Check if output files exist
        expected_files = [
            "flye_out/assembly.fasta",
            "medaka_out/consensus.fasta",
            "metabat2_bins",
            "checkm2_out/quality_report.tsv",
            "bakta_out",
            "gtdbtk_out/classify/gtdbtk.bac120.summary.tsv"
        ]

        for file in expected_files:
            file_path = os.path.join(self.output_dir, "test_sample", file)
            self.assertTrue(os.path.exists(file_path), f"Expected output file not found: {file_path}")

    def test_assembly_pipeline(self):
        logger.debug("Starting test_assembly_pipeline")
        cmd = MMonitorCMD()
        args = cmd.parse_arguments([
            "-a", "assembly",
            "-c", self.config_path,
            "-i", self.test_data_path,
            "-s", "test_sample",
            "-p", "test_project",
            "-u", "test_subproject",
            "-d", "2023-05-01",
            "-t", "1",
            "--overwrite"
        ])
        cmd.initialize_from_args(args)
        
        try:
            cmd.run()
        except subprocess.CalledProcessError as e:
            if "minimap2 is not installed" in str(e):
                self.skipTest("minimap2 is not installed")
            else:
                self.fail(f"Assembly pipeline raised an unexpected exception: {e}")

        # Check if output files exist
        expected_files = [
            "flye_out/assembly.fasta",
            "medaka_out/consensus.fasta",
            "metabat2_bins"
        ]

        for file in expected_files:
            file_path = os.path.join(self.output_dir, "test_sample", file)
            self.assertTrue(os.path.exists(file_path), f"Expected output file not found: {file_path}")

    def test_invalid_input_file(self):
        cmd = MMonitorCMD()
        with self.assertRaises(SystemExit):
            cmd.parse_arguments([
                "-a", "functional",
                "-c", self.config_path,
                "-i", "/path/to/nonexistent/file.fastq.gz",
                "-s", "test_sample",
                "-p", "test_project",
                "-u", "test_subproject",
                "-d", "2023-05-01"
            ])

    def test_missing_required_args(self):
        cmd = MMonitorCMD()
        with self.assertRaises(SystemExit):
            cmd.parse_arguments([
                "-a", "functional",
                "-c", self.config_path
            ])

if __name__ == '__main__':
    unittest.main()