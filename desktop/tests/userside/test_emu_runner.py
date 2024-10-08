import unittest
import os
import shutil
from datetime import date
import sys
from pathlib import Path

# Add the root directory to sys.path
sys.path.append(str(Path(__file__).resolve().parents[3]))
project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))

from mmonitor.userside.emu_runner import EmuRunner
from mmonitor.database.django_db_interface import DjangoDBInterface
from build_mmonitor_pyinstaller import ROOT

class TestEmuRunner(unittest.TestCase):
    def setUp(self):
        self.emu_runner = EmuRunner()
        self.test_file = os.path.join(ROOT, "src", "tests", "userside", "test_sample.fastq.gz")
        self.db_config = os.path.join(ROOT, "src", "resources", "db_config.json")
        self.django_db = DjangoDBInterface(self.db_config)
        
        # Default values for testing
        self.sample_name = "test_sample"
        self.project_name = "test_project"
        self.subproject_name = "test_subproject"
        self.sample_date = date.today().strftime("%Y-%m-%d")
        
        # Create a temporary output directory
        self.output_dir = os.path.join(ROOT, "src", "resources", "pipeline_out", self.sample_name)
        os.makedirs(self.output_dir, exist_ok=True)

    def tearDown(self):
        # Clean up the temporary output directory
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_emu_analysis(self):
        # Run Emu analysis
        emu_success = self.emu_runner.run_emu([self.test_file], self.sample_name, min_abundance=0.01)
        self.assertTrue(emu_success, "Emu analysis failed")

        # Check if output files were created
        abundance_file = os.path.join(self.output_dir, f"{self.sample_name}_rel-abundance.tsv")
        self.assertTrue(os.path.exists(abundance_file), "Abundance file was not created")

        # Update database with results
        try:
            update_success = self.django_db.update_django_with_emu_out(
                self.output_dir, "species", self.sample_name, self.project_name,
                self.sample_date, self.subproject_name, True
            )
            self.assertTrue(update_success, "Failed to update database with Emu results")
        except Exception as e:
            self.fail(f"Database update raised an exception: {str(e)}")

        # You can add more assertions here to check the content of the output files
        # or verify the database entries if needed

if __name__ == '__main__':
    unittest.main()