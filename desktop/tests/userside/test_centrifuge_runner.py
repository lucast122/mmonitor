import unittest
import os
import shutil
from datetime import date
from mmonitor.userside.CentrifugeRunner import CentrifugeRunner
from mmonitor.database.django_db_interface import DjangoDBInterface
from build_mmonitor_pyinstaller import ROOT

class TestCentrifugeRunner(unittest.TestCase):
    def setUp(self):
        self.centrifuge_runner = CentrifugeRunner()
        self.test_data_dir = os.path.join(ROOT, "src", "tests", "test_data")
        self.db_config = os.path.join(ROOT, "src", "resources", "db_config.json")
        self.django_db = DjangoDBInterface(self.db_config)
        
        # Default values for testing
        self.sample_name = "test_centrifuge_sample"
        self.project_name = "test_project"
        self.subproject_name = "test_subproject"
        self.sample_date = date.today().strftime("%Y-%m-%d")
        
        # Create a temporary input directory with a dummy file
        self.input_dir = os.path.join(self.test_data_dir, self.sample_name)
        os.makedirs(self.input_dir, exist_ok=True)
        dummy_file = os.path.join(self.input_dir, f"{self.sample_name}.fastq.gz")
        with open(dummy_file, 'w') as f:
            f.write("Dummy fastq content")
        
        # Create a temporary output directory
        self.output_dir = os.path.join(ROOT, "src", "resources", "pipeline_out", self.sample_name)
        os.makedirs(self.output_dir, exist_ok=True)

    def tearDown(self):
        # Clean up the temporary directories
        if os.path.exists(self.input_dir):
            shutil.rmtree(self.input_dir)
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

    def test_centrifuge_analysis(self):
        # Run Centrifuge analysis
        centrifuge_success = self.centrifuge_runner.run_centrifuge(self.input_dir, self.sample_name)
        self.assertTrue(centrifuge_success, "Centrifuge analysis failed")

        # Check if output files were created
        output_file = os.path.join(self.output_dir, f"{self.sample_name}_centrifuge_out.tsv")
        self.assertTrue(os.path.exists(output_file), "Centrifuge output file was not created")

        # Update database with results
        try:
            update_success = self.django_db.send_nanopore_record_centrifuge(
                output_file, self.sample_name, self.project_name,
                self.subproject_name, self.sample_date, True
            )
            self.assertTrue(update_success, "Failed to update database with Centrifuge results")
        except Exception as e:
            self.fail(f"Database update raised an exception: {str(e)}")

        # You can add more assertions here to check the content of the output files
        # or verify the database entries if needed

if __name__ == '__main__':
    unittest.main()