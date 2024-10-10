import glob
import logging
import multiprocessing
import os.path
import subprocess
import zipfile
from os import path
import sys
import os

import pandas as pd
import requests

from build_mmonitor_pyinstaller import ROOT

"""
This is a runner for a functional analysis pipeline. As input it takes raw nanopore reads, then assembles them with
flye then corrects the assembly medaka, bins with metabat2, assigns taxonomy with gtdb-tk, annotatts MAGs with gtdb-tk
and finally runs KEGG pathway analysis on annotated MAGs using keggcharter
"""


class FunctionalRunner:


    def __init__(self):

        self.logger = logging.getLogger('timestamp')
        self.cpus = multiprocessing.cpu_count()  # get number of cpu cores available

        self.working_directory = os.getcwd()
        # get absolute path of all tools used and safe as class variable
        # os.chdir("../../lib/")
        self.flye_path = os.path.join(ROOT, "lib", "Flye-2.9.5", "bin", "flye")
        self.resources_path = os.path.join(ROOT,"src","resources")
        self.minimap2_path = os.path.join(ROOT, "lib", "minimap2", "minimap2")

        print(self.flye_path)
        print(self.minimap2_path)
        os.chdir(self.working_directory)
        # all kegg ids to input into keggcharter. TODO: put them in a file and parse instead of hardcoding here
        self.kegg_ids = [10, 20, 30, 40, 51, 52, 53, 61, 62, 71, 73, 100, 120, 121, 130, 140, 190, 195, 196, 220, 230,
                         232, 240, 250, 253, 254,
                         260, 261, 270, 280, 281, 290, 300, 310, 311, 330, 331, 332, 333, 340, 350, 360, 361, 362, 363,
                         364, 365, 380, 400,
                         401, 402, 403, 404, 405, 410, 430, 440, 450, 460, 470, 480, 500, 510, 511, 512, 513, 514, 515,
                         520, 521, 522, 523,
                         524, 525, 531, 532, 533, 534, 540, 541, 542, 550, 552, 561, 562, 563, 564, 565, 571, 572, 590,
                         591, 592, 600, 601,
                         603, 604, 620, 621, 622, 623, 624, 625, 626, 627, 630, 633, 640, 642, 643, 650, 660, 670, 680,
                         710, 720, 730, 740,
                         750, 760, 770, 780, 785, 790, 791, 830, 860, 900, 901, 902, 903, 904, 905, 906, 908, 909, 910,
                         920, 930, 940, 941, 942,
                         943, 944, 945, 950, 960, 965, 966, 970, 980, 981, 982, 983, 984, 996, 997, 998, 999, 1010,
                         1040, 1051, 1052, 1053,
                         1054, 1055, 1056, 1057, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1100, 1110, 1120, 1200, 1210,
                         1212, 1220, 1230, 1232, 1240, 1250, 1501, 1502, 1503, 1521, 1522, 1523, 1524, 2010, 2020, 2024,
                         2025, 2026,
                         2030, 2040, 2060, 3008, 3010, 3013, 3015, 3018, 3020, 3022, 3030, 3040, 3050, 3060, 3070, 3230,
                         3240, 3250, 3320, 3410, 3420, 3430, 3440, 3450, 3460, 4010, 4011, 4012, 4013, 4014, 4015, 4016,
                         4020, 4022, 4024, 4060, 4061, 4062, 4064, 4066, 4068, 4070, 4071, 4072, 4075, 4080, 4110, 4111,
                         4112, 4113, 4114, 4115, 4120, 4122, 4130, 4136, 4137, 4138, 4139, 4140, 4141, 4142, 4144, 4145,
                         4146, 4150, 4151, 4152, 4210, 4211, 4212, 4213, 4214, 4215, 4216, 4217, 4218, 4260, 4261, 4270,
                         4310, 4320, 4330, 4340, 4341, 4350, 4360, 4361, 4370, 4371, 4380, 4390, 4391, 4392, 4510, 4512,
                         4514, 4520, 4530, 4540, 4550, 4610, 4611, 4612, 4613, 4614, 4620, 4621, 4622, 4623, 4624, 4625,
                         4626, 4630, 4640, 4650, 4657, 4658, 4659, 4660, 4662, 4664, 4666, 4668, 4670, 4672, 4710, 4711,
                         4712, 4713, 4714, 4720, 4721, 4722, 4723, 4724, 4725, 4726, 4727, 4728, 4730, 4740, 4742, 4744,
                         4745, 4750, 4810, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920, 4921, 4922,
                         4923, 4924, 4925, 4926, 4927, 4928, 4929, 4930, 4931, 4932, 4933, 4934, 4935, 4936, 4940, 4950,
                         4960, 4961, 4962, 4964, 4966, 4970, 4971, 4972, 4973, 4974, 4975, 4976, 4977, 4978, 4979, 5010,
                         5012, 5014, 5016, 5017, 5020, 5022, 5030, 5031, 5032, 5033, 5034, 5100, 5110, 5111, 5120, 5130,
                         5131, 5132, 5133, 5134, 5135, 5140, 5142, 5143, 5144, 5145, 5146, 5150, 5152, 5160, 5161, 5162,
                         5163, 5164, 5165, 5166, 5167, 5168, 5169, 5170, 5171, 5200, 5202, 5203, 5204, 5205, 5206, 5207,
                         5208, 5210, 5211, 5212, 5213, 5214, 5215, 5216, 5217, 5218, 5219, 5220, 5221, 5222, 5223, 5224,
                         5225, 5226, 5230, 5231, 5235, 5310, 5320, 5321, 5322, 5323, 5330, 5332, 5340, 5410, 5412, 5414,
                         5415, 5416, 5417, 5418, 7011, 7012, 7013, 7014, 7015, 7016, 7017, 7018, 7019, 7020, 7021, 7023,
                         7024, 7025, 7026, 7027, 7028, 7029, 7030, 7031, 7032, 7033, 7034, 7035, 7036, 7037, 7038, 7039,
                         7040, 7041, 7042, 7043, 7044, 7045, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055,
                         7056, 7057, 7110, 7112, 7114, 7117, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220,
                         7221, 7222, 7223, 7224, 7225, 7226, 7227, 7228, 7229, 7230, 7231, 7232, 7233, 7234, 7235]
        self.kegg_ids = [str(int) for int in self.kegg_ids]
        self.kegg_ids = [s.zfill(5) for s in self.kegg_ids]

        self.kegg_ids = ",".join(self.kegg_ids)

        self.checkm2_db_path = os.path.expanduser("~/.checkm2db")
        self.gtdbtk_db_path = os.path.expanduser("~/.gtdbtk")

    def check_software_avail(self):
        try:
            subprocess.call([f'python3 {ROOT}/lib/Flye-2.9/bin/flye', '-h'], stdout=open(os.devnull, 'w'),
                            stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error("Flye executable not found.")

        try:
            subprocess.call(['racon', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error(
                "Make sure that racon is installed and on the sytem path. For more info visit https://github.com/lbcb-sci/racon")

        # try:
        #     subprocess.call(['prokka', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        # except FileNotFoundError:
        #     self.logger.error(
        #         "Prokka is not found in PATH. Try installing via conda")
        #     subprocess.call("conda install -c conda-forge -c bioconda -c defaults prokka")

        try:
            subprocess.call(['medaka_cocnsensus', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error(
                "Make sure that medaka is installed and on the sytem path. For more info visit https://github.com/nanoporetech/medaka"
                "If you have troubles installing on ubuntu try running ")
            self.logger.error("bzip2 g++ zlib1g-dev libbz2-dev liblzma-dev libffi-dev libncurses5-dev")
            self.logger.error("libcurl4-gnutls-dev libssl-dev curl make cmake wget python3-all-dev")
            self.logger.error("python-virtualenv")
            self.logger.error("first.")

        try:
            subprocess.call(['diamond', 'help'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error(
                "Make sure that racon is installed and on the system path. For more info visit https://github.com/bbuchfink/diamond")


    def run_flye(self, input_files, sample_name, output_dir, threads):
        flye_out = os.path.join(output_dir, "flye_out")
        os.makedirs(flye_out, exist_ok=True)

        env = os.environ.copy()
        env["PATH"] = f"{os.path.dirname(self.minimap2_path)}:{env['PATH']}"

        cmd = [self.flye_path, "--nano-raw"] + input_files + ["--out-dir", flye_out, "--threads", str(threads)]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
            print(f"Flye stdout: {result.stdout}")
            print(f"Flye stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Flye: {e}")
            print(f"Flye stderr: {e.stderr}")
            raise
        return os.path.join(flye_out, "assembly.fasta")


    def run_metabat2_pipeline(self, contig_file, read_file, dir_path):
        output_dir = os.path.join(dir_path, "metabat_bins")

        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        sam_files = []
        bam_files = []

        # Align reads to contigs using minimap2
        sam_file = os.path.join(dir_path, os.path.basename(read_file)) + ".sam"
        minimap2_cmd = [
            "minimap2",
            "-ax", "map-ont", contig_file,
            read_file,
            "-o", sam_file
        ]
        try:
            print(f"minimap2 cmd: {minimap2_cmd}")
            result = subprocess.run(minimap2_cmd, check=True, text=True, capture_output=True)
            print("Aligned reads to assembly", result.stdout)
            sam_files.append(sam_file)
            print(sam_files)
        except subprocess.CalledProcessError as e:
            print(f"Error executing minimap2 command for {read_file}:", e.stderr)
            return

        # convert to bam, sort and index resulting sam file
        for sam_file in sam_files:
            bam_file = sam_file.replace(".sam", ".bam")
            sorted_bam_file = bam_file.replace(".bam", ".sorted.bam")

            samtools_view_cmd = ["samtools", "view", "-bS", sam_file, "-o", bam_file]
            samtools_sort_cmd = ["samtools", "sort", bam_file, "-o", sorted_bam_file]
            samtools_index_cmd = ["samtools", "index", sorted_bam_file]

            try:
                subprocess.run(samtools_view_cmd, check=True, text=True, capture_output=True)
                subprocess.run(samtools_sort_cmd, check=True, text=True, capture_output=True)
                subprocess.run(samtools_index_cmd, check=True, text=True, capture_output=True)
                bam_files.append(sorted_bam_file)
            except subprocess.CalledProcessError as e:
                print(f"Error processing SAM/BAM files: {e.stderr}")
                return

        # Generate depth file

        depth_file = os.path.join(dir_path, "depth.txt")
        jgi_cmd = ["jgi_summarize_bam_contig_depths", "--outputDepth", depth_file] + bam_files
        try:
            subprocess.run(jgi_cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print("Error executing jgi_summarize_bam_contig_depths command:", e.stderr)
            return

        # Run Metabat2
        metabat2_cmd = [
            "metabat2",
            "-i", contig_file,
            "-a", depth_file,
            "-o", os.path.join(output_dir, "bin")
        ]
        try:
            result = subprocess.run(metabat2_cmd, check=True, text=True, capture_output=True)
            print("Metabat2 command executed successfully:", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error executing Metabat2 command:", e.stderr)

    def run_medaka_consensus(self, assembly_file, read_file, dir_path, model="r941_min_high_g303"):
        output_dir = os.path.join(dir_path, "medaka_corrected")

        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        consensus_output = os.path.join(output_dir, "consensus.hdf")
        print(self.cpus)

        # Run Medaka consensus
        medaka_cmd = [
            "medaka_consensus",
            "-i", read_file,
            "-d", assembly_file,
            "-o", output_dir,
            "-m", model,
            "-t", str(self.cpus)
        ]
        try:
            print(f"medaka cmd: {medaka_cmd}")
            result = subprocess.run(medaka_cmd, check=True, text=True, capture_output=True)
            print("Medaka consensus command executed successfully:", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error executing Medaka consensus command:", e.stderr)

    def run_fastcat(self, input_files, output_dir):
        output_file = os.path.join(output_dir, "concatenated.fastq")
        cmd = ["fastcat", "-o", output_file] + input_files
        subprocess.run(cmd, check=True)
        return output_file

    def run_filtlong(self, input_file, output_dir, min_length, min_quality):
        output_file = os.path.join(output_dir, "filtered.fastq")
        cmd = ["filtlong", "--min_length", str(min_length), "--min_mean_q", str(min_quality), input_file]
        with open(output_file, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)
        return output_file

    def run_metaflye(self, input_file, output_dir):
        cmd = ["flye", "--meta", "--nano-raw", input_file, "--out-dir", output_dir]
        subprocess.run(cmd, check=True)
        return os.path.join(output_dir, "assembly.fasta")

    def run_semibin(self, assembly_file, reads_file, output_dir):
        bins_dir = os.path.join(output_dir, "bins")
        cmd = ["SemiBin", "single_easy_bin", "-i", assembly_file, "-r", reads_file, "-o", bins_dir]
        subprocess.run(cmd, check=True)
        return bins_dir

    def run_gtdbtk(self, bins_dir, output_dir, threads):
        if not os.path.exists(self.gtdbtk_db_path):
            print("Downloading GTDB-TK database...")
            subprocess.run(["gtdbtk", "download", "--all"], check=True)
        
        gtdbtk_out = os.path.join(output_dir, "gtdbtk_out")
        cmd = ["gtdbtk", "classify_wf", "--genome_dir", bins_dir, "--out_dir", gtdbtk_out, "--cpus", str(threads)]
        subprocess.run(cmd, check=True)
        return gtdbtk_out

    def run_checkm2(self, bins_dir, output_dir, threads):
        if not os.path.exists(self.checkm2_db_path):
            print("Downloading CheckM2 database...")
            subprocess.run(["checkm2", "database", "--download"], check=True)
        
        checkm2_out = os.path.join(output_dir, "checkm2_out")
        cmd = ["checkm2", "predict", "--threads", str(threads), "--input", bins_dir, "--output-directory", checkm2_out]
        subprocess.run(cmd, check=True)
        return checkm2_out

    def run_bakta(self, bins_dir, output_dir):
        bakta_out = os.path.join(output_dir, "bakta_out")
        for bin_file in os.listdir(bins_dir):
            if bin_file.endswith(".fa") or bin_file.endswith(".fasta"):
                bin_path = os.path.join(bins_dir, bin_file)
                cmd = ["bakta", "--output", os.path.join(bakta_out, bin_file), bin_path]
                subprocess.run(cmd, check=True)
        return bakta_out

    def create_keggcharter_input(self, path_to_prokka_output):
        keggcharter_sheet = {'taxonomy': ['']}
        df = pd.DataFrame(keggcharter_sheet)
        df_list = []
        df.to_csv(path_to_prokka_output + "keggcharter.tsv", sep='\t')
        for tsv in glob.glob(f"{path_to_prokka_output}/*.tsv"):
            print(tsv)
            data = pd.read_csv(tsv, sep='\t')
            tax = tsv.split('.tsv')[0]
            tax = tax.split('/')[-1]
            try:
                tax = tax[:tax.rfind('_')] + tax[tax.rfind('.'):]
            except:
                print()
            tax = tax.replace('_', ' ')
            if ".fasta" in tax:
                tax = tax.removesuffix(".fasta")
            print(tax)
            data["taxonomy"] = tax
            df_list.append(data)
        df = pd.concat(df_list)
        df.to_csv(path_to_prokka_output + "keggcharter.tsv", sep='\t')

    def env_exists(self, env_name):
        try:
            envs = subprocess.check_output(['conda', 'env', 'list'], text=True)
            return any(env_name in line for line in envs.splitlines())
        except subprocess.CalledProcessError as e:
            print(f"Error checking conda environments: {e.stderr}")
            return False

    def create_conda_env_from_yaml(self, yaml_file="environment.yml"):
        yaml_path = os.path.join(self.resources_path, yaml_file)
        env_name = self.get_env_name_from_yaml(yaml_path)

        if self.env_exists(env_name):
            print(f"Environment '{env_name}' already exists. Updating the environment.")
            try:
                subprocess.run(['conda', 'env', 'update', '-f', yaml_path], check=True)
                print(f"Conda environment '{env_name}' updated successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error updating conda environment '{env_name}': {e.stderr}")
                raise
        else:
            print(f"Creating conda environment '{env_name}' from {yaml_file}.")
            try:
                subprocess.run(['conda', 'env', 'create', '-f', yaml_path], check=True)
                print(f"Conda environment created from {yaml_file} successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error creating conda environment from {yaml_file}: {e.stderr}")
                raise

    def get_env_name_from_yaml(self, yaml_path):
        with open(yaml_path, 'r') as f:
            for line in f:
                if line.startswith('name:'):
                    return line.split(':')[1].strip()
        raise ValueError(f"Could not find environment name in {yaml_path}")

    def run_checkm2(self, bins_dir, dir_path):
        output_dir = os.path.join(dir_path, "checkm2_results")

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        checkm2_cmd = [
            "checkm2", "predict",
            "-x", "fa",
            "-i", bins_dir,
            "-o", output_dir,
            "--threads", str(self.cpus),
            "--force"
        ]

        conda_prefix = subprocess.check_output(['conda', 'info', '--base'], text=True).strip()
        activate_env = os.path.join(conda_prefix, 'etc', 'profile.d', 'conda.sh')
        command = f'source {activate_env} && conda activate mmonitor && {" ".join(checkm2_cmd)}'

        try:
            print(f"checkm2 cmd: {command}")
            subprocess.run(command, shell=True, check=True, text=True, executable='/bin/bash', capture_output=True)
            print("CheckM2 command executed successfully.")
        except subprocess.CalledProcessError as e:
            print("Error executing CheckM2 command:", e.stderr)
            print("Attempting to download CheckM2 database...")
            try:
                download_cmd = f'source {activate_env} && conda activate checkm2 && checkm2 database --download'
                subprocess.run(download_cmd, shell=True, check=True, text=True, executable='/bin/bash', capture_output=True)
                print("CheckM2 database downloaded successfully. Retrying command...")
                subprocess.run(command, shell=True, check=True, text=True, executable='/bin/bash', capture_output=True)
                print("CheckM2 command executed successfully after downloading the database.")
            except subprocess.CalledProcessError as db_e:
                print("Error downloading CheckM2 database or retrying the command:", db_e.stderr)

    def run_bakta(self, input_file, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        bakta_cmd = [
            "bakta",
            "--db", self.resources_path,
            "--input", input_file,
            "--output", output_dir,
            "--threads", str(self.cpus)
        ]

        conda_prefix = subprocess.check_output(['conda', 'info', '--base'], text=True).strip()
        activate_env = os.path.join(conda_prefix, 'etc', 'profile.d', 'conda.sh')
        command = f'source {activate_env} && conda activate mmonitor && {" ".join(bakta_cmd)}'
        try:
            print(f"bakta cmd: {bakta_cmd}")
            result = subprocess.run(command, shell=True, check=True, text=True, executable='/bin/bash', capture_output=True)
            print("Bakta command executed successfully:", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error executing Bakta command:", e.stderr)
            print("Attempting to download Bakta database...")
            try:
                download_cmd = ["bakta_db", "download", "--output", self.resources_path, "--type full"]
                subprocess.run(download_cmd, check=True, text=True, capture_output=True)
                print("Bakta database downloaded successfully. Retrying command...")
                subprocess.run(bakta_cmd, check=True, text=True, capture_output=True)
                print("Bakta command executed successfully after downloading the database.")
            except subprocess.CalledProcessError as db_e:
                print("Error downloading CheckM2 database or retrying the command:", db_e.stderr)

    def run_gtdbtk(self, bins_dir, output_dir, threads):
        gtdbtk_out = os.path.join(output_dir, "gtdbtk_out")
        os.makedirs(gtdbtk_out, exist_ok=True)

        try:
            subprocess.run(["gtdbtk", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            print("GTDB-TK not found. Attempting to install using pip...")
            try:
                subprocess.run([sys.executable, "-m", "pip", "install", "gtdbtk"], check=True)
                print("GTDB-TK installed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Failed to install GTDB-TK: {e}")
                return None

        cmd = [
            "gtdbtk", "classify_wf",
            "--genome_dir", bins_dir,
            "--out_dir", gtdbtk_out,
            "--cpus", str(threads),
            "--extension", "fa"
        ]

        try:
            subprocess.run(cmd, check=True)
            print("GTDB-TK analysis completed successfully.")
            return gtdbtk_out
        except subprocess.CalledProcessError as e:
            print(f"Error running GTDB-TK: {e}")
            return None

    def run_keggcharter(self, kegg_out, keggcharter_input):
        # self.create_keggcharter_input(keggcharter_input)
        cmd = f"python {ROOT}/lib/KEGGCharter-0.3.4/keggcharter.py -o {kegg_out} -f {keggcharter_input} -tc taxonomy -ecc EC_number --input-quantification -mm {self.kegg_ids}"
        os.system(cmd)

    def run_flye(self, input_files, sample_name, output_dir, threads):
        flye_out = os.path.join(output_dir, "flye_out")
        os.makedirs(flye_out, exist_ok=True)

        # Check which minimap2 is being used
        try:
            minimap2_path = subprocess.check_output(["which", "minimap2"]).decode().strip()
            print(f"Using minimap2 from: {minimap2_path}")
        except subprocess.CalledProcessError:
            print("minimap2 not found in PATH")

        cmd = [self.flye_path, "--nano-raw"] + input_files + ["--out-dir", flye_out, "--threads", str(threads)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running Flye: {e}")
            print(f"Flye stderr: {e.stderr}")
            raise
        return os.path.join(flye_out, "assembly.fasta")

    def run_medaka(self, assembly, reads, output_dir, threads):
        medaka_out = os.path.join(output_dir, "medaka_out")
        os.makedirs(medaka_out, exist_ok=True)
        cmd = ["medaka_consensus", "-i", reads, "-d", assembly, "-o", medaka_out, "-t", str(threads)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running Medaka: {e}")
            print(f"Medaka stderr: {e.stderr}")
            raise
        return os.path.join(medaka_out, "consensus.fasta")

    def run_metabat2(self, assembly, reads, output_dir, threads):
        bins_dir = os.path.join(output_dir, "metabat2_bins")
        os.makedirs(bins_dir, exist_ok=True)
        cmd = ["metabat2", "-i", assembly, "-o", os.path.join(bins_dir, "bin"), "-t", str(threads)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running MetaBAT2: {e}")
            print(f"MetaBAT2 stderr: {e.stderr}")
            raise
        return bins_dir

    def run_checkm2(self, bins_dir, output_dir, threads):
        if not os.path.exists(self.checkm2_db_path):
            print("Downloading CheckM2 database...")
            subprocess.run(["checkm2", "database", "--download"], check=True)
        
        checkm2_out = os.path.join(output_dir, "checkm2_out")
        cmd = ["checkm2", "predict", "--threads", str(threads), "--input", bins_dir, "--output-directory", checkm2_out]
        subprocess.run(cmd, check=True)
        return checkm2_out

    
    def run_gtdbtk(self, bins_dir, output_dir, threads):
        if not os.path.exists(self.gtdbtk_db_path):
            print("Downloading GTDB-TK database...")
            subprocess.run(["gtdbtk", "download", "--all"], check=True)
        
        gtdbtk_out = os.path.join(output_dir, "gtdbtk_out")
        cmd = ["gtdbtk", "classify_wf", "--genome_dir", bins_dir, "--out_dir", gtdbtk_out, "--cpus", str(threads)]
        subprocess.run(cmd, check=True)
        return gtdbtk_out

    def run(self, input_file, output_dir, threads=1):
        assembly_dir = os.path.join(output_dir, "assembly")
        os.makedirs(assembly_dir, exist_ok=True)

        cmd = [
            self.flye_path,
            "--nano-raw", input_file,
            "--out-dir", assembly_dir,
            "--threads", str(threads)
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"Flye assembly completed successfully. Results saved in {assembly_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Flye: {e}")
            raise

        # Add additional functional analysis steps here if needed