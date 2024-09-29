import gzip
import logging
import multiprocessing
import os
import subprocess

import pandas as pd

from build_mmonitor_pyinstaller import ROOT
from lib import emu


class EmuRunner:

    def __init__(self):
        self.concat_file_name = None
        self.logger = logging.getLogger('timestamp')
        self.check_emu()
        self.emu_out = ""

    @staticmethod
    def unpack_fastq_list(ls):
        """
        Gets a list of paths and outputs a comma seperated string containing the paths used as input for emu.py
        """
        if len(ls) == 1:
            return f"{ls[0]}"
        elif len(ls) > 1:
            return ",".join(ls)

    def check_emu(self):
        try:
            subprocess.run([f"{ROOT}/lib/emu.py", '-h'], stdout=open(os.devnull, 'w'),
                            stderr=subprocess.STDOUT)
        except FileNotFoundError:
            self.logger.error(
                "Make sure that emu.py is installed and on the sytem path. For more info visit http://www.ccb.jhu.edu/software/centrifuge/manual.shtml")

    def run_emu(self, sequence_list, sample_name, min_abundance):
        print(f"Running emu with min abundance of {min_abundance}")
        self.emu_out = f"{ROOT}/src/resources/pipeline_out/{sample_name}/"

        #remove concatenated files from sequence list to avoid concatenating twice
        sequence_list = [s for s in sequence_list if "concatenated" not in s]
        # print(sequence_list)
        concat_file_name = f"{os.path.dirname(sequence_list[0])}/{sample_name}_concatenated.fastq.gz"
        self.concat_file_name = concat_file_name
        # print(f"concat_file_name: {concat_file_name}")
        if not os.path.exists(concat_file_name):
            self.concatenate_fastq_files(sequence_list, concat_file_name)


        if ".fasta" in sequence_list[0] or ".fa" in sequence_list[0] or ".fastq" in sequence_list[0]\
                or ".fq" in sequence_list[0]:
            # alternative way to call emu. this doesn't work with frozen build from pyinstaller so comment out for now
            # cmd = f"python {ROOT}/lib/emu.py abundance {concat_file_name} --db {ROOT}/src/resources/emu_db/" \
            #       f" --output-dir {self.emu_out} --threads {multiprocessing.cpu_count()} --type map-ont --output-basename {sample_name}"
            # print(cmd)
            # os.system(cmd)
            emu_db = f"{ROOT}/src/resources/emu_db/"

            df_taxonomy = pd.read_csv(os.path.join(emu_db, "taxonomy.tsv"), sep='\t',
                                      index_col='tax_id', dtype=str)
            db_species_tids = df_taxonomy.index
            print(f"Emu out: {self.emu_out}")
            if not os.path.exists(self.emu_out):
                os.makedirs(self.emu_out)

            out_file_base = self.emu_out
            sam_out = f"{out_file_base}/emu_alignments.sam"
            tsv_out = f"{out_file_base}/{sample_name}_rel-abundance"
            # print(f"Out file: {out_file_base}")
            print(f"min abundance: {min_abundance}")
            SAM_FILE = emu.generate_alignments(concat_file_name, sam_out, emu_db, "map-ont",
                                               f"{multiprocessing.cpu_count()}", 50, 500000000)
            log_prob_cigar_op, locs_p_cigar_zero, longest_align_dict = \
                emu.get_cigar_op_log_probabilities(SAM_FILE)
            log_prob_rgs, counts_unassigned, counts_assigned = emu.log_prob_rgs_dict(
                SAM_FILE, log_prob_cigar_op, longest_align_dict, locs_p_cigar_zero)
            f_full, f_set_thresh, read_dist = emu.expectation_maximization_iterations(log_prob_rgs,
                                                                                      db_species_tids,
                                                                                      .01,
                                                                                      input_threshold=min_abundance)
            # print(f_full)
            # print(f_set_thresh)

            emu.freq_to_lineage_df(f_full, tsv_out, df_taxonomy,
                                   counts_assigned, counts_unassigned, True)

            # convert and save frequency to a tsv
            if f_set_thresh:
                emu.freq_to_lineage_df(
                    f_set_thresh,
                    os.path.join(out_file_base, f"{sample_name}_rel-abundance-threshold"),

                    df_taxonomy, counts_assigned, counts_unassigned, True)
        # remove concatenated file after processing
        # TODO: calculate statistics and then remove, for now don't remove
        # os.remove(concat_file_name)

    def get_files_from_folder(self, folder_path):
        """
        Gets a path to a folder, checks if path contains sequencing files with specified endings and returns list
        containing paths to sequencing files.
        """
        files = []
        found = False
        try:
            for file in os.listdir(folder_path):
                # print(file)
                if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fasta") or file.endswith(
                        ".fastq.gz"):
                    files.append(f"{folder_path}/{file}")
                    found = True
            if not found:
                self.logger.error(f"No sequencing files (.fastq, .fq found at {folder_path}")
            return files
        except FileNotFoundError:
            self.logger.error(f"Invalid folder path")

    def concatenate_fastq_files(self, input_files, output_file):
        with gzip.open(output_file, 'wt') as output:
            for input_file in input_files:
                is_gzipped = input_file.endswith(".gz")
                open_func = gzip.open if is_gzipped else open
                mode = 'rt' if is_gzipped else 'r'

                with open_func(input_file, mode) as input2:
                    for line in input2:
                        output.write(line)
