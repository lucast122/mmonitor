#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=consider-using-f-string
"""
@author: kcurry
"""

import os
import argparse
import pathlib
import subprocess
from sys import stdout
from operator import add, mul
from pathlib import Path

import math
import pysam
import numpy as np
import pandas as pd
from flatten_dict import unflatten
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# static global variables
CIGAR_OPS = [1, 2, 4, 10]
CIGAR_OPS_ALL = [0, 1, 2, 4]
TAXONOMY_RANKS = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'clade', 'superkingdom']
RANKS_PRINTOUT = ['tax_id'] + TAXONOMY_RANKS + ['subspecies', 'species subgroup', 'species group']
RANKS_ORDER = ['tax_id'] + TAXONOMY_RANKS[:6] + TAXONOMY_RANKS[7:]


def validate_input(path):
    """Validate input file is either: fasta, fastq, or sam alignement file.

        path(str): path to input file
    """
    # pass if input is sam file
    sam_pysam = None
    try:
        # pylint: disable=maybe-no-member
        sam_pysam = pysam.AlignmentFile(path)
    except (ValueError, OSError):
        pass
    if sam_pysam:
        return

    # fail if input is not fasta/q
    fasta_rd, fastq_rd = None, None
    try:
        fasta_rd = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
        fastq_rd = SeqIO.to_dict(SeqIO.parse(path, "fastq"))
    except (UnicodeDecodeError, ValueError):
        pass
    if not (fasta_rd or fastq_rd):
        raise TypeError("Input must be in desired format: fasta or fastq")


def get_align_stats(alignment):
    """Retrieve list of inquired cigar stats (I,D,S,X) for alignment

        alignment (pysam.AlignmentFile): align of interest
        return (list(int)): list of counts for each cigar operation defined in (I,D,S,X)
    """
    cigar_stats = alignment.get_cigar_stats()[0]
    n_mismatch = cigar_stats[10] - cigar_stats[1] - cigar_stats[2]
    return [cigar_stats[1], cigar_stats[2], cigar_stats[4], n_mismatch]


def get_align_len(alignment):
    """Retrieve number of columns in alignment

        alignment (pysam.AlignmentFile): align of interest
        return (int): number of columns in alignment
    """
    return sum(alignment.get_cigar_stats()[0][cigar_op] for cigar_op in CIGAR_OPS_ALL)


def output_sequences(in_path, seq_output_path, input_type, keep_ids):
    """ Output specified list of sequences from input_file based on sequence id

        in_path (str): path to input fasta or fastq
        seq_output_path (str): output path for fasta/q of unclassified sequences
        input_type (str): fasta or fastq
        keep_ids (set): set of seuqence id strings
    """
    # Filter and write only the sequences with matching read IDs
    with open(in_path, "r", encoding="utf-8") as in_file, \
            open("{}.{}".format(seq_output_path, input_type), "w", encoding="utf-8") as out_seq_file:
        # Parse the FASTA file and filter by read IDs
        filtered_sequences = (seq for seq in SeqIO.parse(in_file, input_type) if seq.id in keep_ids)

        # Write the filtered sequences to the output file
        SeqIO.write(filtered_sequences, out_seq_file, input_type)


def get_cigar_op_log_probabilities(sam_path):
    """P(align_type) for each type in CIGAR_OPS by counting how often the corresponding
            operations occur in the primary alignments and by normalizing over the total
            sum of operations.

        sam_path(str): path to sam file of interest
        return: log probabilities (list(float)) for each cigar operation defined in CIGAR_OPS,
                where p > 0
            zero_locs (list(int)): list of indices (int) where probability == 0
            dict_longest_align (dict[str]:(int)): dict of max alignment length
                for each query read
    """
    cigar_stats_primary = [0] * len(CIGAR_OPS)
    dict_longest_align = {}
    # pylint: disable=maybe-no-member
    sam_pysam = pysam.AlignmentFile(sam_path)
    # add alignment lengths and adjust existing alignment lengths in dict if necessary
    for alignment in sam_pysam.fetch():
        align_len = get_align_len(alignment)
        if align_len not in dict_longest_align:
            dict_longest_align[alignment.query_name] = align_len
        if not alignment.is_secondary and not alignment.is_supplementary \
                and alignment.reference_name:
            cigar_stats_primary = list(map(add, cigar_stats_primary, get_align_stats(alignment)))
            # calculate cigar stats for alignment
            if dict_longest_align[alignment.query_name] < align_len:
                dict_longest_align[alignment.query_name] = align_len
    # check if any probabilities are 0, if so, remove
    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]
    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]
    n_char = sum(cigar_stats_primary)
    return [math.log(x) for x in np.array(cigar_stats_primary)/n_char], zero_locs, \
           dict_longest_align


def compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op, dict_longest_align, align_len):
    """ log(L(r|s)) = log(P(cigar_op)) × n_cigar_op for CIGAR_OPS

        alignment(pysam.AlignmentFile): pysam alignment to score
        cigar_stats(list(int)): list of cigar stats to compute
        log_p_cigar_op(list(float)): list of cigar_op probabilities corresponding to cigar_stats;
                                        computed from primary alignments
        dict_longest_align (dict[str]:(int)): dict of max alignment length for each query read
        align_len (int): number of columns in the alignment
        return: log_score (float): log(L(r|s))
                query_name (str): query name in alignment
                species_tid (int): species-level taxonomy id corresponding to ref_name
    """

    ref_name, query_name = alignment.reference_name, alignment.query_name
    log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats))) * \
                (dict_longest_align[query_name]/align_len)
    species_tid = int(ref_name.split(":")[0])
    return log_score, query_name, species_tid


def log_prob_rgs_dict(sam_path, log_p_cigar_op, dict_longest_align, p_cigar_op_zero_locs=None):
    """dict containing log(L(read|seq)) for all pairwise alignments in sam file

        sam_path(str): path to sam file
        log_p_cigar_op(list(float)): probability for each cigar operation defined in CIGAR_OPS,
                                         where p > 0
        dict_longest_align (dict[str]:(int)): dict of max alignment length for each query read
        zero_locs(list(int)): list of indices (int) where probability == 0
        return ({[str,int]:float}): dict[(query_name,ref_tax_id)]=log(L(query_name|ref_tax_id))
            int: unmapped read count
            int: mapped read count
    """
    # calculate log(L(read|seq)) for all alignments
    log_p_rgs, unmapped_set = {}, set()
    # pylint: disable=maybe-no-member
    sam_filename = pysam.AlignmentFile(sam_path, 'rb')

    if not p_cigar_op_zero_locs:
        for alignment in sam_filename.fetch():
            align_len = get_align_len(alignment)
            if alignment.reference_name and align_len:
                cigar_stats = get_align_stats(alignment)
                log_score, query_name, species_tid = \
                    compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op,
                                        dict_longest_align, align_len)

                if query_name not in log_p_rgs:
                    log_p_rgs[query_name] = ([species_tid], [log_score])
                elif query_name in log_p_rgs:
                    if species_tid not in log_p_rgs[query_name][0]:
                        log_p_rgs[query_name] = (log_p_rgs[query_name][0] + [species_tid],
                                                 log_p_rgs[query_name][1] + [log_score])
                    else:
                        logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
                        if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                            log_p_rgs[query_name][1][logprgs_idx] = log_score

            else:
                unmapped_set.add(alignment.query_name)
    else:
        for alignment in sam_filename.fetch():
            align_len = get_align_len(alignment)
            if alignment.reference_name and align_len:
                cigar_stats = get_align_stats(alignment)
                if sum(cigar_stats[x] for x in p_cigar_op_zero_locs) == 0:
                    for i in sorted(p_cigar_op_zero_locs, reverse=True):
                        del cigar_stats[i]
                    log_score, query_name, species_tid = \
                        compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op,
                                            dict_longest_align, align_len)

                    if query_name not in log_p_rgs:
                        log_p_rgs[query_name] = ([species_tid], [log_score])
                    elif query_name in log_p_rgs and species_tid not in log_p_rgs[query_name][0]:
                        log_p_rgs[query_name] = (log_p_rgs[query_name][0] +[species_tid],
                                                 log_p_rgs[query_name][1] + [log_score])
                    else:
                        logprgs_idx = log_p_rgs[query_name][0].index(species_tid)
                        if log_p_rgs[query_name][1][logprgs_idx] < log_score:
                            log_p_rgs[query_name][1][logprgs_idx] = log_score
            else:
                unmapped_set.add(alignment.query_name)

    mapped_set = set(log_p_rgs.keys())
    unmapped_set = unmapped_set - mapped_set
    unmapped_count = len(unmapped_set)
    stdout.write(f"Unmapped read count: {unmapped_count}\n")

    ## remove low likelihood alignments?
    ## remove if p(r|s) < 0.01
    #min_p_thresh = math.log(0.01)
    #log_p_rgs = {r_map: val for r_map, val in log_p_rgs.items() if val > min_p_thresh}
    return log_p_rgs, unmapped_set, mapped_set


def expectation_maximization(log_p_rgs, freq):
    """One iteration of the EM algorithm. Updates the relative abundance estimation in f based on
    probabilities in log_p_rgs.

    In each iteration, P(s|r), the probabilities of a read given a sequence, is calculated according
    to the current frequency vector. First, L(r|s) * f(s) is calculated, then a fixed multiplier, C,
    is calculated for each read. The fixed multiplier is then multiplied by L(r|s)*f(s),
    L(r|s) * f(s) * c. Then sum(L(r|s) * f(s) * c) for each sequence is calculated.
    Then (L(r|s) * f(s) * c) / (sum(L(r|s) * f(s) * c) for each sequence), is calculated.
    All calculations are done in log space. The frequency vector is then recalculated using the
    P(s|r) values, and the total log likelihood is updated.

    log_p_rgs({str:(int, float)}): dict[query_name]=(ref_tax_id, log(L(query_name|ref_tax_id)))
    freq{int:float}: dict[species_tax_id]:likelihood species is present in sample
    returns: f {int:float}: dict[species_tax_id]:updated likelihood species is present in sample
    total_log_likelihood (float): log likelihood updated f is accurate
    p_sgr {str: {int:float}}: probability of a read given the sequence
    """

    p_sgr_flat = {}
    logpr_sum, n_reads = 0, 0
    for read in log_p_rgs:
        valid_seqs, log_p_rns = [], []
        # check if sequences were found in frequency vector
        for seq in range(len(log_p_rgs[read][0])):
            s_val = log_p_rgs[read][0][seq]
            if s_val in freq and freq[s_val] != 0:
                logprns_val = log_p_rgs[read][1][seq] + math.log(freq[s_val])
                # calculates log(L(r|s))+log(f(s)) for each sequence
                valid_seqs.append(s_val)
                log_p_rns.append(logprns_val)

        if len(valid_seqs) != 0:
            #np.array([log_p_rns], dtype=np.float64)
            logc = -np.max(log_p_rns)  # calculate fixed multiplier, c
            prnsc = np.exp(log_p_rns + logc)  # calculates exp(log(L(r|s) * f(s) * c))
            prc = np.sum(prnsc)  # calculates sum of (L(r|s) * f(s) * c) for each read
            logpr_sum += (np.log(prc) - logc)  # add to sum of log likelihood
            n_reads += 1
            # calculates P(s|r) for each sequence
            for seq in enumerate(valid_seqs):
                p_sgr_flat[(seq[1], read)] = prnsc[seq[0]] / prc
    p_sgr = unflatten(p_sgr_flat)
    # calculates updated frequency vector
    frq = {tax_id: sum(read_id.values()) / n_reads for tax_id, read_id in p_sgr.items()}
    return frq, logpr_sum, p_sgr


def expectation_maximization_iterations(log_p_rgs, db_ids, lli_thresh, input_threshold):
    """Full expectation maximization algorithm for alignments in log_L_rgs dict.
    Calls the expectation_maximization function during each iteration of the algorithm.
    Stops iterations once the log likelihood is calculated to have increased less than threshold.

    log_p_rgs{[str,int]:float}: dict[(query_name,ref_tax_id)]=log(L(query_name|ref_tax_id))
    db_ids(list(int)): list of each unique species taxonomy id present in database
    lli_thresh(float): log likelihood increase minimum to continue EM iterations
    input_threshold(float): minimum relative abundance in output
    return: {int:float}: dict[species_tax_id]:estimated likelihood species is present in sample
            float: min abundance threshold
    """
    n_db = len(db_ids)
    n_reads = len(log_p_rgs)
    stdout.write("Mapped read count: {}\n".format(n_reads))
    # check if there are enough reads
    if n_reads == 0:
        raise ValueError("0 reads mapped")
    freq, counter = dict.fromkeys(db_ids, 1 / n_db), 1

    # set output abundance threshold
    freq_thresh = 1/(n_reads + 1)
    if n_reads > 1000:
        freq_thresh = 10/n_reads

    # performs iterations of the expectation_maximization algorithm
    total_log_likelihood = -math.inf
    while True:
        freq, updated_log_likelihood, _ = expectation_maximization(log_p_rgs, freq)

        # check f vector sums to 1
        freq_sum = sum(freq.values())
        if not .9 <= freq_sum <= 1.1:
            raise ValueError("f sums to {}, rather than 1".format(freq_sum))

        # confirm log likelihood increase
        log_likelihood_diff = updated_log_likelihood - total_log_likelihood
        total_log_likelihood = updated_log_likelihood
        if log_likelihood_diff < 0:
            raise ValueError("total_log_likelihood decreased from prior iteration")

        # exit loop if log likelihood increase less than threshold
        if log_likelihood_diff < lli_thresh:
            stdout.write("Number of EM iterations: {}\n".format(counter))
            freq = {k: v for k, v in freq.items() if v >= freq_thresh}
                # remove tax id if less than the frequency threshold
            freq_full, updated_log_likelihood, p_sgr = expectation_maximization(log_p_rgs, freq)
            freq_set_thresh = None
            if freq_thresh < input_threshold:
                freq = {k: v for k, v in freq_full.items() if v >= input_threshold}
                freq_set_thresh, updated_log_likelihood, p_sgr = \
                    expectation_maximization(log_p_rgs, freq)
            return freq_full, freq_set_thresh, p_sgr

        #output current estimation
        #freq_to_lineage_df(freq, f"{out_file}_{counter}", df_nodes, df_names)
        counter += 1


def lineage_dict_from_tid(taxid, nodes_dict, names_dict):
    """For each given taxid, traverse the node_dict to build the lineage, and record the tax id
            for each taxonomic rank using the names_dict.

        tid(int): tax id to retrieve lineage dict
        nodes_dict{int:(int, str)}: dict of nodes.dmp with 'tax_id' as keys,
                                        tuple ('parent_taxid', 'rank') as values
        names_dict{int:str}: dict of names.dmp with 'tax_id' as keys,
                                        'name_txt' as values
        return (tuple): a tuple containing the scientific name for each taxonomic
            rank for the given taxid
    """
    # initialize list and record tax id
    lineage_list = [taxid] + [""] * (len(RANKS_PRINTOUT) - 1)
    # traverse the nodes to create the lineage
    while names_dict[taxid] != "root":
        tup = nodes_dict[taxid]
        # find the name for each taxonomic rank
        if tup[1] in RANKS_PRINTOUT: # check rank in printout list
            idx = RANKS_PRINTOUT.index(tup[1])
            lineage_list[idx] = names_dict[taxid]
        taxid = tup[0]
    return tuple(lineage_list)


def freq_to_lineage_df(freq, tsv_output_path, taxonomy_df, mapped_count,
                       unmapped_count, mapped_unclassified_count, counts=False):
    """Converts freq to a pandas df where each row contains abundance and tax lineage for
                classified species in f.keys(). Stores df as .tsv file in tsv_output_path.

        freq{int:float}: dict[species_tax_id]:estimated likelihood species is present in sample
        tsv_output_path(str): path to output .tsv file
        taxonomy_df(df): pandas df of all db sequence taxonomy with index 'tax_id'
        mapped_count(int): number of mapped reads
        unmapped_count(int): number of unmapped reads
        mapped_unclassified_count(int): number of that mapped but were assigned due to
                                        low abundant classification
        counts(boolean): True if include estimated counts in output .tsv file
        returns(df): pandas df with lineage and abundances for values in f
    """
    #add tax lineage for values in freq
    results_df = pd.DataFrame(zip(list(freq.keys()) + ['unmapped', 'mapped_unclassified'],
                                  list(freq.values()) + [0, 0]),
                              columns=["tax_id", "abundance"]).set_index('tax_id')
    results_df = results_df.join(taxonomy_df, how='left').reset_index()
        #add in the estimated count values for the mapped and unmapped counts
    if counts:
        classified_count = mapped_count - mapped_unclassified_count
        counts_series = pd.concat([(results_df["abundance"] * classified_count)[:-2],
                                   pd.Series(unmapped_count), pd.Series(mapped_unclassified_count)],
                                    ignore_index=True)
        results_df["estimated counts"] = counts_series
    results_df.to_csv("{}.tsv".format(tsv_output_path), sep='\t', index=False)
    return results_df


def generate_alignments(in_file_list, out_basename, database):
    """ Generate .sam alignment file

        in_file_list(list(str)): list of path(s) to input sequences
        out_basename(str): path and basename for output files
    """
    #indicate filepath and create file
    input_file = " ".join(in_file_list)
    filetype = pathlib.PurePath(args.input_file[0]).suffix
    if filetype == '.sam':
        args.keep_files = True
        return input_file
    sam_align_file = "{}_emu_alignments.sam".format(out_basename)
    db_sequence_file = os.path.join(database, 'species_taxid.fasta')

    # force minimap2 to consider the forward transcript strand only
    if args.mm2_forward_only:
        subprocess.check_output("minimap2 -ax {} -t {} -N {} -p .9 -u f -K {} {} {} -o {}".
                                format(args.type, args.threads, args.N, args.K,
                                    db_sequence_file, input_file, sam_align_file),
                                shell=True)
    else:
        subprocess.check_output("minimap2 -ax {} -t {} -N {} -p .9 -K {} {} {} -o {}".
                                format(args.type, args.threads, args.N, args.K,
                                    db_sequence_file, input_file, sam_align_file),
                                shell=True)

    return sam_align_file


def output_read_assignments(p_sgr, tsv_output_path):
    """ Output file of read assignment distributions for all

        p_sgr({tid:{read_id:probability}}): P(s|r), likelihood read r emanates db seq s
        tsv_output_path(str): path to output .tsv file
        returns(df): pandas df of read assignment distributions
    """
    dist_df = pd.DataFrame(p_sgr)
    dist_df.to_csv("{}.tsv".format(tsv_output_path), sep='\t')
    return dist_df


def create_nodes_dict(nodes_path):
    """convert nodes.dmp file into a dictionary

        nodes_path(str): path to nodes.dmp file
        returns: nodes_dict{int:[int, str]}: dict of nodes.dmp with 'tax_id' as keys,
                tuple ('parent_taxid', 'rank') as values

    """
    node_headers = ['tax_id', 'parent_tax_id', 'rank']
    nodes_df = pd.read_csv(nodes_path, sep='|', header=None, dtype=str)[[0, 1, 2]]
    nodes_df.columns = node_headers
    for col in nodes_df.columns:
        nodes_df[col] = nodes_df[col].str.strip()
    return dict(zip(nodes_df['tax_id'], tuple(zip(nodes_df['parent_tax_id'], nodes_df['rank']))))


def create_names_dict(names_path):
    """convert names.dmp file into a dictionary

        names_path(str): path to names.dmp file
        returns: names_dict{int:str}: pandas dict of names.dmp with 'tax_id' as keys,
                'name_txt' as values
    """
    name_headers = ['tax_id', 'name_txt', 'name_class']
    names_df = pd.read_csv(names_path, sep='|', index_col=False, header=None, dtype=str)\
        .drop([2, 4], axis=1)
    names_df.columns = name_headers
    for col in names_df.columns:
        names_df[col] = names_df[col].str.strip()
    names_df = names_df[names_df["name_class"] == "scientific name"]
    return dict(zip(names_df['tax_id'], names_df['name_txt']))


def get_species_tid(tid, nodes_dict):
    """ Get lowest taxid down to species-level in lineage for taxid [tid]

        tid(int): taxid for species level or more specific
        nodes_dict{int:[int, str]}: dict of nodes.dmp with 'tax_id' as keys,
                                        tuple ('parent_taxid', 'rank') as values
        return(int): species taxid in lineage
    """
    if str(tid) not in nodes_dict.keys():
        raise ValueError("Taxid:{} not found in nodes file".format(tid))
    while nodes_dict[str(tid)][1] not in TAXONOMY_RANKS:
        tid = nodes_dict[str(tid)][0]
    return tid


def create_species_seq2tax_dict(seq2tax_path, nodes_dict):
    """Convert seqid-taxid mapping in seq2tax_path to dict mapping seqid to species level taxid

        seq2tax_path(str): path to seqid-taxid mapping file
        nodes_dict{int:[int, str]}: dict of nodes.dmp with 'tax_id' as keys,
                                        tuple ('parent_taxid', 'rank') as values
        returns {str:int}: dict[seqid] = species taxid
    """
    seq2tax_dict, species_id_dict = {}, {}
    with open(seq2tax_path, encoding="utf8") as file:
        # unpack values in each line of the file
        for line in file:
            (seqid, tid) = line.rstrip().split("\t")
            # retrieve the species taxid from the species_id_dict if already in dictionary
            if tid in species_id_dict:
                species_tid = species_id_dict[tid]
            # find the species taxid if not in the dictionary
            else:
                species_tid = get_species_tid(tid, nodes_dict)
                species_id_dict[tid] = species_tid
            seq2tax_dict[seqid] = species_tid
    return seq2tax_dict


def create_direct_seq2tax_dict(seq2tax_path):
    """Convert seqid-taxid mapping in seq2tax_path to dict mapping
        SILVA seqid to corresponding taxid

        seq2tax_path(str): path to seqid-taxid mapping file
        returns {str:int}: dict[seqid] = species taxid
    """
    with open(seq2tax_path, encoding="utf8") as file:
        seq2_taxid = {}
        for line in file:
            (seqid, taxid) = line.rstrip().split("\t")
            seq2_taxid[seqid] = taxid
    return seq2_taxid


def create_unique_seq_dict(db_fasta_path, seq2tax_dict):
    """ Creates dict of unique sequences to for the sequences connected to each species taxid

        db_fasta_path(str): path to fasta file of database sequences
        seq2tax_dict{str:int}: dict[seqid] = species taxid
        returns {str:{int:[str]}}: dict[seq] = {species_taxid: [list of sequence ids]}
    """
    fasta_dict = {}
    # traverse through the species taxids
    for record in SeqIO.parse(db_fasta_path, "fasta"):
        tid = seq2tax_dict[record.id]
        if tid:
            # if sequence already in the dictionary, add more ids if found
            if record.seq in fasta_dict:
                if tid in fasta_dict[record.seq].keys():
                    fasta_dict[record.seq][tid] += [record.description]
                # create inner species taxid dictionary and add id
                else:
                    fasta_dict[record.seq][tid] = [record.description]
            elif record.seq.reverse_complement() in fasta_dict:
                if tid in fasta_dict[record.seq.reverse_complement()].keys():
                    fasta_dict[record.seq.reverse_complement()][tid] += [record.description]
                else:
                    fasta_dict[record.seq.reverse_complement()][tid] = [record.description]
            else:
                fasta_dict[record.seq] = {tid: [record.description]}
    return fasta_dict


def create_reduced_fasta(fasta_dict, db_name):
    """ Creates fasta file of taxid for each sequences in fasta_dict with id
            'species_taxid:db_name:sequence_id'

        fasta_dict{str:{int:[str]}}: dict[seq] = {species_taxid: [list of sequence ids]}
        db_name(str): name to represent database represented in fasta_dict
        returns (list[Bio.SeqRecord]): list of sequences for output fasta file
    """
    records, count = [], 1
    for seq, tid_dict in fasta_dict.items():
        for taxid, descriptions in tid_dict.items():
            records += [SeqRecord(seq, id="{}:{}:{}".format(taxid, db_name, count),
                                  description="{}".format(descriptions))]
            count += 1
    return records


def build_ncbi_taxonomy(unique_tids, nodes_dict, names_dict, filepath):
    """Creates a tsv file where for each id in unique_tids a tax lineage is written to the file.
            each value of the tax lineage is seperated by a tab.

        unique_tids(set): set of ints of each unique taxid in list of database sequences
        nodes_dict{int:[int, str]}: dict of nodes.dmp with 'tax_id' as keys,
                                        tuple ('parent_taxid', 'rank') as values
        names_dict{int:str}: pandas dict of names.dmp with 'tax_id' as keys,
                                        'name_txt' as values
    """
    with open(filepath, "w", encoding="utf8") as file:
        # write the header to the file
        dummy_str = '\t'.join(['%s',] * len(RANKS_PRINTOUT)) + '\n'
        file.write(dummy_str % tuple(RANKS_PRINTOUT))
        # write each lineage as a row in the file
        for tid in unique_tids:
            lst = lineage_dict_from_tid(tid, nodes_dict, names_dict)
            file.write(dummy_str % lst)


def build_direct_taxonomy(tid_set, lineage_path, taxonomy_file):
    """Create a tsv file that contains the taxid and the corresponding
        lineage for the SILVA database

        tid_set(set): set containing all taxids to add to the database
        lineage_path(str): path to file containing the lineage in the first
                column and the taxid in next column
        taxonomy_file(str): path to output tsv file containing lineages
    """
    with open(taxonomy_file, 'w', encoding="utf8") as tax_output_file:
        with open(lineage_path, encoding="utf8") as file:
            first_line = file.readline()
            tax_output_file.write("{}\t".format(RANKS_PRINTOUT[0]))
            tax_output_file.write(first_line.split("\t", 1)[1])
            for line in file:
                tax_id = line.split("\t", 1)[0]
                if tax_id in tid_set: #only add if taxid is in fasta
                    tax_output_file.write(line)


def collapse_rank(path, rank):
    """ Stores a version of emu-output (path) collapsed at the specified taxonomic rank in same
            folder as input.

        path(str): path to emu output
        rank(str): taxonomic rank for collapsed abundance: ["species", "genus", "family",
            "order", "class", "phylum", "clade", "superkingdom"]
    """
    df_emu = pd.read_csv(path, sep='\t')
    if rank not in TAXONOMY_RANKS:
        raise ValueError("Specified rank must be in list: {}".format(TAXONOMY_RANKS))
    keep_ranks = TAXONOMY_RANKS[TAXONOMY_RANKS.index(rank):]
    for keep_rank in keep_ranks:
        if keep_rank not in df_emu.columns:
            keep_ranks.remove(keep_rank)
    if "estimated counts" in df_emu.columns:
        df_emu_copy = df_emu[['abundance', 'estimated counts'] + keep_ranks]
        df_emu_copy = df_emu_copy.replace({'-': 0})
        df_emu_copy = df_emu_copy.astype({'abundance': 'float', 'estimated counts': 'float'})
    else:
        df_emu_copy = df_emu[['abundance'] + keep_ranks]
        df_emu_copy = df_emu_copy.replace({'-': 0})
        df_emu_copy = df_emu_copy.astype({'abundance': 'float'})
    df_emu_copy = df_emu_copy.groupby(keep_ranks, dropna=False).sum()
    output_path = "{}-{}.tsv".format(os.path.splitext(path)[0], rank)
    df_emu_copy.to_csv(output_path, sep='\t')
    stdout.write("File generated: {}\n".format(output_path))

def combine_outputs(dir_path, rank, split_files=False, count_table=False):
    """ Combines multiple Emu output relative abundance tables into a single table.
        Inlcudes all .tsv files with 'rel-abundance' in the filename in the provided dir_path.

        dir_path(str): path of directory containing Emu output files to combine
        rank(str): taxonomic rank to combine files on
        return(df): Pandas df of the combined relative abundance files
    """
    keep_ranks = RANKS_ORDER[RANKS_ORDER.index(rank):]
    df_combined_full = pd.DataFrame(columns=keep_ranks, dtype=str)
    metric = 'abundance'
    if count_table:
        metric = 'estimated counts'
    for file in os.listdir(dir_path):
        file_extension = pathlib.Path(file).suffix
        if file_extension == '.tsv' and 'rel-abundance' in file:
            name = pathlib.Path(file).stem
            name = name.replace('_rel-abundance', '')
            df_sample = pd.read_csv(os.path.join(dir_path, file), sep='\t', dtype=str)
            df_sample[[metric]] = df_sample[[metric]].apply(pd.to_numeric)
            if rank in df_sample.columns and metric in df_sample.columns:
                #check which keep_ranks are in df_sample
                keep_ranks_sample = [value for value in keep_ranks
                                     if value in set(df_sample.columns)]
                if df_sample.at[len(df_sample.index)-1, 'tax_id'] == 'unmapped':
                    df_sample.at[len(df_sample.index)-1, rank] = 'unmapped'
                df_sample_reduced = df_sample[keep_ranks_sample +
                                              [metric]].rename(columns={metric: name})
                df_sample_reduced = df_sample_reduced.groupby(keep_ranks_sample, dropna=False)\
                    .sum().reset_index() #sum metric within df_sample_reduced if same tax lineage
                df_sample_reduced = df_sample_reduced.astype(object)
                df_sample_reduced[[name]] = df_sample_reduced[[name]].apply(pd.to_numeric)
                df_combined_full = pd.merge(df_combined_full, df_sample_reduced, how='outer')
    df_combined_full = df_combined_full.set_index(rank).sort_index().reset_index()
    filename_suffix = ""
    if count_table:
        filename_suffix = "-counts"
    if split_files:
        abundance_out_path = os.path.join(dir_path, "emu-combined-abundance-{}{}.tsv"
                                          .format(rank, filename_suffix))
        tax_out_path = os.path.join(dir_path, "emu-combined-taxonomy-{}.tsv".format(rank))
        stdout.write("Combined taxonomy table generated: {}\n".format(tax_out_path))
        df_combined_full[keep_ranks].to_csv(tax_out_path, sep='\t', index=False)
        keep_ranks.remove(rank)
        df_combined_full.drop(columns=keep_ranks).to_csv(abundance_out_path, sep='\t', index=False)
        stdout.write("Combined abundance table generated: {}\n".format(abundance_out_path))
    else:
        out_path = os.path.join(dir_path, "emu-combined-{}{}.tsv".format(rank, filename_suffix))
        df_combined_full.to_csv(out_path, sep='\t', index=False)
        stdout.write("Combined table generated: {}\n".format(out_path))
    return df_combined_full

if __name__ == "__main__":
    __version__ = "3.5.0"
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', '-v', action='version', version='%(prog)s v' + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help='sub-commands')
    abundance_parser = subparsers.\
        add_parser("abundance", help="Generate relative abundance estimates")
    abundance_parser.add_argument(
        'input_file', type=str, nargs='+',
        help='filepath to input nt sequence file')
    abundance_parser.add_argument(
        '--type', '-x', choices=['map-ont', 'map-pb', 'sr'], default='map-ont',
        help='short-read: sr, Pac-Bio:map-pb, ONT:map-ont [map-ont]')
    abundance_parser.add_argument(
        '--min-abundance', '-a', type=float, default=0.0001,
        help='min species abundance in results [0.0001]')
    abundance_parser.add_argument(
        '--db', type=str, default=os.environ.get("EMU_DATABASE_DIR"),
        help='path to emu database containing: names_df.tsv, '
             'nodes_df.tsv, species_taxid.fasta, unqiue_taxids.tsv [$EMU_DATABASE_DIR]')
    abundance_parser.add_argument(
        '--N', '-N', type=int, default=50,
        help='minimap max number of secondary alignments per read [50]')
    abundance_parser.add_argument(
        '--K', '-K', type=int, default=500000000,
        help='minibatch size for minimap2 mapping [500M]')
    abundance_parser.add_argument(
        '--mm2-forward-only', action="store_true",
        help='force minimap2 to consider the forward transcript strand only')
    abundance_parser.add_argument(
        '--output-dir', type=str, default="./results",
        help='output directory name [./results]')
    abundance_parser.add_argument(
        '--output-basename', type=str,
        help='basename for all emu output files [{input_file}]')
    abundance_parser.add_argument(
        '--keep-files', action="store_true",
        help='keep working files in output-dir')
    abundance_parser.add_argument(
        '--keep-counts', action="store_true",
        help='include estimated read counts in output')
    abundance_parser.add_argument(
        '--keep-read-assignments', action="store_true",
        help='output file of read assignment distribution')
    abundance_parser.add_argument(
        '--output-unclassified', action="store_true",
        help='output unclassified sequences')
    abundance_parser.add_argument(
        '--threads', type=int, default=3,
        help='threads utilized by minimap [3]')

    build_db_parser = subparsers.add_parser("build-database",
                                            help="Build custom Emu database")
    build_db_parser.add_argument(
        'db_name', type=str,
        help='custom database name')
    build_db_parser.add_argument(
        '--sequences', type=str, required=True,
        help='path to fasta of database sequences')
    build_db_parser.add_argument(
        '--seq2tax', type=str, required=True,
        help='path to tsv mapping species tax id to fasta sequence headers')
    taxonomy_group = build_db_parser.add_mutually_exclusive_group(required=True)
    taxonomy_group.add_argument(
        '--ncbi-taxonomy', type=str,
        help='path to directory containing both a names.dmp and nodes.dmp file')
    taxonomy_group.add_argument(
        '--taxonomy-list', type=str,
        help='path to .tsv file mapping full lineage to corresponding taxid')

    collapse_parser = subparsers.add_parser("collapse-taxonomy",
                                            help="Collapse emu output at specified taxonomic rank")
    collapse_parser.add_argument(
        'input_path', type=str,
        help='emu output filepath')
    collapse_parser.add_argument(
        'rank', type=str,
        help='collapsed taxonomic rank')

    combine_parser = subparsers.add_parser("combine-outputs",
                help="Combine Emu rel abundance outputs to a single table")
    combine_parser.add_argument(
        'dir_path', type=str,
        help='path to directory containing Emu output files')
    combine_parser.add_argument(
        'rank', type=str,
        help='taxonomic rank to include in combined table')
    combine_parser.add_argument(
        '--split-tables', action="store_true",
        help='two output tables:abundances and taxonomy lineages')
    combine_parser.add_argument(
        '--counts', action="store_true",
        help='counts rather than abundances in output table')
    args = parser.parse_args()

    if args.subparser_name == "abundance":
        # check input file is fasta/q or sam alignment file
        # validate_input(args.input_file[0])
        # convert taxonomy files to dataframes
        if not args.db:
            raise ValueError("Database not specified. "
                             "Either 'export EMU_DATABASE_DIR=<path_to_database>' or "
                             "utilize '--db' parameter.")
        df_taxonomy = pd.read_csv(os.path.join(args.db, "taxonomy.tsv"), sep='\t',
                                  index_col='tax_id', dtype=str)
        db_species_tids = df_taxonomy.index

        # set up output paths
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        out_file = os.path.join(args.output_dir, "-".join([Path(v).stem for v in args.input_file]))
        if args.output_basename:
            out_file = os.path.join(args.output_dir, args.output_basename)

        # perform EM algorithm & generate output
        SAM_FILE = generate_alignments(args.input_file, out_file, args.db)
        log_prob_cigar_op, locs_p_cigar_zero, longest_align_dict = \
            get_cigar_op_log_probabilities(SAM_FILE)
        log_prob_rgs, set_unmapped, set_mapped = log_prob_rgs_dict(
            SAM_FILE, log_prob_cigar_op, longest_align_dict, locs_p_cigar_zero)
        f_full, f_set_thresh, read_dist = expectation_maximization_iterations(log_prob_rgs,
                                                                   db_species_tids,
                                                                   .01, args.min_abundance)
        classified_reads = {read_id for inner_dict in read_dist.values() for read_id in inner_dict}
        mapped_unclassified = set_mapped - classified_reads
        stdout.write(f"Unclassified mapped read count: {len(mapped_unclassified)}\n")
        freq_to_lineage_df(f_full, "{}_rel-abundance".format(out_file), df_taxonomy,
                           len(set_mapped), len(set_unmapped), len(mapped_unclassified),
                           args.keep_counts)


        # output read assignment distributions as a tsv
        if args.keep_read_assignments:
            output_read_assignments(read_dist, "{}_read-assignment-distributions".format(out_file))

        # convert and save frequency to a tsv
        if f_set_thresh:
            freq_to_lineage_df(
                f_set_thresh,
                "{}_rel-abundance-threshold-{}".format(out_file, args.min_abundance),
                df_taxonomy, len(set_mapped), len(set_unmapped), len(mapped_unclassified),
                args.keep_counts)

        # gather input sequences that are unmapped according to minimap2
        if args.output_unclassified:
            ext = os.path.splitext(args.input_file[0])[-1]
            INPUT_FILETYPE = "fasta"
            if ext in [".fastq", ".fq"]:
                INPUT_FILETYPE = "fastq"
            output_sequences(args.input_file[0], "{}_unmapped".format(out_file), INPUT_FILETYPE,
                                                set_unmapped)
            output_sequences(args.input_file[0], "{}_unclassified_mapped".format(out_file),
                             INPUT_FILETYPE, mapped_unclassified)

        # clean up extra file
        if not args.keep_files:
            if os.path.exists(SAM_FILE):
                os.remove(SAM_FILE)


    # create a custom database with ncbi taxonomy
    if args.subparser_name == "build-database":
        emu_db_path = os.getcwd()
        custom_db_path = os.path.join(emu_db_path, args.db_name)
        if not os.path.exists(custom_db_path):
            os.makedirs(custom_db_path)
        stdout.write("Emu custom database generating at path: {} ...\n".format(custom_db_path))

        # set up seq2tax dict for either NCBI or direct taxonomy
        if args.ncbi_taxonomy:
            dict_names = create_names_dict(os.path.join(args.ncbi_taxonomy, 'names.dmp'))
            dict_nodes = create_nodes_dict(os.path.join(args.ncbi_taxonomy, 'nodes.dmp'))
            seq2tax = create_species_seq2tax_dict(args.seq2tax, dict_nodes)
        else:
            seq2tax = create_direct_seq2tax_dict(args.seq2tax)

        # print fasta in desired Emu database format
        db_unique_ids = set(seq2tax.values())
        dict_fasta = create_unique_seq_dict(args.sequences, seq2tax)
        fasta_records = create_reduced_fasta(dict_fasta, args.db_name)
        SeqIO.write(fasta_records, os.path.join(custom_db_path, 'species_taxid.fasta'), "fasta")

        # build taxonomy for database
        output_taxonomy_location = os.path.join(custom_db_path, "taxonomy.tsv")
        if args.ncbi_taxonomy:
            build_ncbi_taxonomy(db_unique_ids, dict_nodes, dict_names, output_taxonomy_location)
        else:
            build_direct_taxonomy(db_unique_ids, args.taxonomy_list, output_taxonomy_location)
        stdout.write("Database creation successful\n")


    # collapse Emu results at desired taxonomic rank
    if args.subparser_name == "collapse-taxonomy":
        collapse_rank(args.input_path, args.rank)

    # combine Emu results at desired taxonomic rank
    if args.subparser_name == "combine-outputs":
        combine_outputs(args.dir_path, args.rank, args.split_tables, args.counts)