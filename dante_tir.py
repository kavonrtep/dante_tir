#!/usr/bin/env python3
"""
Script for identification DNA transposons with Terminal Inverted Repeats (TIRs) based
DANTE annotation of conserved domains.
"""

import argparse
import glob
import os
import shutil
import subprocess
from itertools import chain
import random
import hashlib
from version import __version__
from time import time

import dt_utils as dt
from multiprocessing import Pool

def compute_file_checksum(filepath):
    """Compute MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Program  for identification DNA '
                    'transposons with Terminal Inverted '
                    'Repeats (TIRs) basen DANTE annotation of conserved domains of '
                    'transposases'
        )
    parser.add_argument(
        '-g', '--gff3',
        help='GFF3 file with DANTE annotation of conserved domains of transposases',
        required=True
        )
    parser.add_argument(
        '-f', '--fasta', help='FASTA file with genome assembly', required=True
        )
    parser.add_argument(
        '-o', '--output_dir', help='Output directory with TIRs', required=True
        )
    parser.add_argument('-c', '--cpu', help='Number of CPUs to use',
                        type=int, default=1)
    parser.add_argument(
        '--max_class_size', help='Maximum number of sequences per class before splitting for CAP3 assembly',
        type=int, default=None
    )
    parser.add_argument(
        '--version', action='version',
        version='%(prog)s {version}'.format(version=__version__)
        )
    # add debug option
    parser.add_argument(
        '--debug', action='store_true',
        help='Debug mode', default=False
    )
    parser.add_argument(
        '--seed', help='Random seed for reproducibility (used in R scripts)',
        type=int, default=42
    )
    parser.add_argument(
        '--n_beast_iter', help='Number of BEAST iterations with different random seeds for TIR detection (default=1, max=10)',
        type=int, default=1
    )

    print("--------------------------------------------------------")
    print("")
    print("          ---|>>>---- DANTE_TIR---- <<<|---")
    print("")
    print("Domain Based identification of DNA transposons with TIRs")
    print("")
    print("--------------------------------------------------------")




    args = parser.parse_args()

    # Validate n_beast_iter parameter
    if args.n_beast_iter < 1 or args.n_beast_iter > 10:
        parser.error("--n_beast_iter must be between 1 and 10 (default=1)")

    script_dir = os.path.dirname(os.path.realpath(__file__))
    # set random seed for reproducibility
    random.seed(42)

    # create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # create working directory in output directory
    args.working_dir = F'{args.output_dir}/working_dir'
    if not os.path.exists(args.working_dir):
        os.mkdir(args.working_dir)



    # Read GFF3 file with DANTE annotation of conserved domains of transposases,
    # keep only records with Subclass_1
    tir_domains = dt.get_tir_records_from_dante(args.gff3)

    # export tir_domains to file
    dt.save_gff3_dict_to_file(tir_domains, F'{args.working_dir}/tir_domains.gff3')

    # Extract and save amino acid sequences from DANTE records
    print("\nExtracting amino acid sequences from DANTE domains...", end="")
    aa_sequences = dt.extract_aa_sequences_from_dante(tir_domains)
    aa_fasta_files = dt.save_aa_sequences_by_superfamily(aa_sequences, args.working_dir)
    print(" done")

    # Cluster amino acid sequences using mmseqs2
    print("Clustering amino acid sequences with mmseqs2...", end="")
    aa_mmseqs2_output_dirs = dt.cluster_aa_sequences_mmseqs2(
        aa_fasta_files,
        args.working_dir,
        num_threads=args.cpu
    )
    print(" done")

    # Group sequences based on clustering if max_class_size is specified
    aa_sequence_groups = None
    if args.max_class_size:
        print("Grouping sequences based on mmseqs2 clustering...", end="")
        aa_sequence_groups = dt.group_sequences_by_clusters(
            aa_mmseqs2_output_dirs,
            args.max_class_size
        )
        # Save grouping information for debugging
        dt.save_grouping_info(aa_sequence_groups, aa_mmseqs2_output_dirs)
        print(" done")

    # Read FASTA file with genome assembly
    genome = dt.fasta_to_dict(args.fasta)

    # extract sequence for each protein domain including upstream and downstream
    # use strand information to extract upstream and downstream, if strand is
    # negative, then upstream is downstream and vice versa, for negative
    # reverse complement the sequence
    downstream_seq, upstream_seq, coords = dt.extract_flanking_regions(genome,
                                                                       tir_domains)
    # export downstream and upstream sequences to files, one file for each upstream,
    # downstream and classification

    for cls in downstream_seq:
        class_name = cls.replace('/', '_').replace('|', '_')
        dt.save_fasta_dict_to_file(downstream_seq[cls],
                                   F'{args.working_dir}/{class_name}_downstream_regions.fasta')
        # make blast database from downstream sequences
        dt.make_blast_db(F'{args.working_dir}/{class_name}_downstream_regions.fasta')
        dt.save_fasta_dict_to_file(upstream_seq[cls],
                                   F'{args.working_dir}/'
                                   F'{class_name}_upstream_regions.fasta')
        # make blast database from upstream sequences
        dt.make_blast_db(F'{args.working_dir}/{class_name}_upstream_regions.fasta')


    # export coordinaes with header to file
    dt.save_coords_to_file(coords, F'{args.working_dir}/tir_flank_coords.txt')
    # export downstream and upstream sequences to file
    dt.save_fasta_dict_to_file(downstream_seq, F'{args.working_dir}/downstream_seq.fasta')
    dt.save_fasta_dict_to_file(upstream_seq, F'{args.working_dir}/upstream_seq_rc.fasta')

    # Split original sequences BEFORE fragmentation (if needed)
    if args.max_class_size:
        if aa_sequence_groups:
            # Use clustering-based splitting
            print(f"Splitting classes based on AA domain clustering...", end="")
            split_upstream, split_downstream, split_mapping = dt.split_sequences_by_clustering(
                upstream_seq, downstream_seq, aa_sequence_groups
            )
        else:
            # Fallback: random splitting by threshold
            print(f"Splitting large classes (threshold: {args.max_class_size} sequences)...", end="")
            split_upstream, split_downstream, split_mapping = dt.split_sequences_by_class(
                upstream_seq, downstream_seq, args.max_class_size
            )
        print(" done")
    else:
        # No splitting: wrap in part structure for uniform processing
        split_upstream = {cls: {1: upstream_seq[cls]} for cls in upstream_seq}
        split_downstream = {cls: {1: downstream_seq[cls]} for cls in downstream_seq}
        split_mapping = {cls: [1] for cls in upstream_seq}
    # Fragment each part separately and save to files
    frg_names_upstream = {}
    frg_names_downstream = {}

    for cls in split_upstream:
        frg_names_upstream[cls] = {}
        frg_names_downstream[cls] = {}

        for part_num in split_mapping[cls]:
            # Fragment this part
            up_frags = dt.dict_fasta_to_dict_fragments(split_upstream[cls][part_num])
            down_frags = dt.dict_fasta_to_dict_fragments(split_downstream[cls][part_num])

            # Save fragmented parts
            class_name = cls.replace('/', '_').replace('|', '_')
            up_file = f'{args.working_dir}/{class_name}_upstream.part_{part_num:03d}.fasta'
            down_file = f'{args.working_dir}/{class_name}_downstream.part_{part_num:03d}.fasta'

            dt.save_fasta_dict_to_file(up_frags, up_file)
            dt.save_fasta_dict_to_file(down_frags, down_file)

            frg_names_upstream[cls][part_num] = up_file
            frg_names_downstream[cls][part_num] = down_file

    # Build flat lists of files for parallel processing, tracking class/part mapping
    frgs_fasta_both = []
    frgs_class_mapping = []  # Track which class and part each file pair belongs to

    for cls in frg_names_downstream:
        for part_num in frg_names_downstream[cls]:
            up_file = frg_names_upstream[cls][part_num]
            down_file = frg_names_downstream[cls][part_num]

            frgs_fasta_both.append(up_file)
            frgs_fasta_both.append(down_file)
            frgs_class_mapping.append((cls, part_num, 'upstream'))
            frgs_class_mapping.append((cls, part_num, 'downstream'))

    print("Assembling TIR boundaries...", end="")
    with Pool(processes=args.cpu) as pool:
        aln = pool.map(dt.cap3assembly, frgs_fasta_both, chunksize=1)

    print(" done")

    # Log checksums of input/output fragment files for reproducibility verification (debug only)
    if args.debug:
        checksum_log = F'{args.working_dir}/assembly_input_checksums.txt'
        with open(checksum_log, 'w') as f:
            f.write("# Checksums of input fragment files for CAP3 assembly\n")
            f.write("# Use to verify reproducibility across runs\n")
            for idx, fasta_file in enumerate(frgs_fasta_both):
                checksum = compute_file_checksum(fasta_file)
                mapping_info = frgs_class_mapping[idx] if idx < len(frgs_class_mapping) else ("unknown", "unknown")
                f.write(f"{fasta_file}\t{checksum}\t{mapping_info[0]}\t{mapping_info[1]}\n")

        # Log checksums of CAP3 output files
        cap3_checksum_log = F'{args.working_dir}/assembly_output_checksums.txt'
        with open(cap3_checksum_log, 'w') as f:
            f.write("# Checksums of CAP3 assembly output files\n")
            f.write("# Use to verify reproducibility across runs\n")
            for idx, aln_file in enumerate(aln):
                if aln_file and os.path.exists(aln_file):
                    checksum = compute_file_checksum(aln_file)
                    mapping_info = frgs_class_mapping[idx] if idx < len(frgs_class_mapping) else ("unknown", "unknown")
                    f.write(f"{aln_file}\t{checksum}\t{mapping_info[0]}\t{mapping_info[1]}\n")

    # Reconstruct class-based mapping from flat results
    # Now aln contains alternating upstream/downstream results
    # We need to group them by class, and if split, collect multiple files per class
    ctg_upstream = {}
    ctg_downstream = {}

    aln_idx = 0
    for cls in frg_names_downstream:
        # frg_names now contains dicts with part numbers as keys
        downstream_parts_dict = frg_names_downstream[cls]
        upstream_parts_dict = frg_names_upstream[cls]

        # Collect all aln files for this class
        aln_upstream_parts = []
        aln_downstream_parts = []
        upstream_fasta_parts = []
        downstream_fasta_parts = []

        # Iterate through all parts for this class
        for part_num in sorted(downstream_parts_dict.keys()):
            up_file = upstream_parts_dict[part_num]
            down_file = downstream_parts_dict[part_num]

            aln_upstream_parts.append(aln[aln_idx])
            aln_idx += 1
            aln_downstream_parts.append(aln[aln_idx])
            aln_idx += 1
            upstream_fasta_parts.append(up_file)
            downstream_fasta_parts.append(down_file)

        print("Parsing assemblies for", cls, "...", end="")
        # Pass lists to parse_cap3_aln (it handles both single files and lists)
        ctg_upstream[cls] = dt.parse_cap3_aln(aln_upstream_parts, upstream_fasta_parts, ncpus=args.cpu)
        ctg_downstream[cls] = dt.parse_cap3_aln(aln_downstream_parts, downstream_fasta_parts, ncpus=args.cpu)
        print("done")

    for cls in ctg_upstream:
        prefix = os.path.join(
            args.working_dir,
            cls.replace('/', '_').replace('|', '_')
            )
        for ctg_name in ctg_upstream[cls]:
            # Include part number in filename to avoid collisions from split assemblies
            contig_obj = ctg_upstream[cls][ctg_name]
            if contig_obj.part_num is not None:
                filename = prefix + '_upstream_{}_part{:03d}.fasta'.format(ctg_name, contig_obj.part_num)
            else:
                filename = prefix + '_upstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(contig_obj.alignment,
                                       filename, uppercase=False)

        for ctg_name in ctg_downstream[cls]:
            # Include part number in filename to avoid collisions from split assemblies
            contig_obj = ctg_downstream[cls][ctg_name]
            if contig_obj.part_num is not None:
                filename = prefix + '_downstream_{}_part{:03d}.fasta'.format(ctg_name, contig_obj.part_num)
            else:
                filename = prefix + '_downstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(contig_obj.alignment,
                                       filename, uppercase=False)

    # find TIRs in contigs using R script
    cmd = (F'{script_dir}/detect_tirs.R --contig_dir {args.working_dir} --output '
           F'{args.working_dir} --threads {args.cpu} '
           F'--genome {args.fasta} --seed {args.seed} '
           F'--n_beast_iter {args.n_beast_iter}')
    subprocess.check_call(cmd, shell=True)

    # copy output files to output directory
    file_to_copy = ["DANTE_TIR*",
                    "TIR_classification_summary.txt",
                    ]
    for f in file_to_copy:
        flist = glob.glob(F'{args.working_dir}/{f}')
        for file in flist:
            shutil.copy(file, args.output_dir)
    # remove working directory if debug is not set
    if not args.debug:
        shutil.rmtree(args.working_dir)

    # add version information to gff3 files copied to output directory
    gff3_files = glob.glob(F'{args.output_dir}/*.gff3')
    # first line is header - keep it
    # add second line with version information
    # then rest of the file
    for gff3_file in gff3_files:
        with open(gff3_file, 'r') as f:
            lines = f.readlines()
        with open(gff3_file, 'w') as f:
            f.write(lines[0])
            f.write(F'##DANTE_TIR version {__version__}\n')
            f.writelines(lines[1:])

if __name__ == '__main__':
    main()
