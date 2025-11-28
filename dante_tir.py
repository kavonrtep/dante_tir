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

    print("--------------------------------------------------------")
    print("")
    print("          ---|>>>---- DANTE_TIR---- <<<|---")
    print("")
    print("Domain Based identification of DNA transposons with TIRs")
    print("")
    print("--------------------------------------------------------")




    args = parser.parse_args()

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

    # make fragments of upstream and downstream sequences and save to file
    frg_names_downstream, frg_names_upstream = dt.make_fragment_files(
        args.working_dir, downstream_seq, upstream_seq
        )

    # assembly fragments into contigs
    # Handle optional splitting of large classes
    if args.max_class_size:
        print(f"Splitting large classes (threshold: {args.max_class_size} sequences)...", end="")
        # Split upstream and downstream files while keeping them paired
        split_upstream = {}
        split_downstream = {}
        for cls in frg_names_upstream:
            up_parts, down_parts = dt.split_fasta_coordinated(
                frg_names_upstream[cls],
                frg_names_downstream[cls],
                args.max_class_size
            )
            split_upstream[cls] = up_parts
            split_downstream[cls] = down_parts
        frg_names_upstream = split_upstream
        frg_names_downstream = split_downstream
        print(" done")

    # Build flat lists of files for parallel processing, tracking class/part mapping
    frgs_fasta_both = []
    frgs_class_mapping = []  # Track which class each file pair belongs to
    # file are converted to list to be able to use zip function
    # pass list to map function
    for cls in frg_names_downstream:
        downstream_parts = frg_names_downstream[cls] if isinstance(frg_names_downstream[cls], list) else [frg_names_downstream[cls]]
        upstream_parts = frg_names_upstream[cls] if isinstance(frg_names_upstream[cls], list) else [frg_names_upstream[cls]]

        for up_file, down_file in zip(upstream_parts, downstream_parts):
            frgs_fasta_both.append(up_file)
            frgs_fasta_both.append(down_file)
            frgs_class_mapping.append((cls, 'upstream'))
            frgs_class_mapping.append((cls, 'downstream'))

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
        downstream_parts = frg_names_downstream[cls] if isinstance(frg_names_downstream[cls], list) else [frg_names_downstream[cls]]
        upstream_parts = frg_names_upstream[cls] if isinstance(frg_names_upstream[cls], list) else [frg_names_upstream[cls]]

        # Collect all aln files for this class
        aln_upstream_parts = []
        aln_downstream_parts = []
        upstream_fasta_parts = []
        downstream_fasta_parts = []

        for up_file, down_file in zip(upstream_parts, downstream_parts):
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
            filename = prefix + '_upstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(ctg_upstream[cls][ctg_name].alignment,
                                       filename, uppercase=False)

        for ctg_name in ctg_downstream[cls]:
            filename = prefix + '_downstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(ctg_downstream[cls][ctg_name].alignment,
                                       filename, uppercase=False)

    # find TIRs in contigs using R script
    cmd = (F'{script_dir}/detect_tirs.R --contig_dir {args.working_dir} --output '
           F'{args.working_dir} --threads {args.cpu} '
           F'--genome {args.fasta} --seed {args.seed}')
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
