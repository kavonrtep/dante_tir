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
from version import __version__
from time import time

import dt_utils as dt
from multiprocessing import Pool

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
        '--version', action='version',
        version='%(prog)s {version}'.format(version=__version__)
        )
    # add debug option
    parser.add_argument(
        '--debug', action='store_true',
        help='Debug mode', default=False
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
    frgs_fasta_upstream = []
    frgs_fasta_downstream = []
    # file are converted to list to be able to use zip function
    # pass list to map function
    for cls in frg_names_downstream:
        frgs_fasta_downstream.append(frg_names_downstream[cls])
        frgs_fasta_upstream.append(frg_names_upstream[cls])
    frgs_fasta_both = list(chain.from_iterable(zip(frgs_fasta_upstream, frgs_fasta_downstream)))

    #frgs_fasta = list(frg_names_upstream.values()) + list(frg_names_downstream.values())
    print("Assembling TIR boundaries...", end="")
    with Pool(processes=args.cpu) as pool:
        aln = pool.map(dt.cap3assembly, frgs_fasta_both, chunksize=1)

    print(" done")
    aln_upstream = aln[::2]
    aln_downstream = aln[1::2]

    ctg_upstream = {}
    ctg_downstream = {}
    print("Parsing assemblies")
    for cls, x, y in zip(frg_names_downstream, aln_upstream, aln_downstream):
        print("Parsing assemblies for", cls, "...", end="")
        ctg_upstream[cls] = dt.parse_cap3_aln(x, frg_names_upstream[cls], ncpus=args.cpu)
        ctg_downstream[cls] = dt.parse_cap3_aln(y, frg_names_downstream[cls], ncpus=args.cpu)
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
           F'--genome {args.fasta}')
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
