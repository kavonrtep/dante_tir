#!/usr/bin/env python3
"""
Script for identification DNA transpsons with Terminal Inverted Repeats (TIRs) basen
DANTE annotation of conserved domains
"""

import argparse
import os
from itertools import chain

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

    args = parser.parse_args()

    # Read GFF3 file with DANTE annotation of conserved domains of transposases,
    # keep only records with Subclass_1
    tir_domains = dt.get_tir_records_from_dante(args.gff3)

    # Read FASTA file with genome assembly
    genome = dt.fasta_to_dict(args.fasta)

    # extract sequence for each protein domain including upstream and downstream
    # use strand information to extract upstream and downstream, if strand is
    # negative, then upstream is downstream and vice versa, for neagetive
    # reverse complement the sequence
    downstream_seq, upstream_seq = dt.extract_flanking_regions(genome, tir_domains)

    # write sequences to file in output directory
    # create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)


    # make fragments of upstream and downstream sequences and save to file
    frg_names_downstream, frg_names_upstream = dt.make_fragment_files(
        args.output_dir, downstream_seq, upstream_seq
        )

    # assembly fragments into contigs
    frgs_fasta_upstream = []
    frgs_fasta_downstream = []
    for cls in tir_domains:
        frgs_fasta_downstream.append(frg_names_downstream[cls])
        frgs_fasta_upstream.append(frg_names_upstream[cls])
    frgs_fasta_both = list(chain.from_iterable(zip(frgs_fasta_upstream, frgs_fasta_downstream)))

    frgs_fasta = list(frg_names_upstream.values()) + list(frg_names_downstream.values())
    print(frgs_fasta)
    with Pool(processes=3) as pool:
        aln = pool.map(dt.cap3assembly, frgs_fasta_both)

    aln_upstream = aln[::2]
    aln_downstream = aln[1::2]


    ctg_upstream = {}
    ctg_downstream = {}
    for cls, x, y in zip(tir_domains, aln_upstream, aln_downstream):
        ctg_upstream[cls] = dt.parse_cap3_aln(x)
        ctg_downstream[cls] = dt.parse_cap3_aln(y)

    # write contigs to file as multifasta
    for cls in ctg_upstream:
        prefix = os.path.join(
            args.output_dir,
            cls.replace('/', '_').replace('|', '_')
            )
        for ctg_name in ctg_upstream[cls].alignments:
            filename = prefix + '_upstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(ctg_upstream[cls].alignments[ctg_name], filename)
        for ctg_name in ctg_downstream[cls].alignments:
            filename = prefix + '_downstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(ctg_downstream[cls].alignments[ctg_name], filename)

    # adjust alignments:
    for cls, fup, fdown in zip(ctg_upstream,frgs_fasta_upstream, frgs_fasta_downstream):
        prefix = os.path.join(
                args.output_dir,
                cls.replace('/', '_').replace('|', '_')
                )
        for ctg_name in ctg_upstream[cls].alignments:
            aln = ctg_upstream[cls].alignments[ctg_name]
            orientations = ctg_upstream[cls].orientations[ctg_name]
            aln = dt.add_masked_reads_to_alignment(aln,orientations,  fup)
            # export adjusted alignment to file
            filename = prefix + '_upstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(aln, filename)
        for ctg_name in ctg_downstream[cls].alignments:
            aln = ctg_downstream[cls].alignments[ctg_name]
            orientations = ctg_downstream[cls].orientations[ctg_name]
            aln = dt.add_masked_reads_to_alignment(aln, orientations, fdown)
            # export adjusted alignment to file
            filename = prefix + '_downstream_' + ctg_name + '.fasta'
            dt.save_fasta_dict_to_file(aln, filename)




if __name__ == '__main__':
    main()
