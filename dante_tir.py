#!/usr/bin/env python3
"""
Script for identification DNA transpsons with Terminal Inverted Repeats (TIRs) basen
DANTE annotation of conserved domains
"""

import argparse
import os
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
    tir_domains = {}
    with open(args.gff3, 'r') as f:
        for line in f:
            if line[0] != '#':
                gff = dt.Gff3Feature(line)
                classification = gff.attributes_dict['Final_Classification']
                if 'Subclass_1' in classification:
                    if classification not in tir_domains:
                        tir_domains[classification] = []
                    tir_domains[classification].append(gff)
    # print statistics - counts for each classification
    print("Number of protein domain for each superfamily:")
    for classification in tir_domains:
        print(classification, len(tir_domains[classification]))

    # Read FASTA file with genome assembly
    genome = dt.fasta_to_dict(args.fasta)

    # extract sequence for each protein domain including upstream and downstream
    # use strand information to extract upstream and downstream, if strand is
    # negative, then upstream is downstream and vice versa, for neagetive
    # reverse complement the sequence
    upstream_seq = {}
    downstream_seq = {}
    id = 0
    gff_dict = {}
    offset = 6000
    offset2 = 300
    for classification in tir_domains:
        upstream_seq[classification] = {}
        downstream_seq[classification] = {}
        for gff in tir_domains[classification]:
            if gff.strand == '+':
                upstream = genome[gff.seqid][gff.start - offset:gff.start + offset2]
                downstream = genome[gff.seqid][gff.end - offset2:gff.end + offset]
            else:
                upstream = dt.reverse_complement(
                    genome[gff.seqid][gff.end:gff.end + offset]
                    )
                downstream = dt.reverse_complement(
                    genome[gff.seqid][gff.start - offset:gff.start]
                    )
            id += 1
            upstream_seq[classification][id] = upstream
            downstream_seq[classification][id] = downstream
            gff.attributes_dict['ID'] = id
            gff_dict[id] = gff # to better access gff records.

    # write sequences to file in output directory
    # create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # make fragments and save to file
    frg_names_upstream = {}
    frg_names_downstream = {}
    for classification in tir_domains:
        # sanitize classification name so it can be used as a file name
        prefix = os.path.join(
            args.output_dir,
            classification.replace('/', '_').replace('|', '_')
            )
        frg_names_upstream[classification] = prefix + '_upstream.fasta'
        frg_names_downstream[classification] = prefix + '_downstream.fasta'

        fragments = dt.fragment_fasta_dict(upstream_seq[classification])
        dt.save_fasta_dict_to_file(fragments, frg_names_upstream[classification])

        fragments = dt.fragment_fasta_dict(downstream_seq[classification])
        dt.save_fasta_dict_to_file(fragments, frg_names_downstream[classification])


    # assembly fragments into contigs
    frgs_fasta = list(frg_names_upstream.values()) + list(frg_names_downstream.values())
    print(frgs_fasta)
    with Pool(processes=3) as pool:
        assembly = pool.map(dt.cap3assembly, frgs_fasta)
    for res in assembly:
        print(res)

    dt.parse_cap3_aln(assembly[4])

if __name__ == '__main__':
    main()
