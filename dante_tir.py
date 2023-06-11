#!/usr/bin/env python3
"""
Script for identification DNA transpsons with Terminal Inverted Repeats (TIRs) basen
DANTE annotation of conserved domains
"""

import argparse
import dt_utils as dt

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Program  for identification DNA '
        'transposons with Terminal Inverted '
        'Repeats (TIRs) basen DANTE annotation of conserved domains of transposases')
    parser.add_argument(
        '-g', '--gff3', help='GFF3 file with DANTE annotation of conserved domains of transposases', required=True)
    parser.add_argument(
        '-f', '--fasta', help='FASTA file with genome assembly', required=True)
    parser.add_argument(
        '-o', '--output', help='Output file with TIRs', required=True)

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
    # reversecomplement the sequence
    for classification in tir_domains:
        for gff in tir_domains[classification]:
            if gff.strand == '+':
                upstream = genome[gff.seqid][gff.start - 1000:gff.start]
                downstream = genome[gff.seqid][gff.end:gff.end + 1000]
            else:
                upstream = genome[gff.seqid][gff.end:gff.end + 1000].reverse_complement()
                downstream = genome[gff.seqid][gff.start - 1000:gff.start].reverse_complement()
            gff.attributes_dict['upstream'] = upstream
            gff.attributes_dict['downstream'] = downstream


if __name__ == '__main__':
    main()


