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
    tir_domains = {}
    with open(args.gff3, 'r') as f:
        for line in f:
            if line[0] != '#':
                gff = dt.Gff3Feature(line)
                cls = gff.attributes_dict['Final_Classification']
                if 'Subclass_1' in cls:
                    if cls not in tir_domains:
                        tir_domains[cls] = []
                    tir_domains[cls].append(gff)
    # print statistics - counts for each classification
    print("Number of protein domain for each superfamily:")
    for cls in tir_domains:
        print(cls, len(tir_domains[cls]))

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
    for cls in tir_domains:
        upstream_seq[cls] = {}
        downstream_seq[cls] = {}
        for gff in tir_domains[cls]:
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
            upstream_seq[cls][id] = upstream
            downstream_seq[cls][id] = downstream
            gff.attributes_dict['ID'] = id
            gff_dict[id] = gff # to better access gff records.

    # write sequences to file in output directory
    # create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # make fragments and save to file
    frg_names_upstream = {}
    frg_names_downstream = {}
    for cls in tir_domains:
        # sanitize classification name so it can be used as a file name
        prefix = os.path.join(
            args.output_dir,
            cls.replace('/', '_').replace('|', '_')
            )
        frg_names_upstream[cls] = prefix + '_upstream.fasta'
        frg_names_downstream[cls] = prefix + '_downstream.fasta'

        fragments = dt.fragment_fasta_dict(upstream_seq[cls])
        dt.save_fasta_dict_to_file(fragments, frg_names_upstream[cls])

        fragments = dt.fragment_fasta_dict(downstream_seq[cls])
        dt.save_fasta_dict_to_file(fragments, frg_names_downstream[cls])


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
        print(cls, x, y)
        print("-----------------------------")
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



if __name__ == '__main__':
    main()
