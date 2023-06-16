#!/usr/bin/env python
import random
import subprocess
import os
import tempfile
import re
from typing import List, Dict, Tuple
class Gff3Feature:
    """
    Class for gff3 feature
    """

    def __init__(self, line):
        self.line = line
        self.items = line.strip().split('\t')
        self.seqid = self.items[0]
        self.source = self.items[1]
        self.type = self.items[2]
        self.start = int(self.items[3])
        self.end = int(self.items[4])
        self.score = self.items[5]
        self.strand = self.items[6]
        self.frame = self.items[7]
        self.attributes = self.items[8]
        self._attributes_dict = {}
        for item in self.attributes.split(';'):
            if item != '':
                try:
                    key, value = item.split('=')
                    self._attributes_dict[key] = value
                except ValueError:
                    print(item)
                    print(self.attributes)
                    raise

        self._attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self._attributes_dict.items()]
            )

    @property
    def attributes_dict(self):
        """
        store attributes as dictionary
        :return:
        """
        return self._attributes_dict

    @property
    def attributes_str(self):
        """
        store attributes as string
        :return:
        """
        return self._attributes_str

    @attributes_str.getter
    def attributes_str(self):
        """
        store attributes as string
         :return:
        """
        self._attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self._attributes_dict.items()]
            )
        return self._attributes_str

    @attributes_dict.setter
    def attributes_dict(self, value):
        self._attributes_dict = value
        self._attributes_str = ';'.join(
            ['{}={}'.format(key, value) for key, value in self._attributes_dict.items()]
            )

    def __str__(self):
        return '\t'.join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __repr__(self):
        return '\t'.join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __eq__(self, other):
        return self.line_recalculated() == other.line_recalculated()

    def __hash__(self):
        return hash(self.line_recalculated())

    def get_line(self):
        """returns original line"""
        return self.line

    def overlap(self, other):
        """
        Check if two features overlap
        :param other:
        :return:
        """
        if self.start <= other.end and self.end >= other.start:
            return True
        else:
            return False

    def line_recalculated(self):
        """
        :return:
        string with recalculated line
        """
        return '\t'.join(
            [self.seqid, self.source, self.type, str(self.start), str(self.end),
             self.score, self.strand, self.frame, self.attributes_str]
            ) + '\n'

    def __lt__(self, other):
        width = self.end - self.start
        other_width = other.end - other.start
        return width < other_width

    def __gt__(self, other):
        width = self.end - self.start
        other_width = other.end - other.start
        return width > other_width

    def identical_region(self, other):
        """
        Check if two features are identical
        :param other:
        :return:
        """
        if self.start == other.start and self.end == other.end and self.seqid == \
                other.seqid:
            return True
        else:
            return False

    def print_line(self):
        """
        :return:
        string with recalculated line
        """
        columns = [self.seqid, self.source, self.type, str(self.start), str(self.end),
                   self.score, self.strand, self.frame]
        attributes_list = ['{}={}'.format(key, value) for key, value in
                           self.attributes_dict.items()]
        attributes = [';'.join(attributes_list)]
        return '\t'.join(columns + attributes) + '\n'


def read_fasta_sequence_size(fasta_file):
    """Read size of sequence into dictionary"""
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line[0] == '>':
                header = line.strip().split(' ')[0][1:]  # remove part of name after space
                fasta_dict[header] = 0
            else:
                fasta_dict[header] += len(line.strip())
    return fasta_dict


def fasta_to_dict(fasta_file):
    """
    convert fasta file to dictionary
    :param fasta_file: path to fasta file
    :return: dictionary with fasta sequences
    """
    fasta_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:

            if line.startswith(">"):
                seq_name = line.split()[0].replace(">", "")
                fasta_dict[seq_name] = [""]  # initialize list
            else:
                fasta_dict[seq_name].append(line.strip())
    # concatenate list in dictionary
    fasta_dict = {k: "".join(v) for k, v in fasta_dict.items()}
    return fasta_dict


def get_seq_from_fasta(fasta_dict, seq_id, start=None, end=None):
    """
    get sequence from fasta dictionary
    :param fasta_dict: dictionary with fasta sequences
    :param seq_id: sequence ID
    :param start: start position
    :param end: end position
    :return: sequence
    """
    if start and end:
        return fasta_dict[seq_id][start:end]
    return fasta_dict[seq_id]


def gff3_to_fasta(gff3_file, fasta_file, additonal_attribute=None):
    """
    extract fasta sequences from gff3 file
    it is generator, returns one sequence at time and seq ID plus additional attribute
    if provided
    :param additonal_attribute: yield additional attribute from gff3 file
    :param gff3_file: path to gff3 file
    :param fasta_file: path to fasta file
    :return:
    """

    fasta_dict = fasta_to_dict(fasta_file)
    with open(gff3_file, 'r') as f1:
        for line in f1:
            if line.startswith("#"):
                continue
            gff3_feature: Gff3Feature = Gff3Feature(line)
            s = get_seq_from_fasta(
                fasta_dict, gff3_feature.seqid, gff3_feature.start, gff3_feature.end
                )
            if "ID" not in gff3_feature.attributes_dict:
                gff3_feature.attributes_dict["ID"] = (gff3_feature.seqid + "_" + str(
                    gff3_feature.start
                    ) + "_" + str(gff3_feature.end))
            if additonal_attribute:
                yield [gff3_feature.attributes_dict['ID'], s,
                       gff3_feature.attributes_dict[additonal_attribute]]
            else:
                yield [gff3_feature.attributes_dict['ID'], s]


def save_fasta_dict_to_file(fasta_dict, fasta_file):
    """
    save fasta dictionary to file
    :param fasta_dict: dictionary with fasta sequences
    :param fasta_file: path to fasta file
    :return:
    """
    with open(fasta_file, 'w') as f:
        for k, v in fasta_dict.items():
            f.write(">{}\n{}\n".format(k, v.upper()))


def reverse_complement(dna: str) -> str:
    """
    reverse complement of dna sequence, including ambiguous bases
    :param dna:
    :return: reverse complement of dna sequence
    """
    # complement including all IUPAC codes
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}

    return ''.join(complement[base] for base in reversed(dna.upper()))


def dict_fasta_to_dict_fragments(fasta_dict: Dict[str, str],
                                 step=70, length=100, jitter=10) -> Dict[str, str]:
    """
    fragment fasta sequences to overlapping parts
    :param fasta_dict: dictionary with fasta sequences
    :param step: step size
    :param length: length of fragment
    :param jitter: random jitter size - make step size random withing jitter limits
    :return: dictionary with fragments
    """
    fragments_dict = {}
    for seq_id, seq in fasta_dict.items():
        for i in range(0, len(seq), step):
            if i > 0:
                ii = i + random.randint(-jitter, jitter)
            else:
                ii = i
            id = F'{seq_id}_{ii}_{ii + length}'
            fragments_dict[id] = seq[ii:ii + length]
    return fragments_dict


def make_fragment_files(output_dir: str,
                        downstream_seq: Dict[str, Dict[ int, str]],
                        upstream_seq: Dict[str, Dict[ int, str]]
                        ) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    make fragment files for downstream and upstream sequences
    :param output_dir:
    :param downstream_seq:
    :param upstream_seq:
    :return: dictionaries with fragment file names
    """
    frg_names_upstream = {}
    frg_names_downstream = {}
    for cls in upstream_seq:
        # sanitize classification name so it can be used as a file name
        prefix = os.path.join(
                output_dir,
                cls.replace('/', '_').replace('|', '_')
                )
        frg_names_upstream[cls] = prefix + '_upstream.fasta'
        frg_names_downstream[cls] = prefix + '_downstream.fasta'

        fragments = dict_fasta_to_dict_fragments(upstream_seq[cls])
        save_fasta_dict_to_file(fragments, frg_names_upstream[cls])

        fragments = dict_fasta_to_dict_fragments(downstream_seq[cls])
        save_fasta_dict_to_file(fragments, frg_names_downstream[cls])
    return frg_names_downstream, frg_names_upstream



def cap3assembly(fasta_file):
    """
    run cap3 assembly
    :param fasta_file: path to fasta file
    :return: path to cap3 output file
    assume that cap3 is in PATH
    """
    cmd = F'cap3 {fasta_file} -o 40 -p 80 -x cap -w'
    stdout_file = F'{fasta_file}.cap.aln'
    stderr_file = F'{fasta_file}.cap.err'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    with open(stdout_file, 'w') as f:
        f.write(stdout.decode())
    with open(stderr_file, 'w') as f:
        f.write(stderr.decode())
    return stdout_file


def parse_cap3_aln(filename):
    """
    parse cap3 alignment file
    :param cap3_aln_file: path to cap3 alignment file
    :return: dictionary with contigs
    """
    results = []
    with open(filename, 'r') as file:
        contig_name = ''
        reads = {}
        orientations = {}
        header = True
        alignments = {}
        gaps = "-" * 60
        # end gaps are added to beginging and and of all alignments - this will simplify
        # later adding masked parts of reads which would be otherview out ot range
        end_gaps = "-" * 100
        for line in file:
            if header:
                if line.startswith('*******************'):
                    contig_name = line.split()[1] + line.split()[2]
                    reads[contig_name] = []
                    orientations[contig_name] = {}
                elif re.match(r'\d+_\d+_\d+[+-]', line.strip()):
                    parts = line.strip().split()
                    read_name = parts[0][:-1]
                    orientation = parts[0][-1]
                    reads[contig_name].append(read_name)
                    orientations[contig_name][read_name] = orientation
                elif line.startswith('DETAILED DISPLAY OF CONTIGS'):
                    header = False  # end of part 1, breaking to start part 2
                    # create empty dictionary for alignments, with sequence names as keys
                    for contig in reads:
                        contig_name = ''
                        alignments[contig] = {}
                        for read in reads[contig]:
                            alignments[contig][read] = [end_gaps]  # starting gaps of aln
            else:
                if line.startswith('*******************'):
                    contig_name = line.split()[1] + line.split()[2]
                    positions = 0
                    segment_number = 1  # starting with 1, 0 is for starting gaps
                    # fill in gaps for all reads that are asignment to the contig
                    for read in reads[contig_name]:
                        alignments[contig_name][read].append(gaps)
                elif re.match(r'\d+_\d+_\d+[+-]', line):
                    parts = line.strip().split()
                    read_name = parts[0][:-1]
                    seq = line[22:].replace(' ',"-").strip()
                    alignments[contig_name][read_name][segment_number] = seq
                elif line.startswith('consensus'):
                    segment_number += 1
                    positions += 60
                    # fill in gaps for all reads that are asignment to the contig
                    for read in reads[contig_name]:
                        alignments[contig_name][read].append(gaps)

    # concatenate segments for each read
    for contig in alignments:
        for read in alignments[contig]:
            alignments[contig][read] = ''.join(alignments[contig][read]) + end_gaps

    asm = Assembly(reads, orientations, alignments)
    return asm

def first_letter_index(s):
    match = re.search(r'[a-zA-Z]', s)
    if match:
        return match.start()
    else:
        return -1


def replace_substring(s, start, end, replacement):
    return s[:start] + replacement + s[end:]

def add_masked_reads_to_alignment(alignment, orientations, fasta_file):
    """
    add masked reads to alignment
    :param alignment: alignment dictionary
    :param fasta_file: path to fasta file with masked reads
    :return: alignment dictionary with masked reads
    """
    reads = fasta_to_dict(fasta_file)

    letter_pattern = re.compile(r'[a-zA-Z]')
    sum_plus_orientation = sum([1 for x in orientations.values() if x == '+'])
    orientation_ratio = sum_plus_orientation / len(orientations)
    print(F'orientation ratio: {orientation_ratio}')
    print(F'number of reads: {len(orientations)}')
    print("-----------------------------------")
    for read_name in alignment:
        # remove all gaps to get read length
        read_aligned = alignment[read_name].replace('-', '')
        aligned_length = len(alignment[read_name].replace('-', ''))
        read_length_orig = len(reads[read_name])
        # test if read is not fully aligned
        if aligned_length < read_length_orig:
            # use regular expression to get firt leter afte "-"
            matches = [m for m in letter_pattern.finditer(alignment[read_name])]
            start_pos = matches[0].start()
            end_pos = matches[-1].start()
            # get masked part of the read
            # first from beggining
            if orientations[read_name] == '+':
                r = reads[read_name]
            else:
                r = reverse_complement(reads[read_name])
            s_index = r.find(read_aligned)
            e_index = s_index + aligned_length
            left_part = r[:s_index]
            right_part = r[e_index:]
            # add masked part to alignment
            new_aln = replace_substring(alignment[read_name],
                              start_pos - len(left_part),
                              start_pos,
                              left_part.lower())
            new_aln = replace_substring(new_aln,
                                end_pos + 1,
                                end_pos + len(right_part) + 1,
                                right_part.lower())

            alignment[read_name] = new_aln
    return alignment

def get_tir_records_from_dante(gff3_file: str) -> Dict[str, List[Gff3Feature]]:
    tir_domains = {}
    id = 0
    with open(gff3_file, 'r') as f:
        for line in f:
            if line[0] != '#':
                gff = Gff3Feature(line)
                cls = gff.attributes_dict['Final_Classification']
                if 'Subclass_1' in cls:
                    id += 1
                    if cls not in tir_domains:
                        tir_domains[cls] = []
                    # add also unique id to the attributes
                    gff.attributes_dict['ID'] = id
                    tir_domains[cls].append(gff)
    # print statistics - counts for each classification
    print("Number of protein domain for each superfamily:")
    for cls in tir_domains:
        print(cls, len(tir_domains[cls]))
    return tir_domains


def extract_flanking_regions(genome: Dict[str, str],
                             tir_domains: Dict[str, List[Gff3Feature]]
                             ) -> Tuple[
    Dict[str, Dict[int, str]],
    Dict[str, Dict[int, str]]]:

    """
    Extract flanking regions of TIR domains
    :param genome:
    :param tir_domains:
    :return: downstream and upstream sequences for each TIR domain

    If sequence is on the negative strand,  sequences are reverse complemented
    """
    upstream_seq = {}
    downstream_seq = {}
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
                upstream = reverse_complement(
                        genome[gff.seqid][gff.end:gff.end + offset]
                        )
                downstream = reverse_complement(
                        genome[gff.seqid][gff.start - offset:gff.start]
                        )
            ID = gff.attributes_dict['ID']
            upstream_seq[cls][ID] = upstream
            downstream_seq[cls][ID] = downstream

    return downstream_seq, upstream_seq




class Assembly:
    """
    Class to store assembly alignments, reads and orientations
    """

    def __init__(self, reads, orientations, alignments):
        self.reads = reads
        self.orientations = orientations
        self.alignments = alignments