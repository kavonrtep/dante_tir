#!/usr/bin/env python
import math
import random
import subprocess
import os
import shutil
import tempfile
import re
from typing import List, Dict, Tuple, Set
import concurrent.futures
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


def save_fasta_dict_to_file(fasta_dict, fasta_file, uppercase=True):
    """
    save fasta dictionary to file, it will handle nested dictionaries
    :param fasta_dict: dictionary with fasta sequences
    :param fasta_file: path to fasta file
    :return:
    """

    with open(fasta_file, 'w') as f:
        for k, v in fasta_dict.items():
            # check is v is string
            if isinstance(v, str):
                if uppercase:
                    f.write(">{}\n{}\n".format(k, v.upper()))
                else:
                    f.write(">{}\n{}\n".format(k, v))
            else:
                for k2, v2 in v.items():
                    new_id = F'{k2} {k}'
                    if uppercase:
                        f.write(">{}\n{}\n".format(new_id, v2.upper()))
                    else:
                        f.write(">{}\n{}\n".format(new_id, v2))

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

# note : original settings step=70, length=100, jitter=10
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
        seq_len = len(seq)

        for i in range(0, seq_len, step):
            if i > 0:
                ii = i + random.randint(-jitter, jitter)
            else:
                ii = i

            # Ensure offset is within valid bounds
            ii = max(0, min(ii, seq_len - length))

            # Only create fragment if it has content
            if ii + length <= seq_len:
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
    if os.path.exists(stdout_file):
        print(F"File {stdout_file} already exists, skipping assembly")
        return stdout_file
    stderr_file = F'{fasta_file}.cap.err'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    with open(stdout_file, 'w') as f:
        f.write(stdout.decode())
    with open(stderr_file, 'w') as f:
        f.write(stderr.decode())
    return stdout_file


def cap3assembly_with_split(fasta_file, split_threshold=None):
    """
    Run CAP3 assembly with optional splitting for large files.

    If the FASTA file contains more sequences than split_threshold, it will be
    split into multiple parts, each assembled separately. The results are then
    merged into a single output.

    :param fasta_file: Path to FASTA file
    :param split_threshold: Maximum sequences per CAP3 run. If None, no splitting.
    :return: List of CAP3 output files (one element if no split, multiple if split)
    """
    if split_threshold is None:
        # No splitting, just run CAP3 normally
        return [cap3assembly(fasta_file)]

    # Count sequences in file
    fasta_dict = fasta_to_dict(fasta_file)
    num_sequences = len(fasta_dict)

    if num_sequences <= split_threshold:
        # Below threshold, no splitting needed
        return [cap3assembly(fasta_file)]

    # Need to split the file
    seq_ids = sorted(fasta_dict.keys())
    output_files = []
    part_num = 1

    for i in range(0, num_sequences, split_threshold):
        part_seq_ids = seq_ids[i:i + split_threshold]

        # Create part-specific filename
        base_file = fasta_file.rsplit('.', 1)[0] if '.' in fasta_file else fasta_file
        part_file = F'{base_file}.part_{part_num:03d}.fasta'

        # Extract and save sequences for this part
        part_dict = {seq_id: fasta_dict[seq_id] for seq_id in part_seq_ids}
        save_fasta_dict_to_file(part_dict, part_file)

        # Run CAP3 on this part
        output_files.append(cap3assembly(part_file))
        part_num += 1

    return output_files


# Global variable to hold the FASTA dictionary in each worker.
FASTA_DATA = {}

def init_worker(fasta_data):
    """Initializer for worker processes: sets the global FASTA_DATA variable."""
    global FASTA_DATA
    FASTA_DATA = fasta_data

def process_contig_helper(args):
    """
    Processes alignment adjustment for one contig.

    :param args: Tuple (contig, aln_data, orient_data, end_gaps)
    :return: (contig, trimmed) where trimmed is the adjusted alignment for the contig.
    """
    contig, aln_data, orient_data, end_gaps = args
    # 1. Join segments for each read and add end gaps.
    joined = {read: ''.join(aln_data[read]) + end_gaps for read in aln_data}
    max_width = max(len(seq) for seq in joined.values())
    # 2. Pad sequences to have uniform width.
    padded = {
        read: (seq if len(seq) == max_width else seq + ('-' * (max_width - len(seq)))) for
        read, seq in joined.items()}
    # 3. Adjust masked regions using the preloaded FASTA data.
    adjusted = add_masked_reads_to_alignment(padded, orient_data, FASTA_DATA)
    # 4. Determine first and last non-gap columns.
    cols = list(zip(*adjusted.values()))
    aln_start = 0
    aln_end = max_width
    for i, col in enumerate(cols):
        if any(char != '-' for char in col):
            aln_start = i
            break
    for i, col in enumerate(reversed(cols)):
        if any(char != '-' for char in col):
            aln_end = max_width - i
            break
    # 5. Trim each read’s sequence.
    trimmed = {read: seq[aln_start:aln_end] for read, seq in adjusted.items()}
    return contig, trimmed


def parse_cap3_aln(aln_file, input_fasta_file: str, ncpus: int = 1):
    """
    Parse cap3 alignment file(s).

    Handles both single alignment file and multiple alignment files (from splitting).
    When multiple files are provided, contigs from all parts are merged into a single
    output dictionary.

    :param aln_file: Path to cap3 alignment file, or list of paths for split assemblies.
    :param input_fasta_file: Path to input FASTA file used for assembly, or list of paths for splits.
    :param ncpus: Number of CPUs to use for processing.
    :return: Dictionary with contigs (merged from all parts if split).
    """
    # Handle both single file and list of files
    if isinstance(aln_file, str):
        aln_files = [aln_file]
    else:
        aln_files = aln_file

    if isinstance(input_fasta_file, str):
        fasta_files = [input_fasta_file]
    else:
        fasta_files = input_fasta_file

    # Verify file counts match
    if len(aln_files) != len(fasta_files):
        raise ValueError(
            f"Number of alignment files ({len(aln_files)}) does not match "
            f"number of FASTA files ({len(fasta_files)})"
        )

    # Process each part separately, then merge
    all_reads = {}
    all_orientations = {}
    all_alignments = {}
    all_fasta_data = {}

    for aln_file_part, fasta_file_part in zip(aln_files, fasta_files):
        reads, orientations, alignments = _parse_single_cap3_aln(aln_file_part)

        # Merge with all previous parts
        all_reads.update(reads)
        all_orientations.update(orientations)
        all_alignments.update(alignments)

        # Load FASTA data for this part
        fasta_data = fasta_to_dict(fasta_file_part)
        all_fasta_data.update(fasta_data)

    # Filter contigs based on criteria (applied after merging all parts)
    n_min_elements = 4
    contigs_to_exclude = []
    for contig in all_alignments:
        n_elements = len({read.split('_')[0] for read in all_alignments[contig]})
        N_plus = sum(
            1 for read in all_alignments[contig] if all_orientations[contig][read] == '+'
            )
        N_minus = len(all_orientations[contig]) - N_plus
        if (max(N_plus, N_minus) / (N_plus + N_minus) < 0.8) or (
                n_elements < n_min_elements):
            contigs_to_exclude.append(contig)
    for contig in contigs_to_exclude:
        del all_alignments[contig]
        del all_reads[contig]
        del all_orientations[contig]

    # Create argument tuples for per-contig processing
    args_list = [(contig, all_alignments[contig], all_orientations[contig], "-" * 100) for contig
        in all_alignments]
    adjusted_alignments = {}

    # Use ProcessPoolExecutor with initializer to pass the FASTA data to workers.
    with concurrent.futures.ProcessPoolExecutor(
            initializer=init_worker, initargs=(all_fasta_data,)
            ) as executor:
        futures = {executor.submit(process_contig_helper, args): args[0] for args in
                   args_list}
        for future in concurrent.futures.as_completed(futures):
            contig, trimmed = future.result()
            adjusted_alignments[contig] = trimmed

    contigs = {}
    for contig in all_reads:
        contigs[contig] = Contig(
            all_reads[contig], all_orientations[contig], adjusted_alignments[contig]
            )
    return contigs


def _parse_single_cap3_aln(aln_file: str) -> Tuple[Dict, Dict, Dict]:
    """
    Parse a single cap3 alignment file.

    Helper function that parses one alignment file and returns raw data structures
    (before filtering). Used by parse_cap3_aln() to handle multiple split files.

    :param aln_file: Path to cap3 alignment file.
    :return: Tuple of (reads, orientations, alignments) dictionaries.
    """
    # Precompile regex to capture read name and orientation.
    read_line_regex = re.compile(r'^(\d+_\d+_\d+)([+-])')
    reads = {}
    orientations = {}
    alignments = {}
    header = True
    gaps = "-" * 60
    end_gaps = "-" * 100

    contig_name = None
    segment_number = None

    with open(aln_file, 'r') as file:
        for line in file:
            if header:
                if line.startswith('*******************'):
                    parts = line.split()
                    contig_name = parts[1] + parts[2]
                    reads[contig_name] = []
                    orientations[contig_name] = {}
                else:
                    match = read_line_regex.match(line.strip())
                    if match:
                        read_name = match.group(1)
                        orientation = match.group(2)
                        reads[contig_name].append(read_name)
                        orientations[contig_name][read_name] = orientation
                    elif line.startswith('DETAILED DISPLAY OF CONTIGS'):
                        header = False
                        # Initialize alignment data structure.
                        for contig in reads:
                            alignments[contig] = {}
                            for read in reads[contig]:
                                # Each alignment starts with the ending gaps.
                                alignments[contig][read] = [end_gaps]
            else:
                if line.startswith('*******************'):
                    parts = line.split()
                    contig_name = parts[1] + parts[2]
                    segment_number = 1  # Reset segment number (index 0 is starting gaps)
                    for read in reads[contig_name]:
                        alignments[contig_name][read].append(gaps)
                elif read_line_regex.match(line):
                    parts = line.strip().split()
                    # Remove the orientation character.
                    read_name = parts[0][:-1]
                    # Extract alignment segment.
                    seq_segment = line[22:].replace(' ', '-').strip()
                    if len(alignments[contig_name][read_name]) <= segment_number:
                        alignments[contig_name][read_name].append(seq_segment)
                    else:
                        alignments[contig_name][read_name][segment_number] = seq_segment
                elif line.startswith('consensus'):
                    segment_number += 1
                    for read in reads[contig_name]:
                        alignments[contig_name][read].append(gaps)

    return reads, orientations, alignments


def first_letter_index(s):
    match = re.search(r'[a-zA-Z]', s)
    if match:
        return match.start()
    else:
        return -1


def replace_substring(s, start, end, replacement):
    return s[:start] + replacement + s[end:]

def add_masked_reads_to_alignment(alignment: Dict[str, str],
                                  orientations: Dict[str, str],
                                  reads: Dict[str, str]) -> Dict[str, str]:
    """
    add masked reads to alignment
    :param alignment: alignment dictionary
    :param fasta_file: path to fasta file with masked reads
    :return: alignment dictionary with masked reads
    """
    #reads = fasta_to_dict(fasta_file)

    letter_pattern = re.compile(r'[a-zA-Z]')
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

def gff3_quality_ok(gff: Gff3Feature, max_interuptions=3.0, min_identity=0.35,
                    min_similarity=0.45, max_aln_prop=1.2, min_aln_prop=0.8) -> bool:
    """
    Check if quality of gff3 record pass the criteria
    :param gff3_file: path to gff3 file
    :param max_interruptions: maximum number of interruptions in alignment per 100 nt
    :param min_identity: minimum identity
    :param min_similarity: minimum similarity
    :param max_aln_prop: maximum proportion of alignment
    :param min_aln_prop: minimum proportion of alignment
    :return: True if quality is OK
    """
    c1 = min_identity <= float(gff.attributes_dict['Identity'])
    c2 = min_similarity <= float(gff.attributes_dict['Similarity'])
    c3 = max_interuptions >= float(gff.attributes_dict['Relat_Interruptions'])
    c4 = max_aln_prop >= float(gff.attributes_dict['Relat_Length'])
    c5 = min_aln_prop <= float(gff.attributes_dict['Relat_Length'])
    return c1 and c2 and c3 and c4 and c5


def get_tir_records_from_dante(gff3_file: str) -> Dict[str, List[Gff3Feature]]:
    tir_domains = {}
    id = 0
    with open(gff3_file, 'r') as f:
        for line in f:
            # line cannot be empty, must be valid gff3 line
            comment = line[0] == '#'
            empty = len(line.strip()) == 0
            if not comment and not empty:
                gff = Gff3Feature(line)
                cls = gff.attributes_dict['Final_Classification']
                if 'Subclass_1' in cls  and gff3_quality_ok(gff):
                    id += 1
                    if cls not in tir_domains:
                        tir_domains[cls] = []
                    # add also unique id to the attributes
                    gff.attributes_dict['ID'] = id
                    tir_domains[cls].append(gff)
    # print statistics - counts for each classification
    print("\nNumber of protein domain for each superfamily:")
    for cls in tir_domains:
        print(cls, len(tir_domains[cls]))
    print("")
    return tir_domains
def save_gff3_dict_to_file(gff3_dict_list: Dict[str, List[Gff3Feature]], gff3_file: str):
    """
    save list of Gff3Feature to file
    :param gff3_list: list of Gff3Feature
    :param gff3_file: path to gff3 file
    :return:
    """
    with open(gff3_file, 'w') as f:
        # write header
        f.write("##gff-version 3\n")
        for cls in gff3_dict_list:
            for gff in gff3_dict_list[cls]:
                f.write(str(gff))

def extract_flanking_regions(genome: Dict[str, str],
                             tir_domains: Dict[str, List[Gff3Feature]]
                             ) -> Tuple[
    Dict[str, Dict[int, str]],
    Dict[str, Dict[int, str]],
    Dict[str, Dict[str, List[int]]]
]:

    """
    Extract flanking regions of TIR domains
    :param genome:
    :param tir_domains:
    :return: downstream and upstream sequences for each TIR domain
    :return: table with coordinates of donwstream and upstream sequences relative genome

    If sequence is on the negative strand,  sequences are reverse complemented
    """
    upstream_seq = {}
    downstream_seq = {}
    coords = {}
    offset = 6000
    offset2 = 300
    for cls in tir_domains:
        coords[cls] = {}
        upstream_seq[cls] = {}
        downstream_seq[cls] = {}
        for gff in tir_domains[cls]:
            if gff.strand == '+':
                u_s = max(0, gff.start - offset)
                upstream = genome[gff.seqid][u_s:gff.start + offset2]
                downstream = genome[gff.seqid][gff.end - offset2:gff.end + offset]
                # coordinates could be out of range - adjust to sequence length
                offset_u_adjusted = len(upstream) - offset2
                offset_d_adjusted = len(downstream) - offset2
                coords[cls][gff.attributes_dict['ID']] = [
                    gff.seqid,
                    gff.start - offset_u_adjusted,
                    gff.start + offset2,
                    gff.end - offset2,
                    gff.end + offset_d_adjusted,
                    gff.strand]

            else:
                d_s = max(0, gff.start - offset)
                upstream = reverse_complement(
                        genome[gff.seqid][gff.end - offset2:gff.end + offset]
                        )
                downstream = reverse_complement(
                        genome[gff.seqid][d_s:gff.start + offset2]
                        )
                # coordinates could be out of range - adjust to sequence length
                offset_u_adjusted = len(upstream) - offset2
                offset_d_adjusted = len(downstream) - offset2

                coords[cls][gff.attributes_dict['ID']] = [
                    gff.seqid,
                    gff.end - offset2,
                    gff.end + offset_u_adjusted,
                    gff.start - offset_d_adjusted,
                    gff.start + offset2,
                    gff.strand]

            ID = gff.attributes_dict['ID']
            upstream_seq[cls][ID] = upstream
            # to make analysis easier, downstream sequence is reverse complemented - insertions site will be on 5' end
            downstream_seq[cls][ID] = reverse_complement(downstream)

    return downstream_seq, upstream_seq, coords

def split_sequences_by_clustering(upstream_seq: Dict[str, Dict[int, str]],
                                  downstream_seq: Dict[str, Dict[int, str]],
                                  aa_sequence_groups: Dict[str, Dict[int, list]]) -> Tuple[
    Dict[str, Dict[int, Dict[int, str]]],
    Dict[str, Dict[int, Dict[int, str]]],
    Dict[str, List[int]]
]:
    """
    Split upstream and downstream sequences based on AA domain clustering groups.

    Uses the clustering-based grouping from mmseqs2 to partition sequences.
    Each group becomes a separate processing part, respecting cluster boundaries.

    :param upstream_seq: Dict[class][ID] → sequence
    :param downstream_seq: Dict[class][ID] → sequence
    :param aa_sequence_groups: Dict[class][group_id] → [list of sequence_ids] from clustering
    :return: Tuple of (split_upstream, split_downstream, split_mapping)
        - split_upstream: Dict[class][part_num][ID] → sequence
        - split_downstream: Dict[class][part_num][ID] → sequence
        - split_mapping: Dict[class] → list of part_nums
    """
    split_upstream = {}
    split_downstream = {}
    split_mapping = {}

    for cls in upstream_seq:
        split_upstream[cls] = {}
        split_downstream[cls] = {}
        split_mapping[cls] = []

        # Get clustering groups for this class if available
        if cls in aa_sequence_groups and aa_sequence_groups[cls]:
            # Use clustering-based groups
            groups = aa_sequence_groups[cls]

            # Sort group IDs to process in order
            for group_id in sorted(groups.keys()):
                group_member_ids = groups[group_id]

                # Convert member IDs to integers for matching with upstream_seq keys
                part_ids = []
                for member_id in group_member_ids:
                    # member_id might be string, need to convert
                    try:
                        seq_id = int(member_id)
                        if seq_id in upstream_seq[cls]:
                            part_ids.append(seq_id)
                    except (ValueError, TypeError):
                        # Skip if cannot convert
                        continue

                # Create part with paired sequences from this group
                if part_ids:
                    part_num = len(split_mapping[cls]) + 1
                    split_upstream[cls][part_num] = {
                        seq_id: upstream_seq[cls][seq_id] for seq_id in part_ids
                    }
                    split_downstream[cls][part_num] = {
                        seq_id: downstream_seq[cls][seq_id] for seq_id in part_ids
                    }
                    split_mapping[cls].append(part_num)
        else:
            # Fallback: no clustering groups, treat entire class as one part
            upstream_ids = sorted(upstream_seq[cls].keys())
            split_upstream[cls][1] = upstream_seq[cls].copy()
            split_downstream[cls][1] = downstream_seq[cls].copy()
            split_mapping[cls] = [1]

    return split_upstream, split_downstream, split_mapping


def split_sequences_by_class(upstream_seq: Dict[str, Dict[int, str]],
                             downstream_seq: Dict[str, Dict[int, str]],
                             threshold: int) -> Tuple[
    Dict[str, Dict[int, Dict[int, str]]],
    Dict[str, Dict[int, Dict[int, str]]],
    Dict[str, List[int]]
]:
    """
    Split upstream and downstream sequences into parts while maintaining pairing.

    For each class, if the number of sequences exceeds threshold, split into
    multiple parts. Each part contains paired upstream/downstream sequences from
    the SAME original IDs, ensuring no ID mismatches.

    :param upstream_seq: Dict[class][ID] → sequence
    :param downstream_seq: Dict[class][ID] → sequence
    :param threshold: Maximum sequences per part (max_class_size)
    :return: Tuple of (split_upstream, split_downstream, split_mapping)
        - split_upstream: Dict[class][part_num][ID] → sequence
        - split_downstream: Dict[class][part_num][ID] → sequence
        - split_mapping: Dict[class] → list of part_nums
    """
    split_upstream = {}
    split_downstream = {}
    split_mapping = {}

    for cls in upstream_seq:
        # Get sorted IDs for consistent ordering
        upstream_ids = sorted(upstream_seq[cls].keys())

        if len(upstream_ids) <= threshold:
            # No split needed - wrap in part structure for uniform processing
            split_upstream[cls] = {1: upstream_seq[cls].copy()}
            split_downstream[cls] = {1: downstream_seq[cls].copy()}
            split_mapping[cls] = [1]
        else:
            # Split into multiple parts
            split_upstream[cls] = {}
            split_downstream[cls] = {}
            split_mapping[cls] = []
            part_num = 1

            for i in range(0, len(upstream_ids), threshold):
                part_ids = upstream_ids[i:i + threshold]

                # Create part with paired sequences
                split_upstream[cls][part_num] = {
                    seq_id: upstream_seq[cls][seq_id] for seq_id in part_ids
                }
                split_downstream[cls][part_num] = {
                    seq_id: downstream_seq[cls][seq_id] for seq_id in part_ids
                }
                split_mapping[cls].append(part_num)
                part_num += 1

    return split_upstream, split_downstream, split_mapping

def save_coords_to_file(coords: Dict[str, Dict[str, List[int]]], coords_file: str):
    """
    save coordinates to file
    :param coords:
    :param coords_file:
    :return:
    """
    with open(coords_file, 'w') as f:
        # make header
        f.write("SeqID\tClass\tID\tUpstream_start\tUpstream_end\tDownstream_start"
                "\tDownstream_end\tStrand\n")
        for cls in coords:
            for ID in coords[cls]:
                f.write(
                    F'{coords[cls][ID][0]}\t{cls}\t{ID}\t{coords[cls][ID][1]}\t'
                    F'{coords[cls][ID][2]}\t{coords[cls][ID][3]}\t{coords[cls][ID][4]}\t'
                    F'{coords[cls][ID][5]}\n'
                    )

class Contig:
    """
    Class to store contig sequences as alignment, reads and orientations
    """

    def __init__(self, reads: List[str], orientations: Dict[str, str],
        alignment: Dict[str, str]):
        """
        :type alignments: Dict[str, str]
        :type orientations: Dict[str, str]
        :type reads: List[str]
        """
        self.reads = reads
        self.orientations = orientations
        self.alignment = alignment
        self.alignment_total_length = len(alignment[reads[0]])
        self._pssm = None
        self._coverage = None
        self._consensus = None
        self._masked_proportion = None
        self._element_start_end = None
        self._mean_coverage = None
        self._information_content = None

    @property
    def pssm(self) -> List[Dict[str, int]]:
        """
        Calculate position specific frequencies of nucleotides
        :return: position specific frequencies
        """
        if self._pssm is None:
            self._pssm = self.calculate_pssm()
        return self._pssm

    @property
    def coverage(self) -> List[int]:
        """
        coverage at each position of the alignment excluding gaps
        :return: coverage
        """
        if self._coverage is None:
            self._coverage = self.calculate_coverage()
        return self._coverage

    @property
    def masked_proportion(self) -> List[float]:
        """
        Proportion of masked nucleotides at each position of the alignment
        :return:masksed proportion
        """
        if self._masked_proportion is None:
            self._masked_proportion = self.calculate_proportion_of_masked_bases()
        return self._masked_proportion


    def calculate_pssm(self) -> List[Dict[str, int]]:
        """
        Calculate position specific frequences of nucleotides
        :return: position specific frequences
        """

        pssm = []
        for i in range(self.alignment_total_length):
            pssm.append({'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0})
        for read_name in self.alignment:
            for i in range(self.alignment_total_length):
                # skip if N
                if self.alignment[read_name][i].upper() == 'N':
                    continue
                pssm[i][self.alignment[read_name][i].upper()] += 1
        return pssm

    def calculate_coverage(self) -> List[int]:
        """
        Calculate coverage by nucleotide along the alignment, excluding gaps
        :return: coverage
        """
        coverage = [0] * self.alignment_total_length
        # calculate coverage from pssm
        for i in range(self.alignment_total_length):
            coverage[i] = sum(self.pssm[i].values()) - self.pssm[i]['-']
        return coverage

    def calculate_proportion_of_masked_bases(self) -> List[float]:
        """
        Calculate proportion of masked bases at each position of
        the alignment
        Masked bases are in lower case
        :return: proportion of masked bases
        """
        masked_bases = [0] * self.alignment_total_length
        for i in range(self.alignment_total_length):
            for read_name in self.alignment:
                if self.alignment[read_name][i].islower():
                    masked_bases[i] += 1
        # recalculate to proportion using coverage
        for i in range(self.alignment_total_length):
            if self.coverage[i] > 0:
                masked_bases[i] = masked_bases[i] / self.coverage[i]
            else:
                masked_bases[i] = -1
        return masked_bases

    @property
    def information_content(self) -> List[float]:
        """
        return information content at each position of the alignment
        :return: List[float]
        """
        if self._information_content is None:
            self._information_content = self.calculate_information_content()
        return self._information_content


    #recalculater pssf frequency to information content
    def calculate_information_content(self) -> List[float]:
        """
        Calculate information content at each position of the alignment
        ignoring gaps
        :return: information content
        """
        information_content = [0] * self.alignment_total_length
        for i in range(self.alignment_total_length):
            if self.coverage[i] > 0:
                for base in ['A', 'C', 'G', 'T']:
                    if self.pssm[i][base] > 0:
                        p = self.pssm[i][base] / self.coverage[i]
                        information_content[i] += p * math.log2(p)
                information_content[i] = 2 + information_content[i]
            else:
                information_content[i] =  0
        return information_content
    @property
    def consensus(self) -> str:
        """
        Calculate consensus sequence of the alignment
        :return: consensus sequence
        """
        if self._consensus is None:
            self._consensus = self.calculate_consensus_sequence()
        return self._consensus

    def calculate_consensus_sequence(self) -> str:
        """
        Calculate 50% consensus sequence of the alignment,
        :return: consensus sequence
        """
        consensus = []
        for i in range(self.alignment_total_length):
            max_base = '-'
            max_count = 0
            if self.coverage[i] == 0:
                consensus += '-'
                continue
            consensus += 'N'
            for base in ['A', 'C', 'G', 'T']:
                if self.pssm[i][base]/self.coverage[i] > 0.5:
                    consensus[i] = base
                    break
        return "".join(consensus)

    @property
    def element_start_end(self) -> Set[str]:
        """
        Calculate list of elements in the contig
        :return: list of elements
        """
        if self._element_start_end is None:
            elements_start = {}
            elements_end = {}
            element_list = []
            for read in self.reads:
                elements_id, start, end = read.split('_')
                if elements_id not in elements_start:
                    elements_start[elements_id] = int(start)
                    elements_end[elements_id] = int(end)
                else:
                    elements_start[elements_id] = min(elements_start[elements_id], int(start))
                    elements_end[elements_id] = max(elements_end[elements_id], int(end))

            # get keys from dictionary sortef by values (reverse)
            sorted_keys = sorted(elements_start, key=elements_start.get, reverse=True)
            for i in sorted_keys:
                element_list.append([i, elements_start[i], elements_end[i]])
            self._element_start_end = element_list
        return self._element_start_end

    @property
    def mean_coverage(self) -> float:
        """
        Calculate mean coverage of the contig, excluding gaps
        :return: mean coverage
        """
        if self._mean_coverage is None:
            self._mean_coverage = sum(self.coverage) / len(self.coverage)
        return self._mean_coverage


    def find_left_insertion_sites(self, debug=False) -> int:
        """
        Find insertion sites in the contig
        :return: number of insertion sites
        """

        ws = 10
        max_iter = int(min(100, len(self.coverage)/2 + ws))
        C1 = 5 # minimum coverage in masked part
        C2 = 5 # minimum coverage in TIR part
        left_insertion_site = -1
        # check masked region exists
        for i in range(max_iter):
            up_mask = sum(self.masked_proportion[0:i])/(i+1)
            up_coverage = sum(self.coverage[0:i])/(i+1)
            up_bits = sum(self.information_content[0:i])/(i+1)
            tir_0_bits = self.information_content[i]
            tir_window_bits = sum(self.information_content[i:i + ws])/ws
            tir_coverage = sum(self.coverage[i:i + ws])/ws
            tir_mask = sum(self.masked_proportion[i:i + ws])/ws
            cond1 = up_mask > 0.5 and up_coverage > C1 and up_bits < 1
            cond2 = tir_0_bits == 2 and tir_window_bits > 1.8 and tir_coverage > C2 and\
                    tir_mask < 0.1
            if debug:
                print(i, up_mask, up_coverage, up_bits, tir_0_bits, tir_window_bits,
                      tir_coverage, tir_mask, cond1, cond2)

            if cond1 and cond2:
                left_insertion_site = i
                break

        if debug:
            print("coverage")
            print(self.coverage)
            print('masked proportion')
            print(self.masked_proportion)
            print('information content')
            print(self.information_content)
            print('consensus')
            print(self.consensus)

        return left_insertion_site


class BreakIt(Exception): pass


def clean_aa_sequence(aa_seq: str) -> str:
    """
    Remove non-standard characters from amino acid sequence.
    Keeps only standard amino acids (A-Z).
    :param aa_seq: amino acid sequence string
    :return: cleaned amino acid sequence
    """
    # Keep only letters (standard amino acids)
    cleaned = ''.join(c for c in aa_seq if c.isalpha())
    return cleaned


def extract_aa_sequences_from_dante(tir_domains: Dict[str, list]) -> Dict[str, Dict[int, str]]:
    """
    Extract amino acid sequences from DANTE GFF3 records.
    The Region_Seq attribute contains the AA sequence of the domain.
    :param tir_domains: dictionary of tir_domains from get_tir_records_from_dante
    :return: nested dictionary with structure {classification: {id: aa_sequence}}
    """
    aa_sequences = {}

    for cls in tir_domains:
        aa_sequences[cls] = {}
        for gff in tir_domains[cls]:
            record_id = gff.attributes_dict['ID']
            # Get Region_Seq from attributes
            if 'Region_Seq' in gff.attributes_dict:
                aa_seq = gff.attributes_dict['Region_Seq']
                # Clean the sequence (remove non-standard characters)
                cleaned_seq = clean_aa_sequence(aa_seq)
                aa_sequences[cls][record_id] = cleaned_seq

    return aa_sequences


def save_aa_sequences_by_superfamily(aa_sequences: Dict[str, Dict[int, str]],
                                     working_dir: str) -> Dict[str, str]:
    """
    Save amino acid sequences to FASTA files, one file per superfamily.
    IDs match the format used in upstream/downstream_regions.fasta files.
    :param aa_sequences: nested dictionary {classification: {id: aa_sequence}}
    :param working_dir: working directory path
    :return: dictionary mapping classification to AA FASTA file path
    """
    aa_fasta_files = {}

    for cls in aa_sequences:
        class_name = cls.replace('/', '_').replace('|', '_')
        aa_fasta_file = F'{working_dir}/{class_name}_aa_sequences.fasta'

        with open(aa_fasta_file, 'w') as f:
            for record_id, aa_seq in aa_sequences[cls].items():
                f.write(F'>{record_id}\n{aa_seq}\n')

        aa_fasta_files[cls] = aa_fasta_file

    return aa_fasta_files


def cluster_aa_sequences_mmseqs2(aa_fasta_files: Dict[str, str],
                                  working_dir: str,
                                  num_threads: int = 1,
                                  min_seq_id: float = 0.8,
                                  cov_mode: int = 0,
                                  coverage: float = 0.8,
                                  sensitivity: float = 7.0,
                                  e_value: float = 1e-3) -> Dict[str, str]:
    """
    Cluster amino acid sequences using mmseqs2 easy-cluster for each superfamily.
    Each AA FASTA file is clustered independently with default settings.

    :param aa_fasta_files: dictionary mapping classification to AA FASTA file path
    :param working_dir: working directory path
    :param num_threads: number of CPUs to use
    :param min_seq_id: minimum sequence identity for clustering (0.0-1.0) [default: 0.9]
    :param cov_mode: coverage mode
                     0: coverage of query and target
                     1: coverage of target
                     2: coverage of query
                     3: target seq. length >= x% of query length
                     4: query seq. length >= x% of target length
                     5: short seq. >= x% of the other seq. length [default: 0]
    :param coverage: fraction of aligned (covered) residues [default: 0.8]
    :param sensitivity: sensitivity: 1.0 faster; 4.0 fast; 7.5 sensitive [default: 4.0]
    :param e_value: E-value threshold [default: 1e-3]
    :return: dictionary mapping classification to mmseqs2 output directory
    """
    mmseqs_output_dirs = {}

    for cls, aa_fasta_file in aa_fasta_files.items():
        class_name = cls.replace('/', '_').replace('|', '_')
        mmseqs_output_prefix = F'{working_dir}/mmseqs2_aa_{class_name}/clusters'

        # Create output directory for mmseqs2 results
        os.makedirs(os.path.dirname(mmseqs_output_prefix), exist_ok=True)

        # Create a temporary directory for mmseqs2 intermediate files
        tmp_dir = F'{working_dir}/.mmseqs2_tmp_{class_name}'
        os.makedirs(tmp_dir, exist_ok=True)

        try:
            # Run mmseqs2 easy-cluster with default parameters
            cmd = (
                F'mmseqs easy-cluster {aa_fasta_file} {mmseqs_output_prefix} '
                F'{tmp_dir} '
                F'--min-seq-id {min_seq_id} '
                F'--cov-mode {cov_mode} '
                F'-c {coverage} '
                F'-s {sensitivity} '
                F'-e {e_value} '
                F'--threads {num_threads}'
            )

            subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

            # Store the directory containing the output files
            mmseqs_output_dirs[cls] = os.path.dirname(mmseqs_output_prefix)

        except subprocess.CalledProcessError as e:
            print(F'Error clustering {class_name} with mmseqs2: {e}')
            continue
        finally:
            # Clean up temporary directory
            if os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)

    return mmseqs_output_dirs


def parse_mmseqs2_clusters(cluster_tsv_file: str) -> Dict[str, list]:
    """
    Parse mmseqs2 cluster output file (adjacency list format).
    Format: cluster_id (representative) \t member_id (one per line)

    Example:
      7545    7545   (representative 7545 contains member 7545)
      8910    44     (representative 8910 contains member 44)
      8910    133    (representative 8910 contains member 133)

    :param cluster_tsv_file: path to clusters_cluster.tsv file from mmseqs2
    :return: dictionary {cluster_id: [list of member_ids]}
    """
    clusters = {}

    if not os.path.exists(cluster_tsv_file):
        return clusters

    with open(cluster_tsv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                cluster_id = parts[0]  # representative sequence ID (cluster ID)
                member_id = parts[1]   # member sequence ID

                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(member_id)

    return clusters


def group_sequences_by_clusters(mmseqs_output_dirs: Dict[str, str],
                                max_class_size: int) -> Dict[str, Dict[int, list]]:
    """
    Group sequences based on mmseqs2 clustering results, respecting max_class_size.

    Algorithm:
    1. Parse mmseqs2 cluster assignments from clusters_cluster.tsv
    2. Sort clusters by size (largest first)
    3. For large clusters (>= 2x max_class_size): split into equal-sized groups
    4. For smaller clusters: merge together to reach ~max_class_size
    5. Avoid splitting small clusters; instead merge them

    :param mmseqs_output_dirs: dictionary mapping classification to mmseqs2 output directory
    :param max_class_size: maximum number of sequences per group
    :return: dictionary {classification: {group_id: [sequence_ids]}}
    """
    sequence_groups = {}

    for cls in mmseqs_output_dirs:
        mmseqs_dir = mmseqs_output_dirs[cls]
        cluster_file = F'{mmseqs_dir}/clusters_cluster.tsv'

        # Parse cluster assignments
        clusters = parse_mmseqs2_clusters(cluster_file)

        if not clusters:
            # If no clusters found, return empty groups
            sequence_groups[cls] = {}
            continue

        # Convert to list of (cluster_id, [members]) sorted by cluster size (descending)
        cluster_list = [(cluster_id, members) for cluster_id, members in clusters.items()]
        cluster_list.sort(key=lambda x: len(x[1]), reverse=True)

        groups = {}
        current_group = 1
        current_group_members = []

        # First pass: handle large clusters that need splitting
        remaining_clusters = []

        for cluster_id, members in cluster_list:
            cluster_size = len(members)

            # If cluster is >= 2x max_class_size, split it
            if cluster_size >= 2 * max_class_size:
                # Split large cluster into equal-sized groups
                num_subgroups = (cluster_size + max_class_size - 1) // max_class_size
                subgroup_size = cluster_size // num_subgroups

                for i in range(num_subgroups):
                    start_idx = i * subgroup_size
                    if i == num_subgroups - 1:
                        # Last subgroup gets remaining members
                        end_idx = cluster_size
                    else:
                        end_idx = (i + 1) * subgroup_size

                    subgroup = members[start_idx:end_idx]
                    groups[current_group] = subgroup
                    current_group += 1
            else:
                # Keep smaller clusters intact for merging
                remaining_clusters.append((cluster_id, members))

        # Second pass: merge smaller clusters
        for cluster_id, members in remaining_clusters:
            cluster_size = len(members)

            # Check if adding this cluster would exceed max_class_size
            if len(current_group_members) + cluster_size <= max_class_size:
                # Add to current group
                current_group_members.extend(members)
            else:
                # Save current group if it has members and start a new one
                if current_group_members:
                    groups[current_group] = current_group_members
                    current_group += 1
                    current_group_members = []

                # Add cluster to new group
                current_group_members.extend(members)

        # Save the last group if it has members
        if current_group_members:
            groups[current_group] = current_group_members

        sequence_groups[cls] = groups

    return sequence_groups


def save_grouping_info(sequence_groups: Dict[str, Dict[int, list]],
                       mmseqs_output_dirs: Dict[str, str]) -> None:
    """
    Save grouping information to mmseqs2 output directory for debugging.

    :param sequence_groups: dictionary {classification: {group_id: [sequence_ids]}}
    :param mmseqs_output_dirs: dictionary mapping classification to mmseqs2 output directory
    :return: None
    """
    for cls in sequence_groups:
        mmseqs_dir = mmseqs_output_dirs[cls]
        grouping_file = F'{mmseqs_dir}/grouping_info.txt'

        with open(grouping_file, 'w') as f:
            f.write("# Grouping information based on mmseqs2 clustering\n")
            f.write("# Format: Group_ID\tNumber_of_sequences\tSequence_IDs\n\n")

            for group_id in sorted(sequence_groups[cls].keys()):
                members = sequence_groups[cls][group_id]
                member_str = ','.join(str(m) for m in members)
                f.write(F'{group_id}\t{len(members)}\t{member_str}\n')


def make_blast_db(fasta_file: str, dbtype: str = 'nucl'):
    """
    return None
    make blast database from fasta file
    :param fasta_file: path to fasta file
    :param dbtype: type of database
    :return: None
    """
    cmd = F'makeblastdb -in {fasta_file} -dbtype {dbtype}'
    # capture output
    p = subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    return None
