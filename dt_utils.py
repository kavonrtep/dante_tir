#!/usr/bin/env python
import tempfile


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
