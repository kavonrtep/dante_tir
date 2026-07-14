#!/usr/bin/env python3
"""Equivalence test for the memory-frugal extract_flanking_regions rewrite.

extract_flanking_regions used to take a whole-genome dict from fasta_to_dict
(~genome_bp of RAM, OOM on large assemblies). It now takes the FASTA path and
streams one sequence at a time (peak = largest single sequence). This test
proves the streaming version returns EXACTLY the same upstream/downstream/coords
as the original dict-based logic, across:
  - both strands,
  - windows clamped at the sequence start (start < offset) and end
    (end + offset > length),
  - a line-wrapped multi-sequence FASTA,
  - domains on more than one sequence.

Run: python3 tests/test_extract_flanking_regions.py
"""
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import dt_utils as dt  # noqa: E402


def reference_extract(genome, tir_domains):
    """Verbatim copy of the ORIGINAL dict-based implementation (the reference)."""
    upstream_seq, downstream_seq, coords = {}, {}, {}
    offset, offset2 = 6000, 300
    for cls in tir_domains:
        coords[cls], upstream_seq[cls], downstream_seq[cls] = {}, {}, {}
        for gff in tir_domains[cls]:
            if gff.strand == '+':
                u_s = max(0, gff.start - offset)
                upstream = genome[gff.seqid][u_s:gff.start + offset2]
                downstream = genome[gff.seqid][gff.end - offset2:gff.end + offset]
                offset_u_adjusted = len(upstream) - offset2
                offset_d_adjusted = len(downstream) - offset2
                coords[cls][gff.attributes_dict['ID']] = [
                    gff.seqid, gff.start - offset_u_adjusted, gff.start + offset2,
                    gff.end - offset2, gff.end + offset_d_adjusted, gff.strand]
            else:
                d_s = max(0, gff.start - offset)
                upstream = dt.reverse_complement(
                    genome[gff.seqid][gff.end - offset2:gff.end + offset])
                downstream = dt.reverse_complement(
                    genome[gff.seqid][d_s:gff.start + offset2])
                offset_u_adjusted = len(upstream) - offset2
                offset_d_adjusted = len(downstream) - offset2
                coords[cls][gff.attributes_dict['ID']] = [
                    gff.seqid, gff.end - offset2, gff.end + offset_u_adjusted,
                    gff.start - offset_d_adjusted, gff.start + offset2, gff.strand]
            ID = gff.attributes_dict['ID']
            upstream_seq[cls][ID] = upstream
            downstream_seq[cls][ID] = dt.reverse_complement(downstream)
    return downstream_seq, upstream_seq, coords


def make_seq(n, salt=0):
    return "".join("ACGT"[(i * 7 + salt) % 4] for i in range(n))


def write_wrapped_fasta(path, seqs, width=70):
    with open(path, "w") as fh:
        for name, seq in seqs:
            fh.write(">%s some description here\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def gff(seqid, start, end, strand, fid, cls):
    line = "\t".join([seqid, "dante_tir", "tir", str(start), str(end),
                      ".", strand, ".", "ID=%s;Class=%s" % (fid, cls)])
    return dt.Gff3Feature(line)


def main():
    seqs = [("chr1", make_seq(9000, 1)), ("chr2", make_seq(1000, 2))]
    tmp = tempfile.mkdtemp(prefix="tir_flank_test_")
    fasta = os.path.join(tmp, "genome.fasta")
    write_wrapped_fasta(fasta, seqs)

    tir_domains = {
        "hAT": [
            gff("chr1", 7000, 7400, "+", "d_mid_plus", "Class_II|Subclass_1|TIR|hAT"),
            gff("chr1", 500, 900, "-", "d_mid_minus", "Class_II|Subclass_1|TIR|hAT"),
            gff("chr1", 100, 400, "+", "d_start_clamp", "Class_II|Subclass_1|TIR|hAT"),
        ],
        "MuDR_Mutator": [
            gff("chr2", 600, 950, "-", "d_end_clamp", "Class_II|Subclass_1|TIR|MuDR_Mutator"),
            gff("chr2", 300, 500, "+", "d_chr2_plus", "Class_II|Subclass_1|TIR|MuDR_Mutator"),
        ],
    }

    genome = dt.fasta_to_dict(fasta)
    ref_down, ref_up, ref_coords = reference_extract(genome, tir_domains)
    new_down, new_up, new_coords = dt.extract_flanking_regions(fasta, tir_domains)

    assert new_down == ref_down, "downstream sequences differ from reference"
    assert new_up == ref_up, "upstream sequences differ from reference"
    assert new_coords == ref_coords, "coords differ from reference"

    n = sum(len(v) for v in ref_up.values())
    print("  extract_flanking_regions: %d domains, streaming output == dict output" % n)

    # sanity: the generator reconstructs wrapped sequences exactly
    got = dict(dt.fasta_record_generator(fasta))
    for name, seq in seqs:
        assert got[name] == seq, "fasta_record_generator mangled %s" % name
    print("  fasta_record_generator: wrapped multi-seq FASTA reconstructed exactly")
    print("test_extract_flanking_regions: PASSED")


if __name__ == "__main__":
    main()
