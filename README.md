# DANTE_TIR

**DANTE_TIR** is a bioinformatics tool designed for the identification of DNA transposable elements containing Terminal Inverted Repeats (TIRs). This tool leverages prior annotations from [**DANTE**](https://github.com/kavonrtep/dante), a separate software that identifies conserved protein domains, specifically transposases.

## Overview

**DANTE_TIR** identifies DNA transposons with Terminal Inverted Repeats (TIRs) based on pre-existing annotations of conserved transposase domains provided by the **DANTE** tool. It extends the annotation workflow by extracting and assembling sequences flanking these conserved domains to discover potential TIR motifs through comparative analysis.




### Process Overview:

1. **Domain Annotation Input:**
   - Utilizes output from **DANTE**, which identifies conserved transposase protein domains in genomic sequences.
2. **Sequence Extraction:**
   - Extracts genomic sequences upstream and downstream from the identified transposase domains.
3. **Contig Assembly:**
   - Upstream and donwstream flankig sequences from each superfamily are assembled into contigs.
4. **TIR Detection:**
   - Analyzes assembled contigs to identify conserved sequences potentially representing TIR motifs.
   - Compares the 5' and 3' ends of contigs to confirm inverted repeat structures.
5. **Target Site Duplication (TSD) Analysis:**
   - Upon TIR confirmation, the tool examines flanking regions to identify potential Target Site Duplications (TSDs).
6. **Iterative Refinement:**
   - Confirmed TIR sequences are reused in subsequent rounds of detection to improve sensitivity.



Each DNA transposon superfamily is analyzed separately, with parameters optimized specifically for the characteristics of each superfamily. Detection sensitivity depends on the abundance of transposon elements within the genome. Low-copy superfamilies may be difficult to detect due to insufficient sequence conservation while sensitivity is maximized for superfamilies that are well-represented in the genome. 

## Installation

Using conda:
```
conda install -c conda-forge -c r -c bioconda -c petrnovak  dante_tir
```

## Usage

```
dante_tir.py -g annotation.gff3 -f genome.fasta -o output_directory [-c number_of_CPUs]
```

### Arguments:
- `-g, --gff3`: GFF3 file with DANTE annotation of conserved domains of transposases (**required**).
- `-f, --fasta`: FASTA file with genome assembly (**required**).
- `-o, --output_dir`: Output directory where TIR results will be stored (**required**).
- `-c, --cpu`: Number of CPUs to use (optional, default = 1).

## Dependencies

- DANTE (required for initial conserved domain annotation)
- (Other dependencies to be added)

## Citation

TODO: Add citation information


