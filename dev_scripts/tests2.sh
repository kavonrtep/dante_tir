#!/bin/bash
# Developer-only HPC harness — not run in CI. See dev_scripts/README.md.
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
echo "SCRIPT_DIR: $SCRIPT_DIR"

if [ -d "$SCRIPT_DIR/tmp/test2" ]; then
    echo "Output directory $SCRIPT_DIR/tmp/test2 already exists. Do you want to delete it? [y/n]"
    read answer
    if [ "$answer" == "y" ]; then
        rm -r "$SCRIPT_DIR/tmp/test2"
    else
        echo "Exiting..."
        exit 1
    fi
fi


"$SCRIPT_DIR/dante_tir.py" --gff3 /mnt/raid/454_data/DToL/Henderson_paper/tests/container_repeat_annot/JunEffu_240903/DANTE/DANTE.gff3 \
--fasta /mnt/raid/454_data/DToL/Henderson_paper/tests/container_repeat_annot/JunEffu_240903/genome.fasta \
--output_dir "$SCRIPT_DIR/tmp/test2" --cpu 10
