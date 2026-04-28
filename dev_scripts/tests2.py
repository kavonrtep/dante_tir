#!/usr/bin/env python
'''Test dante_tir.py on various inputs.'''
import os
import subprocess
from time import time

# identify files DANTE.gff3 and genome.fasta list of directories

main_dir_path = '/mnt/ceph/454_data/DToL/Henderson_paper/tests/container_repeat_annot/'

# structure is main_dir_path/species_name/genome.fasta and main_dir_path/species_name/DANTE/DANTE.gff3
# if either of the files is missing, the species is skipped

species_dirs = [d for d in os.listdir(main_dir_path) if os.path.isdir(os.path.join(main_dir_path, d))]
output_dir = "tmp/tests_outputs"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
print(species_dirs)

species_table = F'{output_dir}/species_table.txt'

for species in species_dirs:
    genome_file = os.path.join(main_dir_path, species, 'genome.fasta')
    gff3_file = os.path.join(main_dir_path, species, 'DANTE', 'DANTE.gff3')
    if not os.path.exists(genome_file) or not os.path.exists(gff3_file):
        print(f"Skipping {species} because genome or DANTE file is missing.")
        continue

    print(f"Running dante_tir.py on {species}...")
    # run dante_tir.py on the species
    cmd = f"dante_tir.py -f {genome_file} -g {gff3_file} -o {output_dir}/{species} -c 20"
    print(cmd)
    # we need to capture the output of the command and all errors to log
    # log name is output_dir/species.log
    #  record the time of execution
    start_time = time()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    with open(F'{output_dir}/{species}.log', 'w') as log:
        log.write(out.decode())
        log.write(err.decode())
    end_time = time()
    time_of_execution = (end_time - start_time)/60
    # get exit code of the command
    exit_code = p.returncode
    # append species name, exit code, and output directory to species_table
    with open(species_table, 'a') as st:
        st.write(f"{species}\t{exit_code}\t{output_dir}/{species}\t{time_of_execution}\n")




