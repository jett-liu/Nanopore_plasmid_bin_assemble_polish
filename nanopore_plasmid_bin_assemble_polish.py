#!/usr/bin/env python3

import os
import gzip
import glob
import argparse
from Bio import SeqIO
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def main(output_folder, input_fastq, num_threads, canu_binary_path):
    # Check if output folder exists, create if not
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Add backslash to output folder path if not present
    if output_folder[-1] != '/':
        output_folder += '/'

    #check if input fastq is gzipped, if not, gzip it
    if not input_fastq.endswith('.gz'):
        os.system(f'gzip {input_fastq}')
        input_fastq = input_fastq + '.gz'

    # Read in fastQ and collect read lengths
    read_lengths = []
    with gzip.open(input_fastq, "rt") as f:
        for record in SeqIO.parse(f, "fastq"):
            read_lengths.append(len(record.seq))

    # Convert read_lengths to numpy array for analysis
    read_lengths = np.array(read_lengths)
    bins = np.array(range(read_lengths.min(), read_lengths.max(), 200))

    # Average and standard deviation of number of reads per bin
    avg_density = np.mean([len(read_lengths[(read_lengths >= bins[i]) & (read_lengths < bins[i+1])]) for i in range(len(bins)-1)])
    std_density = np.std([len(read_lengths[(read_lengths >= bins[i]) & (read_lengths < bins[i+1])]) for i in range(len(bins)-1)])
    
    # three sigma cutoff
    two_sigma = avg_density + 3 * std_density
    
    # Find bins where the number of reads exceeds the two sigma cutoff
    bins_two_sigma = [bins[i] for i in range(len(bins)-1) if len(read_lengths[(read_lengths >= bins[i]) & (read_lengths < bins[i+1])]) > two_sigma]

    # Extract reads for each bin that is greater than two sigma
    for bin in bins_two_sigma:
        with gzip.open(input_fastq, "rt") as f:
            reads = [record for record in SeqIO.parse(f, "fastq") if len(record.seq) >= (bin-200) and len(record.seq) <= (bin+200)]
            with open(f'{output_folder}{bin.round(0)}.fastq', 'w') as f_out:
                SeqIO.write(reads, f_out, "fastq")

    # Plot the histogram of read lengths
    sns.histplot(read_lengths, binwidth=200)
    plt.axhline(y=two_sigma, color='gray', linestyle='--', label='Average + 2*std')
    plt.text(max(read_lengths) + 0.05*max(read_lengths), two_sigma, f'binning cutoff', fontsize=10, color='gray')
    sns.despine()
    plt.savefig(f'{output_folder}read_length_hist.png', bbox_inches='tight', dpi=600)

    # Check the number of reads extracted for each bin
    for file in glob.glob(f'{output_folder}*.fastq'):
        print(f'{file}: {len(list(SeqIO.parse(file, "fastq")))} reads')

    # Run Canu for each extracted bin
    for read_file in glob.glob(f"{output_folder}*.fastq"):
        # Create Canu output folder
        canu_output_folder = f'{output_folder}canu_output_{os.path.basename(read_file).replace(".fastq", "")}'
        if not os.path.exists(canu_output_folder):
            os.makedirs(canu_output_folder)

        basename = os.path.basename(read_file).replace(".fastq", "").split(".")[0]
        genome_size = int(round((int(basename)/1000), 0))

        # Run Canu command
        canu_command = f"{canu_binary_path} -p {basename} -d {canu_output_folder} genomeSize={genome_size}k -nanopore {read_file}"
        os.system(canu_command)

    # Trim and polish the Canu outputs
    trimmed_paths = []
    for canu_output_folder in glob.glob(f'{output_folder}canu_output_*'):
        # Get the contigs file
        contigs_file = glob.glob(f'{canu_output_folder}/*.contigs.fasta')[0]
        basename = os.path.basename(canu_output_folder).replace("canu_output_", "")
        count = 0

        for record in SeqIO.parse(contigs_file, 'fasta'):
            count += 1
            trimmed_contigs_file = f'{canu_output_folder}/{basename}_{count}.trimmed_contigs.fasta'

            if "suggestCircular=yes" in record.description:
                print(f"This contig is circular: {basename}_{count}")
                trim_left = int(record.description.split("trim=")[1].split("-")[0])
                trim_right = int(record.description.split("trim=")[1].split("-")[1])
                trimmed_seq = record.seq[trim_left:trim_right]
                with open(trimmed_contigs_file, 'w') as f_out:
                    f_out.write(">" + f"{basename}_{count}" + "\n")
                    f_out.write(str(trimmed_seq) + "\n")
                trimmed_paths.append(trimmed_contigs_file)
            else:
                print(f"This contig is not circular: {basename}_{count}")
                with open(trimmed_contigs_file, 'w') as f_out:
                    f_out.write(">" + f"{basename}_{count}" + "\n")
                    f_out.write(str(record.seq) + "\n")
                trimmed_paths.append(trimmed_contigs_file)

    # Run Medaka polishing
    for contig_path in trimmed_paths:
        polished_output_folder = f'{output_folder}polished_output_{os.path.basename(contig_path).replace(".fasta", "")}'
        medaka_command = f"medaka_consensus -i {input_fastq} -d {contig_path} -o {polished_output_folder} -t {num_threads} --bacteria"
        os.system(medaka_command)

if __name__ == '__main__':
    # Setup argument parsing
    parser = argparse.ArgumentParser(description="Process Nanopore reads and polish contigs")
    parser.add_argument('--output_folder','-o', required=True, help="Output folder path")
    parser.add_argument('--input_fastq','-i', required=True, help="Input FastQ file path")
    parser.add_argument('--num_threads','-t', required=False, default=1, type=int, help="Number of threads (default: 1)")
    parser.add_argument('--canu_binary_path','-c', required=True, help="Path to Canu binary")

    args = parser.parse_args()

    # Call the main function with arguments
    main(args.output_folder, args.input_fastq, args.num_threads, args.canu_binary_path)
