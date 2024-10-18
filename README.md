# Nanopore_plasmid_bin_assemble_polish

Install the conda environment with 
```
conda env create -n Nanopore --file=Nanopore.yaml
```

Packages are listed in 
```
Nanopore.txt
```

Canu also needs to be install from the binary release: https://github.com/marbl/canu/releases

To run the script:
```
conda activate Nanopore
python nanopore_plasmid_bin_assemble_polish.py -i <input_fastq_file> -o <output_folder> --canu_binary_path <path_to_canu_binary> -t <number_of_threads_for_processing>
```

Output:
1) histogram of read_lengths. The cutoff for binning different plasmids is marked by a gray bar.
2) polished assemblies are named "consensus.fasta" in the "polished_output" folders.

Pipeline description:
1) pull reads from each bin in the histogram above the cutoff
2) assemble each bin of reads separately (Canu)
3) polish each assembly with all reads (medaka)

Notes:  
  
Bin width is 200 bp  
Binning cutoff is 3 std above the mean read length  
This seems to break when the number of reads in a bin is < 100  
