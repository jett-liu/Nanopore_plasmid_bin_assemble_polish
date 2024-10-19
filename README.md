# Nanopore_plasmid_bin_assemble_polish

Converting whole plasmid sequencing data to consensus fasta files. Even if multiplexed

## Local installation instructions
Install the conda environment with 
```
conda env create -n Nanopore --file=Nanopore.yaml
```

Packages are listed in 
```
Nanopore.txt
```

Canu also needs to be installed from the binary distribution: https://github.com/marbl/canu/releases

To run the script:
```
conda activate Nanopore
python nanopore_plasmid_bin_assemble_polish.py -i <input_fastq_file> -o <output_folder> --canu_binary_path <path_to_canu_binary> -t <number_of_threads_for_processing>
```

Output:
1) histogram of read lengths. The cutoff for binning different plasmids is marked by a gray bar.
2) polished assemblies are named "consensus.fasta" in the "polished_output" folders.

Pipeline description:
1) pull reads from each bin in the histogram above the cutoff.
2) assemble each bin of reads separately (Canu).
3) polish each assembly with all reads (medaka).

Notes:  
  
Bin width is 200 bp.  
Binning cutoff is 3 std above the mean read length.  
This seems to break when the number of reads in a bin is < 100.  



## Docker instructions
Make sure "Use Virtualization framework" is enabled in Docker Settings > General

Also make sure that Rosetta emulation is turned off in Docker Settings > Features in development

1. Build the docker image using:
```
# docker build -t {name of image} .
docker build -t nanopore .
```

2. Create and run the docker container using:
```
# docker run -it --rm --platform linux/amd64 --name {container name} -v {Path to your github repo that is being mounted}:/{where it's being mounted inside} {image name}
docker run -it --rm --platform linux/amd64 --name 241019_test -v /Users/joncchen/Documents/GitHub/Nanopore_plasmid_bin_assemble_polish:/app/ nanopore
```

3. Activate the conda environment inside the container
```
conda activate Nanopore
```

4. Run the command to initiate assembly + polishing.
```
# python nanopore_plasmid_bin_assemble_polish -i {input file location} -o {output folder} --canu_binary_path canu-2.2/bin/canu -t {# of CPU threads}
python nanopore_plasmid_bin_assemble_polish.py -i input/001_347A_reads.fastq.gz -o output --canu_binary_path canu-2.2/bin/canu -t 8
```
